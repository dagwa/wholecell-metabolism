'''
Reproducing the Matlab evolveState logic.

@author: Matthias Koenig
@date: 2015-03-26
'''
import numpy as np
from pandas import DataFrame
from libsbml import readSBML

class FluxBoundCalculator(object):
    
    def __init__(self, sbml, reaction_index, species_index, enzymes_index,
                        substrateIndexs, state):
        '''
            In the end the full calculation has to be possible solely 
            based on the given SBML.
        '''
        # necessary to have the species/reaction/enzymes to index mapping for lookup
        # and mapping between ids and indices
        # this will later be created from the sbml
        self.r_index = reaction_index
        self.s_index = species_index
        self.e_index = enzymes_index
        self.substrateIndexs = substrateIndexs
        self.s_acp_index = 570;         # handle Acetyl-Carrier Protein
        
        doc = readSBML(sbml)
        self.model = doc.getModel()
        self.protein_reactions = self.find_protein_reactions()
        
        # set reused attributes from state dict
        # the complete information should be parsed from the SBML
        # TODO: step by step reduce the state dependency.
        for key in [
                    # in SBML
                    'fbaReactionCatalysisMatrix',
                    'fbaReactionBounds',
                    'fbaEnzymeBounds',
                    # get the right reactions (TODO: code type in SBML)
                    'fbaReactionIndexs_metabolicConversion',
                    'fbaReactionIndexs_metaboliteInternalExchange',
                    'fbaReactionIndexs_biomassExchange',
                    'fbaReactionIndexs_biomassProduction',
                    
                    'fbaReactionIndexs_metaboliteExternalExchange',
                    'fbaReactionIndexs_metaboliteInternalLimitedExchange',
                    
                    'compartmentIndexs_cytosol',
                    'compartmentIndexs_extracellular',
                    
                    
                    # substrate for reactions (in SBML)
                    'substrateIndexs_externalExchangedMetabolites',
                    'substrateIndexs_internalExchangedLimitedMetabolites',
                    'substrateIndexs_limitableProteins',
                    
                    ]:
            setattr(self, key, state[key])
        
                
        # self.Nr = self.fbaReactionBounds.shape[0]  # 504 reactions (lower and upper for every necessary)
        self.Nr = len(self.model.getListOfReactions())
        # TODO: set as parameter and read
        self.stepSizeSec = 1
        
        
    def find_protein_reactions(self, base=True):
        '''
        All reactions which have at least one protein as substrate or product.
        Proteins are ProteinMonomers or ProteinComplexes and start with MG_ in the model.
        All enzymes are also species in the the sbml, so that the search has to be 
        performed based on the reactions.
        '''
        import re
        if base:
            # find the base proteins
            re_protein = "^MG_\d+_(.+?)MER"
        else:
            # find the base/modified proteins
            re_protein = "^MG_\d+_(.+?)"
        
        def add_proteins(list_of):
            for s in list_of:
                sid = s.getSpecies()
                m = re.search(re_protein, sid)
                if m:
                    proteins.add(m.group(0))
        
        protein_reactions = dict()
        for r in self.model.getListOfReactions():
            proteins = set()
            add_proteins(r.getListOfReactants())
            add_proteins(r.getListOfProducts())
            if len(proteins) > 0:
                protein_reactions[r.getId()] = proteins
        return protein_reactions
     
    def initBounds(self):
        # initialize bounds
        lowerBounds =  -np.inf * np.ones([self.Nr, 1])       # 504
        upperBounds =   np.inf * np.ones([self.Nr, 1])       # 504
        return (lowerBounds, upperBounds)
        
    
    def calcFluxBounds(self, substrates, enzymes, cellDryMass, state_fb=None,
                   applyEnzymeKineticBounds=True, applyEnzymeBounds=True, applyDirectionalityBounds=True,
                   applyExternalMetaboliteBounds=True, applyInternalMetaboliteBounds=True, applyProteinBounds=True):
        '''
        Compute reaction flux upper and lower bounds based on
        1. Enzyme kinetics         (fbaEnzymeBounds)
        2. Enzyme availability     (enzymes)
        3. Transport rates         (fbaReactionBounds)
        4. Metabolite availability (substrates)
        5. Protein availability    (enzymes)
        '''
        
        def print_difference(var, name):
            if state_fb:
                control = state_fb[name]
                print 'Equal {} : {}'.format(name, (var==control).all())
        
        (lowerBounds, upperBounds) = self.initBounds()
    
        # numbers of enzymes catalyzing each reaction, enzyme kinetics
        # every reaction has 1 or zero associated proteins (enzymes)
        # how many enzymes for every reaction
        rxnEnzymes = np.dot(self.fbaReactionCatalysisMatrix, enzymes);  # [504x104]*[104x1]=[504x1]
        print_difference(rxnEnzymes, 'rxnEnzymes')
        #okay
        
        # fix NAs in EnzymeBounds 
        self.fbaEnzymeBounds[np.isnan(self.fbaEnzymeBounds[:, 0]), 0] = -np.inf
        self.fbaEnzymeBounds[np.isnan(self.fbaEnzymeBounds[:, 1]), 1] = np.inf
        if state_fb:
            state_fb.fbaEnzymeBounds[np.isnan(state_fb.fbaEnzymeBounds[:, 0]), 0] = -np.inf
            state_fb.fbaEnzymeBounds[np.isnan(state_fb.fbaEnzymeBounds[:, 1]), 1] = np.inf
        print_difference(self.fbaEnzymeBounds, 'fbaEnzymeBounds')
        
        if applyEnzymeKineticBounds:
            for k in range(self.Nr):
                # BUG: This should not work in Matlab !
                # matlab: -inf * NaN = -inf !  (everybody else: -inf * NaN = NaN)
                # so we have to check for the NaNs resulting from [-inf*0] and [inf*0] (WTF matlab)
                lb = self.fbaEnzymeBounds[k, 0]*rxnEnzymes[k]
                if not np.isnan(lb):
                    lowerBounds[k] = max(lowerBounds[k], lb)
                ub = self.fbaEnzymeBounds[k, 1]*rxnEnzymes[k]
                if not np.isnan(ub):
                    upperBounds[k] = min(upperBounds[k], ub)
            print_difference(lowerBounds, "lowerBounds_applyEnzymeKineticBounds")
            print_difference(upperBounds, "upperBounds_applyEnzymeKineticBounds")
        
        # numbers of enzymes catalyzing each reaction, unkown enzyme kinetics
        if applyEnzymeBounds:
            # reaction depends on enzyme and enzyme is zero
            reaction_has_enzyme = np.sum(self.fbaReactionCatalysisMatrix, 1)  # any(fbaReactionCatalysisMatrix, 2), works due to exactly 1 protein per reaction
            for reaction_k in range(self.Nr):
                if (reaction_has_enzyme[reaction_k] and rxnEnzymes[reaction_k] <= 0):
                    lowerBounds[reaction_k] = 0
                    upperBounds[reaction_k] = 0
            print_difference(lowerBounds, "lowerBounds_applyEnzymeBounds")
            print_difference(upperBounds, "upperBounds_applyEnzymeBounds")
        
        # reaction directionality / thermodynamics (for subset of reactions)
        if applyDirectionalityBounds:
            # all indices with exception of metaboliteExternalExchange
            indices = [self.fbaReactionIndexs_metabolicConversion, 
                   self.fbaReactionIndexs_metaboliteInternalExchange,
                   self.fbaReactionIndexs_biomassExchange,
                   self.fbaReactionIndexs_biomassProduction]
 
            for index in indices:
                for k in index:
                    reaction_k = k-1 # zero indexing
                    lowerBounds[reaction_k] = max(lowerBounds[reaction_k], self.fbaReactionBounds[reaction_k, 0]);
                    upperBounds[reaction_k] = min(upperBounds[reaction_k], self.fbaReactionBounds[reaction_k, 1]);
            print_difference(lowerBounds, "lowerBounds_applyDirectionalityBounds")
            print_difference(upperBounds, "upperBounds_applyDirectionalityBounds")
        
        # external metabolite availability
        if applyExternalMetaboliteBounds:
            for k in range(len(self.fbaReactionIndexs_metaboliteExternalExchange)):
                reaction_k = self.fbaReactionIndexs_metaboliteExternalExchange[k]-1
                substrate_k = self.substrateIndexs_externalExchangedMetabolites[k]-1
                compartment_k = self.compartmentIndexs_extracellular-1
                
                # find the external metabolite concentration and limit via c/stepSizeSec
                upperBounds[reaction_k] = min(upperBounds[reaction_k],
                                                substrates[substrate_k, compartment_k]/self.stepSizeSec);
                
                lowerBounds[reaction_k] = max(lowerBounds[reaction_k], self.fbaReactionBounds[reaction_k, 0]*cellDryMass);
                upperBounds[reaction_k] = min(upperBounds[reaction_k], self.fbaReactionBounds[reaction_k, 1]*cellDryMass);
            print_difference(lowerBounds, "lowerBounds_applyExternalMetaboliteBounds")
            print_difference(upperBounds, "upperBounds_applyExternalMetaboliteBounds")
        
        # internal metabolite availability
        print_difference(substrates, "substrates")
        
        if applyInternalMetaboliteBounds:
            for (k, reaction_k) in enumerate(self.fbaReactionIndexs_metaboliteInternalLimitedExchange):
                reaction_k = reaction_k-1 # reactions are correct
                
                # substrate indexes are defined on the [585x3] substrate matrix :/
                row = self.substrateIndexs[self.substrateIndexs == self.substrateIndexs_internalExchangedLimitedMetabolites[k][0]]
                substrate_k = row[0]-1

                # lookup the FBA substrate index
                # print 'reaction_k, substrate_k:',reaction_k, substrate_k
                # print self.r_index.index[reaction_k], row.index
                
                # reshape to index via the linear index (! Check if reshaped via correct dimension)
                lin_substrates = np.reshape(substrates.transpose(), (substrates.size))
                lowerBounds[reaction_k] = max(lowerBounds[reaction_k], -lin_substrates[substrate_k]/self.stepSizeSec)

            print_difference(lowerBounds, "lowerBounds_applyInternalMetaboliteBounds")
            print_difference(upperBounds, "upperBounds_applyInternalMetaboliteBounds")
        
            
        # protein monomers and complexes
        if applyProteinBounds:
            # some species are proteins, these can limit the reactions
            
            # Matlab nightmare: indexing is not performed on the FBA matrices but full matrices in combination with permutation & logical indexing.
            # No possibility to rebuild this logic 1 to 1:
            
            # limitedReactions = any(any(...
            #        reactionStoichiometryMatrix([substrateMonomerLocalIndexs; substrateComplexLocalIndexs], reactionIndexs_fba, :) & ...
            #        ~permute(repmat(proteinLimitableProteinComposition * substrates(substrateIndexs_limitableProteins, :), [1 1 numel(reactionIndexs_fba)]), [1 3 2]), 3), 1);
            
            # [substrateMonomerLocalIndexs; substrateComplexLocalIndexs] [2; 15] = 17 (in full network), reduces to 12 in full set
            # proteinLimitableProteinComposition [17x5]
            
            
            # Implemented Interpretation: limited reactions are reactions which 
            #     [1] have at least one protein/modified protein as reactant or product
            #     [2] at least one respective unmodified (base) protein has a protein count of zero
            
            # WTF this part is just craziness
            '''
            Mapping of reactions to base proteins occuring as reactants/products
            'Aas4': {'MG_287_MONOMER'},
            'Aas5': {'MG_287_MONOMER'},
            'Aas7': {'MG_287_MONOMER'},
            'NrdEF_ADP': {'MG_229_231_TETRAMER'},
            'NrdEF_CDP': {'MG_229_231_TETRAMER'},
            'NrdEF_GDP': {'MG_229_231_TETRAMER'},
            'NrdEF_maintenance': {'MG_124_MONOMER', 'MG_229_231_TETRAMER'},
            'Ohr_H2O2': {'MG_454_DIMER'},
            'Ohr_rx': {'MG_124_MONOMER', 'MG_454_DIMER'},
            'OsmC_H2O2': {'MG_427_DIMER'},
            'OsmC_rx': {'MG_124_MONOMER', 'MG_427_DIMER'},
            'PlsC4': {'MG_287_MONOMER'},
            'PlsC5': {'MG_287_MONOMER'},
            'PlsC7': {'MG_287_MONOMER'},
            'PlsX4': {'MG_287_MONOMER'},
            'PlsX5': {'MG_287_MONOMER'},
            'PlsX7': {'MG_287_MONOMER'},
            'TrxB': {'MG_124_MONOMER'}}
            '''
            # In the end 5 proteins are limiting
            # thioredoxin - MG_124_MONOMER (FBA substrate) 
            # 'ribonucleoside-diphosphate reductase' - MG_229_231_TETRAMER (FBA substrate)
            # osmotically inducible peroxidase' - MG_427_DIMER (FBA substrate)
            # thiol-dependent peroxidase - MG_454_DIMER (FBA substrate
            # THIS CREATES PROBLEM: neither in enzymes, nor in substrates, i.e. not represented in FBA problem, 
            # but necessary for update
            # Special case has to be handeled
            # 'acyl carrier protein - MG_287_MONOMER (not FBA substrate/protein)

            
            for (reaction_id, proteins) in self.protein_reactions.iteritems():                
                # check if any of the base proteins in the reaction (reactant/product) is zero
                for protein_id in proteins:
                    if protein_id == 'MG_287_MONOMER':  # 'acyl carrier protein'
                        substrate_k = self.s_acp_index
                    else:
                        substrate_k = self.substrateIndexs[protein_id]
                    
                    # substrate indexes are defined on the [585x3] substrate matrix :/
                    # reshape to index via the linear index (! Check if reshaped via correct dimension)
                    lin_substrates = np.reshape(substrates.transpose(), (substrates.size))
                    N_protein = lin_substrates[substrate_k-1]
                            
                    if  N_protein <= 0:
                        # lookup the index
                        reaction_k = self.r_index.k[reaction_id]
                        lowerBounds[reaction_k] = 0;
                        upperBounds[reaction_k] = 0;
                        break; # break for this reaction
            print_difference(lowerBounds, "lowerBounds_applyProteinBounds")
            print_difference(upperBounds, "upperBounds_applyProteinBounds")       
        
        # return bounds
        bounds = np.concatenate((lowerBounds, upperBounds), axis=1)     # [504x2]
        bounds_df = DataFrame(bounds, index=self.r_index.index, columns=['lowerBounds', 'upperBounds'])
        return bounds_df


import fba.cobra.cobra_tools as ct
class GrowthCalculator(object):
    '''
        fluxBounds.shape
        Out[30]: (504, 2)

    '''
    realmax = 1e6
    
    def __init__(self, cobra_model, reaction_index):
        self.cobra_model = cobra_model
        self.r_index = reaction_index
        
    def calcGrowthRate(self, fluxBounds):
        ''' 
        Calculates the growth rate with given fluxBounds via FBA. 
        TODO: remove the copying of the fluxbounds.
        '''
        # flux bounds
        lowerBounds = fluxBounds.lowerBounds
        upperBounds = fluxBounds.upperBounds
            
        # real-valued linear programming with realmax bounds
        for k in range(len(lowerBounds)):
            lowerBounds[k] = max(lowerBounds[k], -self.realmax);
            upperBounds[k] = min(upperBounds[k],  self.realmax);
        
        # Set flux bounds in model
        bounds_df = DataFrame({'lowerBounds':lowerBounds, 'upperBounds':upperBounds}, index=self.r_index.index)
        ct.set_flux_bounds(self.cobra_model, bounds_df)

        # perform the FBA calculation
        # solution is created   
        self.cobra_model.optimize()
        self.solution = self.cobra_model.solution
        
        fbaReactionFluxs = self.solution.x 

        '''
    [fbaReactionFluxs, lambda, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, fbaSMat, 
                fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', linearProgrammingOptions);
    if errFlag:
        warning('WholeCell:warning', 'Linear programming error: %s. Returning feasible, but possibly non-optimal solution x=0.', errMsg);
                fbaReactionFluxs = zeros(size(loFluxBounds));
        '''
    
        # Solution has to be in flux bounds
        # ????? This should not be necessary if the solution was calculated properly ???
        # fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
        
    
        # Readout growth reaction
        growth = self.solution.x_dict['biomassProduction']
        
        # for idx in this.fbaReactionIndexs_metabolicConversion
        # reactionFluxs(this.reactionIndexs_fba) = fbaReactionFluxs(this.fbaReactionIndexs_metabolicConversion);
        reactionFluxs = np.zeros(len(self.r_index));
        idx = 0
        for reaction_id, row in self.r_index.iterrows():
            if row.type == "metabolicConversion":
                reactionFluxs[idx] = self.solution.x_dict[reaction_id]
            else:
                pass
            idx += 1
        
        return (growth, reactionFluxs, fbaReactionFluxs)
    