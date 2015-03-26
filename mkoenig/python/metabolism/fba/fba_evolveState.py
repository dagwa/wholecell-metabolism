'''
Reproducing the Matlab evolveState logic.

@author: Matthias Koenig
@date: 2015-03-26
'''
import numpy as np
from libsbml import readSBML

class FluxBoundCalculator(object):
    
    def __init__(self, sbml, state):
        '''
            In the end the full calculation has to be possible solely 
            based on the given SBML.
        '''
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
                    
                    ]:
            setattr(self, key, state[key])
                
        # self.Nr = self.fbaReactionBounds.shape[0]  # 504 reactions (lower and upper for every necessary)
        self.Nr = len(self.model.getListOfReactions())
        # TODO: set as parameter and read
        self.stepSizeSec = 1
        self.realmax = 1e6
        
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
    
        
    def calcFluxBounds(self, substrates, enzymes, cellDryMass,
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
        # initialize bounds
        lowerBounds =  -np.inf * np.ones([self.Nr, 1])       # 504
        upperBounds =   np.inf * np.ones([self.Nr, 1])       # 504
    
        # numbers of enzymes catalyzing each reaction, enzyme kinetics
        # every reaction has 1 or zero associated proteins (enzymes)
        # how many enzymes for every reaction
        rxnEnzymes = np.dot(self.fbaReactionCatalysisMatrix, enzymes);  # [504x104]*[104x1]=[504x1]
        if applyEnzymeKineticBounds:
            for k in range(self.Nr):
                lowerBounds[k] = max(lowerBounds[k], self.fbaEnzymeBounds[k, 0]*rxnEnzymes[k])
                upperBounds[k] = min(upperBounds[k], self.fbaEnzymeBounds[k, 1]*rxnEnzymes[k])
        
        # numbers of enzymes catalyzing each reaction, unkown enzyme kinetics
        if applyEnzymeBounds:
            # reaction depends on enzyme and enzyme is zero
            reaction_has_enzyme = np.sum(self.fbaReactionCatalysisMatrix, 1)  # any(fbaReactionCatalysisMatrix, 2), works due to exactly 1 protein per reaction
            for reaction_k in range(self.Nr):
                if (reaction_has_enzyme[reaction_k] and rxnEnzymes[reaction_k] <= 0):
                    lowerBounds[reaction_k] = 0
                    upperBounds[reaction_k] = 0
        
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
            
        
        # internal metabolite availability
        if applyInternalMetaboliteBounds:
            for k in self.fbaReactionIndexs_metaboliteInternalLimitedExchange:
                reaction_k = k-1
                substrate_k = self.substrateIndexs_internalExchangedLimitedMetabolites[k]-1
                compartment_k = self.compartmentIndexs_cytosol
                lowerBounds[reaction_k] = max(lowerBounds[reaction_k],
                                - substrates[substrate_k, compartment_k]/self.stepSizeSec);
        
        # protein monomers and complexes
        if applyProteinBounds:
            # some species are proteins, these can limit the reactions
            # complex
            
            # TODO: calculate indices for limited reactions
            # Matlab nightmare: indexing not done on FBA matrices, in combination with permutation & logical indexing
            # limitedReactions = any(any(...
            #        reactionStoichiometryMatrix([substrateMonomerLocalIndexs; substrateComplexLocalIndexs], reactionIndexs_fba, :) & ...
            #        ~permute(repmat(proteinLimitableProteinComposition * substrates(substrateIndexs_limitableProteins, :), [1 1 numel(reactionIndexs_fba)]), [1 3 2]), 3), 1);
            
            # [substrateMonomerLocalIndexs; substrateComplexLocalIndexs] [2; 15] = 17 (in full network), reduces to 12 in full set
            # for instance MG_454_DIMER_ox ?
            # proteinLimitableProteinComposition [17x5]
            
            
            # Interpretation: limited reactions are reactions which 
            # [1] have at least one protein/modified protein as reactant or product
            # [2] the respective unmodified (base) protein has a protein count of zero
            # So necessary to find the protein reactions & the respective base protein for the modified proteins
            '''
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
            for reaction_id in self.protein_reactions.keys():
                # TODO: necessary to put everything from indices to ids
                # Some things are easier with the ids, some with the indices
                reaction_k = None
            
            pass
            # TODO      
            # lowerBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
            # upperBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
            
        
        # return bounds
        bounds = np.concatenate((lowerBounds, upperBounds), axis=1) # [504x2]
        return bounds, rxnEnzymes
