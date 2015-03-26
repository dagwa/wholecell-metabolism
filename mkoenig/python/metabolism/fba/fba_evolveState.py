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
            # TODO: calculate indices for limited reactions
            pass
            '''
            limitedReactions = any(any(...
                    reactionStoichiometryMatrix([substrateMonomerLocalIndexs; substrateComplexLocalIndexs], reactionIndexs_fba, :) & ...
                    ~permute(repmat(proteinLimitableProteinComposition * substrates(substrateIndexs_limitableProteins, :), [1 1 numel(reactionIndexs_fba)]), [1 3 2]), 3), 1);
                    
            lowerBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
            upperBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
            '''
        
        # return bounds
        bounds = np.concatenate((lowerBounds, upperBounds), axis=1) # [504x2]
        return bounds, rxnEnzymes
