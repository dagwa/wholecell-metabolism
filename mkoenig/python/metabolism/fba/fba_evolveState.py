'''
Reproducing the Matlab evolveState logic.

@author: Matthias Koenig
@date: 2015-03-12
'''
import numpy as np

class FluxBoundCalculator(object):
    
    def __init__(self, state):
        self.fbaReactionCatalysisMatrix = state['fbaReactionCatalysisMatrix']
        self.fbaReactionBounds = state['fbaReactionBounds']
        self.fbaEnzymeBounds = state['fbaEnzymeBounds']
        
        self.fbaReactionIndexs_metabolicConversion = state['fbaReactionIndexs_metabolicConversion'] 
        self.fbaReactionIndexs_metaboliteInternalExchange = state['fbaReactionIndexs_metaboliteInternalExchange']
        self.fbaReactionIndexs_biomassExchange = state['fbaReactionIndexs_biomassExchange']
        self.fbaReactionIndexs_biomassProduction = state['fbaReactionIndexs_biomassProduction'] 
        
        self.substrateIndexs_externalExchangedMetabolites
        
        self.Nr = self.fbaReactionBounds.shape[0]  # 504 reactions (lower and upper for every necessary)
        
    def calcFluxBounds(self, substrates, enzymes,
                   applyEnzymeKineticBounds=True, applyEnzymeBounds=True, applyDirectionalityBounds=True,
                   applyExternalMetaboliteBounds=True, applyInternalMetaboliteBounds=True, applyProteinBounds=True):
        '''
        Compute reaction flux upper and lower bounds based on
        1. Enzyme kinetics         (fbaEnzymeBounds)
        2. Enzyme availability     (enzymes)
        3. Transport rates         (fbaReactionBounds)
        4. Metabolite availability (substrates)
        5. Protein availability    (substrates)
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
        '''
        if applyEnzymeBounds:
            # some logical indexing
            # TODO: understand & translate
            lowerBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
            upperBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
        '''
                
        # reaction directionality / thermodynamics (for subset of reactions)
        if applyDirectionalityBounds:
            # all indices with exception of metaboliteExternalExchange
            indices = [self.fbaReactionIndexs_metabolicConversion, 
                   self.fbaReactionIndexs_metaboliteInternalExchange,
                   self.fbaReactionIndexs_biomassExchange,
                   self.fbaReactionIndexs_biomassProduction]
        
        for index in indices:
            for k in index:
                k = k-1 # zero indexing
                lowerBounds[k] = max(lowerBounds[k], self.fbaReactionBounds[k, 0]);
                upperBounds[k] = min(upperBounds[k], self.fbaReactionBounds[k, 1]);
        
        # external metabolite availability
        if applyExternalMetaboliteBounds:
            index = self.fbaReactionIndexs_metaboliteExternalExchange
            for k in index:
                k = k-1 # zero indexing
            # find the external metabolite concentration and limit via c/stepSizeSec
            upperBounds[k] = min(upperBounds[k],
                                    substrates[self.substrateIndexs_externalExchangedMetabolites[k], self.compartmentIndexs_extracellular]/self.stepSizeSec);
                
            lowerBounds[index] = max(lowerBounds[index], self.fbaReactionBounds[index, 0]*self.cellDryMass);
            upperBounds[index] = min(upperBounds[index], self.fbaReactionBounds[index, 1]*self.cellDryMass);
            
        '''
    # internal metabolite availability
    if applyInternalMetaboliteBounds:
        index = fbaReactionIndexs_metaboliteInternalLimitedExchange
        lowerBounds[index] = max(lowerBounds[index],
                                - substrates[substrateIndexs_internalExchangedLimitedMetabolites]/stepSizeSec );
        
    # protein monomers and complexes
    if applyProteinBounds:
        # TODO: calculate indices for limited reactions
        limitedReactions = any(any(...
                    reactionStoichiometryMatrix([substrateMonomerLocalIndexs; substrateComplexLocalIndexs], reactionIndexs_fba, :) & ...
                    ~permute(repmat(proteinLimitableProteinComposition * substrates(substrateIndexs_limitableProteins, :), [1 1 numel(reactionIndexs_fba)]), [1 3 2]), 3), 1);
                    
        lowerBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
        upperBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
        '''
        bounds = np.concatenate((lowerBounds, upperBounds), axis=1) # [504x2]
        return bounds, rxnEnzymes
