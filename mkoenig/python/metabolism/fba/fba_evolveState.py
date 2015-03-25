'''
Reproducing the Matlab evolveState logic.

@author: Matthias Koenig
@date: 2015-03-12
'''
import numpy as np


def calcFluxBounds(substrates, enzymes, fbaReactionBounds, fbaEnzymeBounds,
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

    # initialize
    Nr = fbaReactionBounds.shape[0]            # 504
    lowerBounds =  -np.inf * np.ones([Nr, 1])  # 504
    upperBounds =   np.inf * np.ones([Nr, 1])  # 504
    
    # iterate over the reactions
    for k, 
    # get the modifer protein:
    

    
    # numbers of enzymes catalyzing each reaction, enzyme kinetics
    
    rxnEnzymes = np.dot(fbaReactionCatalysisMatrix, enzymes);  # [504x104]*[104x1]=[504x1]
    print rxnEnzymes.shape
    print np.isnan(rxnEnzymes).any()
    if applyEnzymeKineticBounds:
        lowerBounds = np.maximum(lowerBounds, np.multiply(fbaEnzymeBounds[:, 0].reshape(Nr,1), rxnEnzymes));
        upperBounds = np.minimum(upperBounds, np.multiply(fbaEnzymeBounds[:, 1].reshape(Nr,1), rxnEnzymes));
    '''    
    # numbers of enzymes catalyzing each reaction, unkown enzyme kinetics
    if applyEnzymeBounds:
        # some logical indexing
        # TODO: understand & translate
        lowerBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
        upperBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
                
    # reaction directionality / thermodynamics
    if applyDirectionalityBounds:
        indices = [fbaReactionIndexs_metabolicConversion, 
                      fbaReactionIndexs_metaboliteInternalExchange,
                      fbaReactionIndexs_biomassExchange,
                      fbaReactionIndexs_biomassProduction]
        for index in indices:
            lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 0]);
            upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 1]);
                
    # external metabolite availability
    if applyExternalMetaboliteBounds:
        index = fbaReactionIndexs_metaboliteExternalExchange
        upperBounds[index] = min(upperBounds[index],
                                    substrates[substrateIndexs_externalExchangedMetabolites, compartmentIndexs_extracellular]/stepSizeSec);
                
        lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 0]*cellDryMass);
        upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 1]*cellDryMass);

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
    return bounds