'''
Reproducing the Matlab evolveState logic.

@author: Matthias Koenig
@date: 2015-03-12
'''

# TODO: read state


##############################################################################
# Initialize current state
##############################################################################
# additional index information for the updating of the state

# defined indexes (only created once)
# properties of the Metabolism process
substrateIndexs_externalExchangedMetabolites = state['substrateIndexs_externalExchangedMetabolites']  # [124x1]
substrateIndexs_internalExchangedLimitedMetabolites = state['substrateIndexs_internalExchangedLimitedMetabolites']  # [35x1]
substrateIndexs_limitableProteins = state['substrateIndexs_limitableProteins'] # [5x1]
substrateMonomerLocalIndexs = state['substrateMonomerLocalIndexs'] # [2x1]
substrateComplexLocalIndexs = state['substrateComplexLocalIndexs'] # [15x1]
proteinLimitableProteinComposition = state['proteinLimitableProteinComposition']  # [17x5]

# for the update has to be reshaped to fit to substrates
metabolismNewProduction = state['metabolismNewProduction'] # [585x3]


##############################################################################

# The properties substrates and enzymes represent the counts of metabolites
# and metabolic enzymes.

# TODO: store these values

# The substrate allocated for this time step
substrates = state['substrates']   # [585x3]

# Enzyme availability for time step
enzymes = state['enzymes']         # [104x1]

# TODO: missing (dryWieght ?)
# cellDryMass = sum(mass.cellDry);
cellDryMass = state['cellDryMass'] # [1x1]

##############################################################################
# Evolve state
##############################################################################
# import edu.stanford.covert.util.ConstantUtil;
# TODO: constants imported in original code, but unclear for what?

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
    print lowerBounds.shape
    print upperBounds.shape    
    lowerBounds
    
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


fluxBounds = calcFluxBounds(substrates, enzymes, fbaReactionBounds, fbaEnzymeBounds)

#-------------------------------------------------------------------------------

def calcGrowthRate(fluxBounds, fbaObj=fbaObjective, fbaSMat=fbaReactionStoichiometryMatrix):
    # import edu.stanford.covert.util.ComputationUtil;
    # import edu.stanford.covert.util.ConstantUtil;

    # flux bounds
    loFluxBounds = fluxBounds[:, 1];
    upFluxBounds = fluxBounds[:, 2];
            
    # real-valued linear programming
    loFluxBounds = max(loFluxBounds, -realmax);
    upFluxBounds = min(upFluxBounds,  realmax);

    # perform the FBA calculation    
    [fbaReactionFluxs, lambda, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                'maximize', fbaObj, fbaSMat, 
                fbaRightHandSide, loFluxBounds, upFluxBounds, ...
                'S', 'C', linearProgrammingOptions);
    if errFlag:
        warning('WholeCell:warning', 'Linear programming error: %s. Returning feasible, but possibly non-optimal solution x=0.', errMsg);
                fbaReactionFluxs = zeros(size(loFluxBounds));
            
    # growth
    fbaReactionFluxs = max(min(fbaReactionFluxs, upFluxBounds), loFluxBounds);
    growth = fbaReactionFluxs(fbaReactionIndexs_biomassProduction);
    reactionFluxs = zeros(size(this.reactionStoichiometryMatrix, 2), 1);
    reactionFluxs(this.reactionIndexs_fba) = fbaReactionFluxs(this.fbaReactionIndexs_metabolicConversion);
    
    # additional calculation of reducedCosts and Shadow prices
            if nargout > 3
                fbaReducedCosts = lambda.reducedCosts;
                reducedCosts = zeros(size(this.reactionStoichiometryMatrix, 2), 1);
                reducedCosts(this.reactionIndexs_fba) = fbaReducedCosts(this.fbaReactionIndexs_metabolicConversion);
                
                fbaShadowPrices = lambda.shadowPrices;
                shadowPrices = zeros(size(this.substrates));
                shadowPrices(this.substrateIndexs_fba) = fbaShadowPrices(this.fbaSubstrateIndexs_substrates);
            end
    return [growth, reactionFluxs, fbaReactionFluxs, ...
                reducedCosts, fbaReducedCosts, ...
                shadowPrices, fbaShadowPrices] 

def stochasticRound(value):
    '''
    Rounding of floats to integers weighted by decimal part
    '''
    roundUp = np.random.rand(value.size) < np.mod(value,1);
    value[roundUp] = np.ceil(value[roundUp]);
    value[~roundUp] = np.floor(value[~roundUp]);
    return value



def evolveState(substrates, enzymes):
    '''
    Evolves the given state.
    Check which values really are needed as input.
    '''
    # Calculate flux bounds
    fluxBounds = calcFluxBounds(substrates, enzymes, fbaReactionBounds, fbaEnzymeBounds)
    
    # Calculate growth rate
    [metabolicReaction.growth, metabolicReaction.fluxs, fbaReactionFluxs] = calcGrowthRate(fluxBounds);

    # Compute real-valued amount of nutrient imported, biomass (DNA,
    # RNA, protein, membrane, etc) precursors produced and consumed,
    # and monomers modified. Stochastically round to:
    #     1. Instantaneously, maintain integer-valued amounts of metabolites
    #     2. Over time, maintain experimentally determined ratios of biomass
    #        components (stored in this.metabolismProduction)

    # nutrient uptake
    substrates[substrateIndexs_externalExchangedMetabolites, compartmentIndexs_extracellular] = 
            substrates[substrateIndexs_externalExchangedMetabolites, compartmentIndexs_extracellular]
            - stochasticRound(fbaReactionFluxs[fbaReactionIndexs_metaboliteExternalExchange] * stepSizeSec);

    # recycled metabolites
    substrates[substrateIndexs_internalExchangedMetabolites] =
            substrates[substrateIndexs_internalExchangedMetabolites]
                + stochasticRound(fbaReactionFluxs[fbaReactionIndexs_metaboliteInternalExchange]);  
    
    # new metabolites
    substrates = substrates
                + stochasticRound(metabolismNewProduction * metabolicReaction.growth * stepSizeSec);
    
    # unaccounted energy consumption
    substrates[substrateIndexs_atpHydrolysis] = substrates[substrateIndexs_atpHydrolysis]
                + np.array([-1, -1, 1, 1, 1]) * stochasticRound(unaccountedEnergyConsumption * metabolicReaction.growth * stepSizeSec);
            
    ## make metabolites counts positive (in case stochastic rounding ##
    # made them slightly negative eg -1) except H2O, H+
    # TODO: find the ones smaller zero and set to zero
    # substrates[substrateMetaboliteLocalIndexs[:, 1], :] < 0 
            
    ## update cell mass & volume ##
    # TODO
    # mass_calcMass() - state cellMass
    # recalculate based on metabolite amount
    # geometry_calculateVolume() - state cellGeometry

    mass = CellMass() # get the state CellMass
    # mass.waterWt               %water weight (g)
    # mass.metaboliteWt          %metabolite weight (g)
    # mass.dnaWt                 %DNA weight (g)
    # mass.rnaWt                 %RNA weight (g)
    # mass.proteinWt             %protein weight (g)
    # mass.total                 %total weight of cell and media (g)
    # mass.cell                  %total weight of cell (g)
    # mass.cellDry               %total dry weight of cell (g)
    # mass.media                 %total weight of media (g)

