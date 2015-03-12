'''
Creating the SBML based on the actual FBA performed in the Karr model
based on the read stoichiometric matrix. 

@author: Matthias Koenig
@date: 2015-03-12


   Knowledge Base
   ===============
   M. genitalium metabolism was reconstructed from a variety of sources,
   including flux-balance analysis metabolic models of other bacterial species
   and the reaction kinetics database SABIO-RK, and was organized into 504
   reactions in the knowledge base. These reactions are loaded into this process
   by the initializeConstants method.

      Object       No.
      ==========   ===
      substrates   585
      enzymes      104
      reactions    504
        chemical   ?
        transport  ?

   Representation
   ===============
   The properties substrates and enzymes represent the counts of metabolites
   and metabolic enzymes.

   fbaReactionStoichiometryMatrix represents the stoichiometry and compartments
   of metabolites and biomass in each of the 641 chemical/transport reactions,
   exchange pseudoreactions, and biomass production pseudoreaction.

   fbaReactionCatalysisMatrix represents the enzyme which catalyzes each
   reaction. 
   fbaEnzymeBounds represents the foward and backward kcat of the
   catalyzing enzyme of each reaction. 
   fbaReactionBounds represents the maximal import and export rates of each
   metabolite. 
   fbaObjective indicates which reaction represents the biomass production pseudoreaction. 
   fbaRightHandSide is a vector of zeros representing the change 
   in concentration over time of each metabolite and biomass. 
   metabolismProduction is redundant with the
   biomass production reaction in fbaReactionStoichiometryMatrix.
   metabolismProduction is calculated by summing the metabolic demands of all
   the other processes over the entire cell cycle. The table below lists the
   units of several properties of this process.

%      Property                       Units
%      ===========================    ==============================
%      fbaEnzymeBounds                molecules/enzyme/s
%      fbaReactionBounds              molecules/(gram dry biomass)/s
%      substrates                     molecules
%      enzymes                        molecules
%      stepSizeSec                    s
%      lowerBounds                    reactions/s
%      upperBounds                    reactions/s
%      growth                         cell/s
%      biomassComposition             molecules/cell
%      metabolismProduction           molecules/cell
%      chamberVolume                  L
%      setValues                      molecules/chamber
%      growthAssociatedMaintanence
'''

##############################################################################
# constant definitions in Karr model
##############################################################################
stepSizeSec = 1                         # defined in Process
realmax = 1e6                           # maximal flux bound
compartmentIndexs_cytosol       = 1;    # defined in Metabolism
compartmentIndexs_extracellular = 2;    # defined in Metabolism
compartmentIndexs_membrane      = 3;    # defined in Metabolism

# Load state data from matlab dump
import scipy.io
import numpy as np
state = scipy.io.loadmat('this.mat')
for key, value in sorted(state.iteritems()):
    if isinstance(value, np.ndarray):
        print key, value.shape 
    else:
        print key

###############################################
# FULL SET
###############################################
# 645 Reactions
# 585 Substrates
# 104 Enzymes
# 3 Compartments
# reactionStoichiometryMatrix (585, 645, 3)
###############################################
rids = state["reactionWholeCellModelIDs"]
sids = state["substrateWholeCellModelIDs"]
eids = state["enzymeWholeCellModelIDs"]

def cleanWholeCellIds(ids, col_id):
    df = DataFrame(ids, columns=[col_id])
    for k in df.index:
        df.ix[k] = df.ix[k][0]
    df = df.set_index(df[col_id])
    return df

r_df = cleanWholeCellIds(rids, 'rid')
s_df = cleanWholeCellIds(sids, 'sid')
e_df = cleanWholeCellIds(eids, 'eid')    
print r_df
print s_df
print e_df

###############################################
# FBA SET
###############################################
# 504 Reactions
# 376 Substrates
# 104 Enzymes

# It is clear where the different components are stored, but it is not clear
# which reactions and substrates correspond to the fba substrates and reactions.

# <Reactions>
# 336 Conversion
# 124 ExternalExchange
# 42 InternalExchange
# 1 biomassProduction
# 1 biomassExchange
# sum = 504
fbaReactionIndexs_metabolicConversion = state['fbaReactionIndexs_metabolicConversion'] # [336x1]
fbaReactionIndexs_metaboliteExternalExchange = state['fbaReactionIndexs_metaboliteExternalExchange']  # [124x1]
fbaReactionIndexs_metaboliteInternalExchange = state['fbaReactionIndexs_metaboliteInternalExchange'] # [42x1]
fbaReactionIndexs_biomassProduction = state['fbaReactionIndexs_biomassProduction'] # [1x1] (index 503)
fbaReactionIndexs_biomassExchange = state['fbaReactionIndexs_biomassExchange'] # [1x1] (index 504)

reactionIndexs_fba = state["reactionIndexs_fba"] # [336,1]

# Reaction identifier
r_fba_df = DataFrame(range(0,504), columns=['rid'])
# metabolic conversion
for k, item in enumerate(reactionIndexs_fba):
    index = item[0]-1
    print index
    r_fba_df.loc[k] = r_df['rid'][index]
r_fba_df

# metabolite external Exchange
for k, index in enumerate(fbaReactionIndexs_metaboliteExternalExchange):
    index = index-1
    r_fba_df.loc[index] = 'metaboliteExternalExchange_ex_{}'.format(k)

# metabolite internal Exchange
for k, index in enumerate(fbaReactionIndexs_metaboliteInternalExchange):
    index = index-1
    r_fba_df.loc[index] = 'metaboliteInternalExchange_ix_{}'.format(k)

# biomass production
index = fbaReactionIndexs_biomassProduction[0]-1
r_fba_df.loc[index] = 'biomassProduction'.format(k)

# biomass consumption
index = fbaReactionIndexs_biomassExchange[0]-1
r_fba_df.loc[index] = 'biomassConsumption'.format(k)


# <Substrates>
# 368 fba substrate
# 7 internal exchange
# 1 biomass
# sum = 376
fbaSubstrateIndexs_substrates = state['fbaSubstrateIndexs_substrates']# [368x1]
fbaSubstrateIndexs_metaboliteInternalExchangeConstraints = state['fbaSubstrateIndexs_metaboliteInternalExchangeConstraints'] # [7x1]
fbaSubstrateIndexs_biomass = state['fbaSubstrateIndexs_biomass'] # [1x1]
substrateIndexs_fba = state["substrateIndexs_fba"] # [376,1]

# substrate identifier
s_fba_df = DataFrame(range(0,376), columns=['sid'])
s_fba_df
# substrates
for k, item in enumerate(substrateIndexs_fba):
    index = item[0]-1
    compartment = 'c'
    if (index >= 585):
        compartment = 'e'
    if (index >= 2*585):
        compartment = 'm'
    # get the columns via modulo
    index = index % 585
    
    # add the correct compartment    
    sid ='{}_{}'.format(s_df['sid'][index], compartment)
    
    if compartment != 'c':
        k = k + 8    
    s_fba_df.loc[k] = sid
    
# internal exchange
for k, index in enumerate(fbaSubstrateIndexs_metaboliteInternalExchangeConstraints):
    index = index-1
    s_fba_df.loc[index] = 'metaboliteInternalExchangeConstraints_ix_{}'.format(k)
s_fba_df
# biomass
for k, index in enumerate(fbaSubstrateIndexs_biomass):
    index = index-1
    s_fba_df.loc[index] = 'biomass'.format(k)

s_fba_df.set_index(s_fba_df['sid'])
s_fba_df


# TODO: annotate the matrixes with the actual names

# Reaction Stoichiometry Matrxix [nMetabolites?, nReactions]
# fbaReactionStoichiometryMatrix represents the stoichiometry and compartments
#   of metabolites and biomass in each of the 641 chemical/transport reactions,
#   exchange pseudoreactions, and biomass production pseudoreaction.
fbaReactionStoichiometryMatrix = state['fbaReactionStoichiometryMatrix']  # [376x504]
fbaReactionStoichiometryMatrix.shape

# fbaReactionCatalysisMatrix represents the enzyme which catalyzes each
# reaction. 
fbaReactionCatalysisMatrix = state['fbaReactionCatalysisMatrix'] # [504x104]
fbaReactionCatalysisMatrix.shape


# maximal import and export rates (upper, lower)
fbaReactionBounds = state['fbaReactionBounds'] # [504x2]
np.isnan(fbaReactionBounds).any()

# enzyme kinetics kcat (upper, lower)
fbaEnzymeBounds = state['fbaEnzymeBounds']     # [504x2]
np.isnan(fbaEnzymeBounds).any()
fbaEnzymeBounds[np.isnan(fbaEnzymeBounds[:,0]),0] = -np.inf
fbaEnzymeBounds[np.isnan(fbaEnzymeBounds[:,1]),1] = np.inf
np.isnan(fbaEnzymeBounds).any()
print fbaEnzymeBounds

# fbaObjective indicates which reaction represents the biomass production pseudoreaction. 
fbaObjective = state['fbaObjective']            # [504x1]
fbaObjective.shape

# fbaRightHandSide is a vector of zeros representing the change 
# in concentration over time of each metabolite and biomass. 
fbaRightHandSide = state['fbaRightHandSide'] # [376x1]
fbaRightHandSide.shape


# defined indexes (only created once)
# properties of the Metabolism process
substrateIndexs_externalExchangedMetabolites = state['substrateIndexs_externalExchangedMetabolites']  # [124x1]
substrateIndexs_externalExchangedMetabolites.shape
substrateIndexs_internalExchangedLimitedMetabolites = state['substrateIndexs_internalExchangedLimitedMetabolites']  # [35x1]
substrateIndexs_internalExchangedLimitedMetabolites.shape
substrateIndexs_limitableProteins = state['substrateIndexs_limitableProteins'] # [5x1]
substrateIndexs_limitableProteins.shape

substrateMonomerLocalIndexs = state['substrateMonomerLocalIndexs'] # [2x1]
substrateMonomerLocalIndexs.shape

substrateComplexLocalIndexs = state['substrateComplexLocalIndexs'] # [15x1]
substrateComplexLocalIndexs.shape

reactionIndexs_fba = state['reactionIndexs_fba'] # [336x1]
reactionIndexs_fba.shape 

proteinLimitableProteinComposition = state['proteinLimitableProteinComposition']  # [17x5]
proteinLimitableProteinComposition.shape

metabolismNewProduction = state['metabolismNewProduction'] # [585x3]
metabolismNewProduction.shape
print metabolismNewProduction[:,0]




##############################################################################
# Initialize current state
##############################################################################







##############################################################################
# Initialize current state
##############################################################################
# The properties substrates and enzymes represent the counts of metabolites
# and metabolic enzymes.

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

