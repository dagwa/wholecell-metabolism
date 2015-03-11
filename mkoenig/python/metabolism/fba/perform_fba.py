'''
Created on Mar 11, 2015

@author: mkoenig
'''
import os
from metabolism_settings import VERSION, RESULTS_DIR 
sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_{}_L3V1.xml".format(VERSION))

#-------------------------------------
# FBA with model
#-------------------------------------
import cobra
# Read model
model = cobra.io.read_sbml_model(sbml)

# Perform FBA
model.optimize()

# Solution status
model.solution.status

# objective function
{reaction: reaction.objective_coefficient for reaction in model.reactions
if reaction.objective_coefficient != 0}

# Print flux bounds for all reactions
print '*** Reactions ***'
info = []
for r in model.reactions:
    info.append([r.id, r.upper_bound, r.lower_bound])
from pandas import DataFrame
df = DataFrame(info, columns=['id', 'lw', 'ub'])
print df
print '*'*80


# Calculate new flux bounds
# input copy numbers for the flux bounds
n_Metabolite = None
n_Protein = None
n_ProteinComplex = None

# Reaction information necessary to set flux bounds
# [reaction, n_Protein, n_ProteinComplex, kcat_f, kcat_b, kex_f, kex_b]

# Update metabolite copy numbers

##############################################################################
# Load data to evolve state
##############################################################################
import scipy.io
state = scipy.io.loadmat('/home/mkoenig/Desktop/matlab.mat')
print state.keys()


# Global state variables
cellDryMass = state['cellDryMass']

# The substrate allocated for this time step
substrates = state['substrates']   # [585, 3]


##############################################################################
# Evolve state
##############################################################################
# Dimemsions
N_reactions = 504
N_compartments = 3
N_metabolites = 585

# index in substrates
substrateIndexs_externalExchangedMetabolites = None # [?]
substrateIndexs_internalExchangedLimitedMetabolites
compartmentIndexs_extracellular


# Enzyme availability
enzymes = None      # [104x1]
# Transport rates
fbaReactionBounds = None
# Enzyme kinetics
fbaEnzymeBounds = None # [504, 2] (upper, lower)


substrateMonomerLocalIndexs = None
substrateComplexLocalIndexs = None

# Reaction Stoichiometry Matrxix [nMetabolites?, nReactions]
fbaReactionStoichiometryMatrix = None

fbaReactionCatalysisMatrix = None ?

# Index definition
fbaReactionIndexs_metabolicConversion
fbaReactionIndexs_metaboliteInternalExchange
fbaReactionIndexs_biomassExchange
fbaReactionIndexs_biomassProduction
fbaReactionIndexs_metaboliteExternalExchange

reactionIndexs_fba
proteinLimitableProteinComposition
substrateIndexs_limitableProteins

stepSizeSec = 1  # defined in Process

# 
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
    
    # import edu.stanford.covert.util.ConstantUtil;
    # TODO: create the constant file - What needed for ?
             
    # numbers
    nReactions = fbaReactionStoichiometryMatrix.shape[1]
            
    # initialize
    from numpy import inf
    lowerBounds =  -inf * ones([nReactions, 1])
    upperBounds =   inf * ones([nReactions, 1])
            
    # numbers of enzymes catalyzing each reaction, enzyme kinetics
    rxnEnzymes = fbaReactionCatalysisMatrix * enzymes;
    if applyEnzymeKineticBounds:
        lowerBounds = max(lowerBounds, fbaEnzymeBounds[:, 1] * rxnEnzymes);
        upperBounds = min(upperBounds, fbaEnzymeBounds[:, 2] * rxnEnzymes);
        
    # numbers of enzymes catalyzing each reaction, unkown enzyme kinetics
    if applyEnzymeBounds:
        # some logical indexing
        # TODO: understand
        lowerBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
        upperBounds(any(fbaReactionCatalysisMatrix, 2) & rxnEnzymes <= 0) = 0;
        
            
    # reaction directionality / thermodynamics
    if applyDirectionalityBounds:
        indices = [fbaReactionIndexs_metabolicConversion, 
                      fbaReactionIndexs_metaboliteInternalExchange,
                      fbaReactionIndexs_biomassExchange,
                      fbaReactionIndexs_biomassProduction]
        for index in indices:
            lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 1]);
            upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 2]);
                
        
    # external metabolite availability
    if applyExternalMetaboliteBounds:
        index = fbaReactionIndexs_metaboliteExternalExchange
        upperBounds[index] = min(upperBounds[index],
                                    substrates[substrateIndexs_externalExchangedMetabolites, compartmentIndexs_extracellular]/stepSizeSec);
                
        cellDryMass = sum(mass.cellDry);
        lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 1]*cellDryMass);
        upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 2]*cellDryMass);

    # internal metabolite availability
    if applyInternalMetaboliteBounds:
        index = fbaReactionIndexs_metaboliteInternalLimitedExchange
        lowerBounds[index] = max(lowerBounds[index], ...
                                - substrates[substrateIndexs_internalExchangedLimitedMetabolites]/stepSizeSec );
        
            
    # protein monomers and complexes
    if applyProteinBounds:
        # TODO: get indices for limited reactions
        limitedReactions = any(any(...
                    reactionStoichiometryMatrix([substrateMonomerLocalIndexs; substrateComplexLocalIndexs], reactionIndexs_fba, :) & ...
                    ~permute(repmat(proteinLimitableProteinComposition * substrates(substrateIndexs_limitableProteins, :), [1 1 numel(reactionIndexs_fba)]), [1 3 2]), 3), 1);
                    
        lowerBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
        upperBounds[fbaReactionIndexs_metabolicConversion[limitedReactions]] = 0;
            
            
    # TODO: create proper numpy array bound
    bounds = [lowerBounds upperBounds];
    return bounds


def calcGrowthRate(fluxBounds):
    pass


def evolveState():
    # Calculate flux bounds
    fluxBounds = calcFluxBounds(substrates, enzymes, fbaReactionBounds, fbaEnzymeBounds)
    
    # Calculate growth rate
    [metabolicReaction.growth, metabolicReaction.fluxs, fbaReactionFluxs] = calcGrowthRate(fluxBounds);

    ## nutrient uptake ##
    # TODO

    ## recycled metabolites ##
    # TODO    
        
    ## new metabolites ##
    # TODO     
    
    ## unaccounted energy consumption ##
    # TODO
            
    ## make metabolites counts positive (in case stochastic rounding ##
    # made them slightly negative eg -1) except H2O, H+
    # TODO
            
    ## update cell volume ##
    # TODO
    # mass_calcMass() - state cellMass
    # recalculate based on metabolite amount
    # geometry_calculateVolume() - state cellGeometry


