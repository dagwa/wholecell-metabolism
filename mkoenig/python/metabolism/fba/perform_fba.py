'''
Created on Mar 11, 2015

@author: mkoenig
'''
import os
from metabolism_settings import VERSION, RESULTS_DIR 
sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_{}_L3V1.xml".format(VERSION))

##############################################################################
# FBA with model
##############################################################################
import cobra

# Read model
# TODO: Encode the additional information used by COBRA in the model, i.e.
# things like Charge and composition (-> check mass balance)
model = cobra.io.read_sbml_model(sbml)
model.compartments
len(model.reactions)    # 1274  (504)
len(model.metabolites)  # 1779  (585*3=1755)
len(model.genes)        # 142   (104)

# Perform FBA
# linear programming
# linearProgrammingOptions = struct(...
#            'solver', 'glpk',...
#            'solverOptions', struct(...
#                'glpk', struct('lpsolver', 1, 'presol', 1, 'scale', 1, 'msglev', 0, 'tolbnd', 10e-7),...
#                'linprog', struct('Display','off'),...
#                'lp_solve', struct('verbose', 0, 'scaling', 3 + 64 + 128, 'presolve', 0), ...
#                'qsopt', struct()));
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

# TODO: necessary to syncronize the SBML model with the actual 
#       FBA problem which is performed. I.e. is necessary to get the 
#       ids of the objects.
#       I.e. the ids of substrates, reactions and genes.
# This can be done with the actual content of the matlab matrices.

'''
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
# constant problem data
##############################################################################
stepSizeSec = 1                         # defined in Process
realmax = 1e6                           # maximal flux bound
compartmentIndexs_cytosol       = 1;    # defined in Metabolism
compartmentIndexs_extracellular = 2;    # defined in Metabolism
compartmentIndexs_membrane      = 3;    # defined in Metabolism



# Load state data
import scipy.io
state = scipy.io.loadmat('this.mat')
for key, value in sorted(state.iteritems()):
    if isinstance(value, np.ndarray):
        print key, value.shape 
    else:
        print key
    
    
    
    

# TODO: get the ids from the names
reactionNames = state['reactionNames']  # [645, 1]
reactionNames.shape
print reactionNames
print reactionNames[0,0][0]

substrateNames = state['substrateNames']  # [585, 1]
substrateNames.shape 
# substrateNames

enzymeNames = state['enzymeNames']  # [104, 1]
enzymeNames.shape


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
# enzyme kinetics kcat (upper, lower)
fbaEnzymeBounds = state['fbaEnzymeBounds']     # [504x2]

# fbaObjective indicates which reaction represents the biomass production pseudoreaction. 
fbaObjective = state['fbaObjective'] # [504x1]
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


fbaReactionIndexs_metabolicConversion = state['fbaReactionIndexs_metabolicConversion'] # [336x1]
fbaReactionIndexs_metabolicConversion.shape

fbaReactionIndexs_metaboliteInternalExchange = state['fbaReactionIndexs_metaboliteInternalExchange'] # [42x1]
fbaReactionIndexs_metaboliteInternalExchange.shape

fbaReactionIndexs_biomassExchange = state['fbaReactionIndexs_biomassExchange'] # [1x1] (index 504)
fbaReactionIndexs_biomassExchange

fbaReactionIndexs_biomassProduction = state['fbaReactionIndexs_biomassProduction'] # [1x1] (index 503)
fbaReactionIndexs_biomassProduction
fbaReactionIndexs_metaboliteExternalExchange = state['fbaReactionIndexs_metaboliteExternalExchange']  # [124x1]
fbaReactionIndexs_metaboliteExternalExchange.shape
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
    N_Reactions = fbaReactionBounds.shape[0]            # 504
    lowerBounds =  -np.inf * np.ones([N_Reactions, 1])  # 504
    upperBounds =   np.inf * np.ones([N_Reactions, 1])  # 504
            
    # numbers of enzymes catalyzing each reaction, enzyme kinetics
    rxnEnzymes = fbaReactionCatalysisMatrix * enzymes;
    if applyEnzymeKineticBounds:
        lowerBounds = max(lowerBounds, fbaEnzymeBounds[:, 1] * rxnEnzymes);
        upperBounds = min(upperBounds, fbaEnzymeBounds[:, 2] * rxnEnzymes);
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
            lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 1]);
            upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 2]);
                
    # external metabolite availability
    if applyExternalMetaboliteBounds:
        index = fbaReactionIndexs_metaboliteExternalExchange
        upperBounds[index] = min(upperBounds[index],
                                    substrates[substrateIndexs_externalExchangedMetabolites, compartmentIndexs_extracellular]/stepSizeSec);
                
        lowerBounds[index] = max(lowerBounds[index], fbaReactionBounds[index, 1]*cellDryMass);
        upperBounds[index] = min(upperBounds[index], fbaReactionBounds[index, 2]*cellDryMass);

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

