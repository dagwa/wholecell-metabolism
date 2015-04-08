'''
Created on Mar 11, 2015

@author: mkoenig
'''
import os
import cobra
import fba.cobra.cobra_tools as ct
from metabolism_settings import VERSION, RESULTS_DIR, DATA_DIR
import scipy.io
import numpy as np
import pandas as pd
import fba.matlab.state_tools as state_tools

##############################################################################
# Read FBC Model
##############################################################################
# objective functions, gene associations, flux bounds

# sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_{}_L3V1.xml".format(VERSION))
sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_matrices_{}_L3V1.xml".format(VERSION))
# sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_4-l3-fbc.xml".format(VERSION))

reload(ct)
model = cobra.io.read_sbml_model(sbml)
model = ct.read_sbml_fbc_model(sbml) # with additional FBC v1 information


print len(model.reactions)    # 504
print len(model.metabolites)  # 479  (336 + 104 -1) species + proteins - protein_species
print len(model.genes)        # 142  (115)
print model.compartments

# Count the genes
full_genes = []
for reaction in model.reactions:
    full_genes.extend([g.id for g in reaction.genes])
print full_genes
print len(full_genes)
print len(set(full_genes))

# print objective coefficients
for key, value in ct.get_objective_coefficients(model).iteritems():
    print key, value

# mass balance & charge balance
# looks good for all reactions
# ix and ex are not balanced, also reactions involving proteins are not balanced
for r in model.reactions:
    mb = r.check_mass_balance()
    if mb:
        print mb
        print r.reaction
        
# Solution status
model.optimize()
model.solution.status
model.solution



# This cuts off some of the flux bounds
maxreal = 1E6
ct.set_max_bound(model, maxreal)
ct.print_flux_bounds(model)
model.optimize()
model.solution.status
sol = model.solution.x
# The Model.optimize() function will return a Solution object, which will also be stored at model.solution. A solution object has several attributes:
#    f: the objective value
#    status: the status from the linear programming solver
#    x_dict: a dictionary of {reaction_id: flux_value} (also called "primal")
#    x: a list for x_dict
#    y_dict: a dictionary of {metabolite_id: dual_value}.
#    y: a list for y_dict


# For lookup
# load matrices (most information can be read from the SBML)
# The goal is to fully encode the necessary information in the SBML.
matrix_dir = os.path.join(RESULTS_DIR, 'fba_matrices')

# handle the sodium NA id (not parsing as NaN)    
s_fba_df = pd.read_csv(os.path.join(matrix_dir, 's_fba.csv'), sep="\t", 
                           keep_default_na=False, na_values=('nan'))
s_fba_df.set_index('sid', inplace=True)

r_fba_df = pd.read_csv(os.path.join(matrix_dir, 'r_fba.csv'), sep="\t")
r_fba_df.set_index('rid', inplace=True)
r_fba_df

e_df = pd.read_csv(os.path.join(matrix_dir, 'e_fba.csv'), sep="\t")
e_df.set_index('eid', inplace=True)

##############################################################################
# Data for fluxbound and growth calculation
##############################################################################
# The properties substrates and enzymes represent the counts of metabolites
# and metabolic enzymes.
state_input = state_tools.read_state(os.path.join(DATA_DIR, 'matlab_dumps', 'input.mat'))
# this object from Process Metabolism
state = state_tools.read_state(os.path.join(DATA_DIR, 'matlab_dumps', 'this.mat'))
# state for calculating fluxbounds
state_fb = state_tools.read_state(os.path.join(DATA_DIR, 'matlab_dumps', 'fluxbounds.mat'))
# state for comparing growth calculation
state_output = state_tools.read_state(os.path.join(DATA_DIR, 'matlab_dumps', 'output.mat'))



##############################################################################
# Flux Bounds
##############################################################################
# init once with the necessary indexes
from pandas import DataFrame
import pandas as pd
Nr = len(r_fba_df)
Ns = len(s_fba_df)
Ne = len(e_df)
reaction_index = DataFrame({'k' : range(Nr), 'type': r_fba_df.type}, 
                            index=r_fba_df.index)
species_index = DataFrame(range(Ns), index=s_fba_df.index, columns=['k'])
enzymes_index = DataFrame(range(Ne), index=e_df.index, columns=['k'])
substrateIndexs = s_fba_df['substrateIndexs']

# Debug calculation against matlab dump
# Matlab dumps for flux

print 'Difference substrates:', np.sum(state_input.substrates - state_fb.substrates)        # check that same inputs
print 'Difference enzymes:', np.sum(state_input.enzymes - state_fb.enzymes)                 # check that same inputs
print 'Difference cellDryMass:', np.sum(state_input.cellDryMass - state_fb.cellDryMass)     # check that same inputs

# substrates (metabolite counts) allocated for time step
substrates = state_input.substrates   # [585x3]
# enzymes (protein counts) available for time step
enzymes = state_input.enzymes         # [104x1]
# cellDryMass
cellDryMass = state_input.cellDryMass

import fba.fba_evolveState as evolve 
reload(evolve)
fb_calc = evolve.FluxBoundCalculator(sbml, reaction_index, species_index, enzymes_index, 
                                     substrateIndexs, state)
fluxBounds = fb_calc.calcFluxBounds(substrates, enzymes, cellDryMass, state_fb=state_fb)


# calc growth rate
reload(ct)
reload(evolve)
growth_calc = evolve.GrowthCalculator(model, reaction_index)
(growth, reactionFluxs, fbaReactionFluxs) = growth_calc.calcGrowthRate(fluxBounds)

# Matlab growth output for comparison
growth - state_output.growth # growth is identical

for k in xrange(len(fbaReactionFluxs)):
    print reaction_index.ix[k].name,  fbaReactionFluxs[k] - state_output.fbaReactionFluxs[k]




# write the lp problem
from cobra.solvers import solver_dict
solver = 'cglpk'
problem = solver_dict[solver].create_problem(model)
problem.solve_problem()
problem.write('test.lp')





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

