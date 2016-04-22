"""
FBA test simulation with the model.
"""
from __future__ import print_function, division
import os
import cobra
import fba.cobra.cobra_tools as ct
from metabolism_settings import VERSION, RESULTS_DIR, DATA_DIR

##############################################################################
# Read FBC Model
##############################################################################
sbml = os.path.join(RESULTS_DIR, "Metabolism_matrices_{}_L3V1.xml".format(VERSION))
model = cobra.io.read_sbml_model(sbml)
# model = ct.read_sbml_fbc_model(sbml) # with additional FBC v1 information

print(len(model.reactions))  # 504
print(len(model.metabolites))  # 479  (336 + 104 -1) species + proteins - protein_species
print(len(model.genes))  # 142  (115)

# Count the genes
full_genes = []
for reaction in model.reactions:
    full_genes.extend([g.id for g in reaction.genes])
print(full_genes)
print(len(full_genes))
print(len(set(full_genes)))

# print objective coefficients
for key, value in ct.get_objective_coefficients(model).iteritems():
    print(key, value)

# mass balance & charge balance
# looks good for all reactions
# ix and ex are not balanced, also reactions involving proteins are not balanced
for r in model.reactions:
    mb = r.check_mass_balance()
    if mb:
        print(mb)
        print(r.reaction)

##############################################################################
# Simulate
##############################################################################
# The Model.optimize() function will return a Solution object, which will also be
# stored at model.solution. A solution object has several attributes:
#    f: the objective value
#    status: the status from the linear programming solver
#    x_dict: a dictionary of {reaction_id: flux_value} (also called "primal")
#    x: a list for x_dict
#    y_dict: a dictionary of {metabolite_id: dual_value}.
#    y: a list for y_dict

# Run optimization
print('#' * 80)
print('# Simulate')
print('#' * 80)

# Necessary to set upper bounds manually, because cobra is not working with +/-inf
maxreal = 1000
ct.set_max_bound(model, maxreal)
# ct.print_flux_bounds(model)

model.optimize()
model.solution.status
print('Solution status:', model.solution.status)
print('Solution:', model.solution)
model.summary()
print("biomassProduction", model.solution.x_dict['biomassProduction'])