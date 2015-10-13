"""
Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

@author: Matthias Koenig
"""

#################################
# create toy models
#################################
from toymodel_settings import fba_file
from toymodel_factory import create_fba
create_fba(fba_file)



#################################
# simulate the kinetic model
#################################


#################################
# simulate simple FBA
#################################
import cobra
cobra_model = cobra.io.read_sbml_model(fba_file)

# [1] Simple FBA
cobra_model.optimize()
print cobra_model.solution.status
# Output:
# 'optimal'
print cobra_model.solution.f
{reaction: reaction.objective_coefficient for reaction in cobra_model.reactions
 if reaction.objective_coefficient > 0}


from fba.cobra.cobra_tools import print_flux_bounds
print_flux_bounds(cobra_model)


# [2] Dynamic FBA
# Change the boundaries dynamically, i.e. changing a parameter/concentration which is used
# in the calculation of the FBA boundaries.
# => recalculate the boundaries

