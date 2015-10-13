"""
Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

@author: Matthias Koenig
"""

#################################
# create toy models
#################################
from toymodel_settings import fba_file, ode_bounds_file
from toymodel_factory import create_fba
create_fba(fba_file)

#################################
# simulate kinetic flux bound model
#################################
import roadrunner
print roadrunner.__version__
rr = roadrunner.RoadRunner(ode_bounds_file)
rr.selections = ['time', 'r1']
# changing parameters
rr.model.k1 = 0.2

result = rr.simulate(0, 10, steps=100)
rr.plot()
rr.reset()
print result

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

#################################
# simulate metabolite update
#################################
# TODO


#################################
# dynamic fba
#################################
# drive the FBA based on changed flux bounds and use
# results to update the species composition



# [2] Dynamic FBA
# Change the boundaries dynamically, i.e. changing a parameter/concentration which is used
# in the calculation of the FBA boundaries.
# => recalculate the boundaries

