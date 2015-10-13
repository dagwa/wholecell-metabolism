"""
Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

@author: Matthias Koenig
"""

#################################
# create toy models
#################################
from toymodel.toymodel_settings import fba_file, ode_bounds_file, ode_update_file
from toymodel.toymodel_factory import create_fba, create_ode_bounds, create_ode_update
create_fba(fba_file)
create_ode_bounds(ode_bounds_file)
create_ode_update(ode_update_file)

#################################
# simulate kinetic flux bound model
#################################
# dynamically changing bounds
import roadrunner
print roadrunner.__version__
rr_bounds = roadrunner.RoadRunner(ode_bounds_file)
rr_bounds.reset()
# rr_bounds.selections = ['time', 'r1']
# changing parameters
rr_bounds.k1 = -0.2

result = rr_bounds.simulate(0, 10, steps=100)
rr_bounds.plot()
print result

#################################
# simulate simple FBA
#################################
# set bounds in FBA model from bounds calculation
rr_fba = roadrunner.RoadRunner(fba_file)
rr_fba.ub_R1 = rr_bounds.ub_R1

# sets the bounds in the cobra model
import cobra
from fba.cobra.cobra_tools import print_flux_bounds
cobra_fba = cobra.io.read_sbml_model(fba_file)
cobra_R1 = cobra_fba.reactions.get_by_id("R1")
cobra_R1.upper_bound = rr_bounds.ub_R1
print_flux_bounds(cobra_model)

# calculate FBA
cobra_model.optimize()
print cobra_model.solution.status
# Output:
# 'optimal'
print cobra_model.solution.f
{reaction: reaction.objective_coefficient for reaction in cobra_model.reactions
 if reaction.objective_coefficient > 0}

# set solution fluxes in rr_fba
# constant fluxes
for (rid, flux) in cobra_model.solution.x_dict.iteritems():
    pid = "v_{}".format(rid)
    rr_fba[pid] = flux

rr_fba.v_R1




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

