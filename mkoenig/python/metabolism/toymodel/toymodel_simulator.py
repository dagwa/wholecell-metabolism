"""
Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

@author: Matthias Koenig
"""

#################################
# create toy models
#################################
from toymodel.toymodel_settings import fba_file, ode_bounds_file, ode_update_file, ode_model_file
from toymodel.toymodel_factory import create_fba, create_ode_bounds, create_ode_update, create_ode_model
create_fba(fba_file)
create_ode_bounds(ode_bounds_file)
create_ode_update(ode_update_file, fba_file=fba_file)
create_ode_model(ode_model_file)

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
print_flux_bounds(cobra_fba)

# calculate FBA
cobra_fba.optimize()
print cobra_fba.solution.status
# Output:
# 'optimal'
print cobra_fba.solution.f
{reaction: reaction.objective_coefficient for reaction in cobra_fba.reactions
 if reaction.objective_coefficient > 0}

# set solution fluxes in rr_fba
# constant fluxes
for (rid, flux) in cobra_fba.solution.x_dict.iteritems():
    pid = "v_{}".format(rid)
    rr_fba[pid] = flux

print rr_fba.v_R1

#################################
# simulate metabolite update
#################################
rr_update = roadrunner.RoadRunner(ode_update_file)
# set the fluxes
for rid in cobra_fba.solution.x_dict:
    pid = "v_{}".format(rid)
    rr_update[pid] = rr_fba[pid]

# simulate
result = rr_update.simulate(0, 10, steps=100)
rr_update.plot()
print result

#################################
# simulate kinetic model
#################################
rr_model = roadrunner.RoadRunner(ode_model_file)
# synchronize concentrations
rr_model.model["init(C)"] = rr_update.C

result = rr_model.simulate(0, 10, steps=100)
rr_model.plot()
print result