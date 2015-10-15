"""
Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

Simulating the model.

@author: Matthias Koenig
"""
import roadrunner
import cobra
from fba.cobra.cobra_tools import print_flux_bounds

print roadrunner.__version__
print cobra.__version__

#################################
# load the models
#################################
from toymodel_settings import ode_bounds_file, ode_model_file, ode_update_file, fba_file

# roadrunner
rr_bounds = roadrunner.RoadRunner(ode_bounds_file)
rr_model = roadrunner.RoadRunner(ode_model_file)
rr_update = roadrunner.RoadRunner(ode_update_file)

# The complete kinetic part has to be simulated as one model.
rr_fba = roadrunner.RoadRunner(fba_file)
# cobra
cobra_fba = cobra.io.read_sbml_model(fba_file)

def simulate(time_start=0, time_end=10):
    # dynamically changing bounds
    result = rr_bounds.simulate(0, 10, steps=100)
    pass



rr_bounds.plot()
print result

#################################
# simulate simple FBA
#################################
# set new bounds in FBA model
rr_fba = roadrunner.RoadRunner(fba_file)
rr_fba.ub_R1 = rr_bounds.ub_R1

# sets the bounds in the cobra model

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

result = rr_model.simulate(0, 10, steps=10)
rr_model.plot()
print result