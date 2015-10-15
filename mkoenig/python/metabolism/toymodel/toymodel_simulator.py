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
# load ode and fba model
#################################
from toymodel_settings import fba_file, comp_file

# fba model
cobra_fba = cobra.io.read_sbml_model(fba_file)
# ode model
rr_comp = roadrunner.RoadRunner(comp_file)
sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getFloatingSpeciesIds()] \
        + rr_comp.model.getReactionIds()
rr_comp.timeCourseSelections = sel

def simulate(tend=10, step_size=0.1):

    time = 0.0
    while time <= tend:
        print "-" * 80
        print "Time: {}".format(time)
        # --------------------------------------
        # FBA
        # --------------------------------------

        # set bounds in cobra model
        cobra_R1 = cobra_fba.reactions.get_by_id("R1")
        cobra_R1.upper_bound = rr_comp.submodel_bounds__ub_R1
        print_flux_bounds(cobra_fba)

        # optimize
        cobra_fba.optimize()
        print cobra_fba.solution.status
        print cobra_fba.solution.x_dict

        # set solution fluxes in rr_comp
        # constant fluxes
        for (rid, flux) in cobra_fba.solution.x_dict.iteritems():
            pid = "submodel_update__v_{}".format(rid)
            rr_comp[pid] = flux

        print "-" * 80
        # --------------------------------------
        # ODE
        # --------------------------------------
        # simulate (1 step)
        if step_size:
            # constant step size
            result = rr_comp.simulate(0, end=step_size, steps=1)
        else:
            # variable step size
            result = rr_comp.simulate(0, steps=1, variableStep=True)
        print result
        
        # store simulation values & get time step
        delta_time = step_size
        time = time + delta_time
        

if __name__ == "__main__":
    simulate(tend=1.0)
