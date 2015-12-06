"""
Simulator to run the models consisting of multiple
modeling frameworks.

Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

Simulating the model.

TODO: Create one comp file with all the ports and
    figure out the FBA and ODE subparts based on the annotations.
    The manual assignment which is done below has to be performed
    automatically based on the provided SBML information.
TODO: Fix the zero time point of the simulation
"""

from __future__ import print_function
import pandas as pd
from pandas import DataFrame
import roadrunner
import cobra
from settings import fba_file, comp_ode_file

print(roadrunner.__version__)
print(cobra.__version__)


def print_flux_bounds(model):
    """ Prints flux bounds for all reactions. """
    info = []
    for r in model.reactions:
        info.append([r.id, r.lower_bound, r.upper_bound])
    df = DataFrame(info, columns=['id', 'lb', 'ub'])
    pd.set_option('display.max_rows', len(df))
    print(df)
    pd.reset_option('display.max_rows')


#################################
# load ode and fba model
#################################
def simulate_manual(fba_sbml, comp_ode_sbml, tend=10.0, step_size=0.01, debug=True):
    """
    Performs the model integration.

    Manual connection of the fba model and the combined comp_ode model.
    The connections are hard coded. Proof-of-principle.

    :param tend: end time of the simulation
    :param step_size: step size for the integration, if None variable step size will be used
    :param debug: additional information
    :return: pandas solution data frame
    """

    # load fba model
    cobra_fba = cobra.io.read_sbml_model(fba_sbml)
    # ode model
    rr_comp = roadrunner.RoadRunner(comp_ode_sbml)
    sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getFloatingSpeciesIds()] \
        + rr_comp.model.getReactionIds()
    rr_comp.timeCourseSelections = sel
    rr_comp.reset()

    # store results
    all_results = [] 
    all_time = []
    result = None
    time = 0.0
    while time <= tend:
        if debug:
            print("-" * 80)
            print("Time: {}".format(time))
        
        # --------------------------------------
        # FBA
        # --------------------------------------
        # set bounds in cobra model
        cobra_R1 = cobra_fba.reactions.get_by_id("R1")
        cobra_R1.upper_bound = rr_comp.submodel_bounds__ub_R1
        
        # optimize
        cobra_fba.optimize()

        # set solution fluxes in rr_comp
        # constant fluxes
        for (rid, flux) in cobra_fba.solution.x_dict.iteritems():
            pid = "submodel_update__v_{}".format(rid)
            rr_comp[pid] = flux
            
        if debug:
            print_flux_bounds(cobra_fba)
            print(cobra_fba.solution.status)
            print(cobra_fba.solution.x_dict)
            print("-" * 80)
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
        
        # store results
        all_results.append(result[1])
        all_time.append(time)
        
        # store simulation values & get time step
        delta_time = result['time'][1]
        time = time + delta_time
        
        if debug:        
            print(result)

    # create result matrix    
    df = pd.DataFrame(data=all_results, columns=result.colnames)    
    df.time = all_time
    print(df)
    return df

def simulate(mixed_sbml, tend=10.0, step_size=0.01, debug=True):
    """
    Performs the model integration.

    The simulator has to figure out based on the replacement and SBO
    annotation which submodels to simulate with which algorithms.

    :param tend: end time of the simulation
    :param step_size: step size for the integration, if None variable step size will be used
    :param debug: additional information
    :return: pandas solution data frame
    """
    # TODO: implement
    raise Exception("NOT IMPLEMENTED")


if __name__ == "__main__":
    
    df1 = simulate_manual(fba_sbml=fba_file, comp_ode_sbml=comp_ode_file,
                          tend=50.0, step_size=0.1, debug=False)
    # df2 = simulate(tend=10.0, step_size=None, debug=False)

    df1.plot(x='time', y=['submodel_update__R1',
                          'submodel_update__R2',
                          'submodel_update__R3',
                          'submodel_model__R4'])
    df1.plot(x='time', y=['[submodel_update__A]',
                          '[submodel_update__B1]',
                          '[submodel_update__B2]',
                          '[C]',
                          '[submodel_model__D]'])

    # TODO: save figures as files

