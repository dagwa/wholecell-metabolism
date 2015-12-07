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
import libsbml
import pandas as pd
from pandas import DataFrame
import roadrunner
import cobra
import multiscale.sbmlutils.comp as comp

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



def simulate(mixed_sbml, tend=10.0, step_size=0.1, debug=False):
    """
    Performs the model integration.

    The simulator has to figure out based on the replacement and SBO
    annotation which submodels to simulate with which algorithms.

    Necessary to get the FBA submodels out and perform FBA on it, and clearly
    define the interfaces between the submodels.

    :param tend: end time of the simulation
    :param step_size: step size for the integration, if None variable step size will be used
    :param debug: additional information
    :return: pandas solution data frame
    """

    # ode model
    # the roadrunner ode file is the flattend comp file, the FBA part does not change any of the
    # kinetics, so it can just be integrated as is
    rr_comp = roadrunner.RoadRunner(mixed_sbml)
    sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr_comp.model.getFloatingSpeciesIds()] \
        + rr_comp.model.getReactionIds()
    rr_comp.timeCourseSelections = sel
    rr_comp.reset()

    # load fba sub models
    # get the dictionary of submodels and look based on the submodel SBOterms
    # which of them have to be simulated with FBA
    doc = libsbml.readSBMLFromFile(mixed_sbml)
    model_frameworks = comp.get_submodel_frameworks(doc)
    model = doc.getModel()

    # get the individual files of the submodels
    fba_models = {}
    for key, value in model_frameworks.iteritems():
        if value["sbo"] == 624:
            print('FBA model')
            # get the sbml file
            modelRef = value["modelRef"]
            mdoc = doc.getPlugin("comp")
            emd = mdoc.getExternalModelDefinition(modelRef)
            source = emd.getSource()
            print(source)
            fba_models[key] = {'cobra': cobra.io.read_sbml_model(source),
                               'doc': libsbml.readSBMLFromFile(source)}
        elif value['sbo'] == 62:
            print('ODE model')
            # get the sbml file
            modelRef = value["modelRef"]
            mdoc = doc.getPlugin("comp")
            emd = mdoc.getExternalModelDefinition(modelRef)
            source = emd.getSource()
            print(source)
        else:
            pass

    # this are the submodels handled via FBA:
    print(fba_models)

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
        # all fba submodels have to be simulated with FBA
        # and perform their updates
        for fba_key, fba_info in fba_models.iteritems():
            cobra_model = fba_info['cobra']
            sbml_model = fba_info['doc'].getModel()
            # [1] now the replacements have to be used to figure out what has to
            # be passed between the models.

            # <parameter id="ub_R1" name="ub_r1" value="1" units="item_per_s" constant="false">
            # <comp:listOfReplacedElements>
            #     <comp:replacedElement comp:idRef="ub_R1" comp:submodelRef="fba"/>
            #     <comp:replacedElement comp:idRef="ub_R1" comp:submodelRef="bounds"/>
            # </comp:listOfReplacedElements>
            # </parameter>

            # find all parameters replaced in the fba submodel

            # if any of them is either lower or upper bound of a reaction
            # get the reaction and set the upper bound in the model

            # the calculation of the replacement is only necessary once, than can be reused.
            # Done every time for laziness right now.

            # set bounds in fba model

            # <reaction id="R2" name="B1 &lt;-&gt; B2 (R2)" reversible="true" fast="false" fbc:lowerFluxBound="lb" fbc:upperFluxBound="ub">
            # rid: parameter

            # which parameters are upper or lower bounds
            from collections import defaultdict
            ub_parameters = defaultdict(list)
            lb_parameters = defaultdict(list)
            for r in sbml_model.getListOfReactions():
                mr = r.getPlugin("fbc")
                rid = r.getId()
                if mr.isSetUpperFluxBound():
                    ub_parameters[mr.getUpperFluxBound()].append(rid)
                if mr.isSetLowerFluxBound():
                    lb_parameters[mr.getLowerFluxBound()].append(rid)
            print(ub_parameters)
            print(lb_parameters)

            # search in global replacements
            for p in model.getListOfParameters():
                pid = p.getId()
                mp = p.getPlugin("comp")
                for rep_element in mp.getListOfReplacedElements():
                    # the submodel of the replacement belongs to the current model
                    if rep_element.getSubmodelRef() == fba_key:
                        # replace upper and lower bounds
                        for rid in ub_parameters.get(pid, []):
                            print(rid, ': (upper) -> ', pid)
                            cobra_reaction = cobra_model.reactions.get_by_id(rid)
                            cobra_reaction.upper_bound = rr_comp[pid]
                        for rid in lb_parameters.get(pid, []):
                            print(rid, ': (lower) -> ', pid)
                            cobra_reaction = cobra_model.reactions.get_by_id(rid)
                            cobra_reaction.lower_bound = rr_comp[pid]


            # R1 = cobra_model.reactions.get_by_id("R1")
            # R1.upper_bound = rr_comp['ub_R1']

            # [2] optimize
            cobra_model.optimize()

            # also via the replacement the things have to be written back in the ODE
            # subpart

            # set solution fluxes in rr_comp
            # constant fluxes
            for (rid, flux) in cobra_model.solution.x_dict.iteritems():
                pid = "v_{}".format(rid)
                rr_comp[pid] = flux

            if debug:
                print_flux_bounds(cobra_model)
                print(cobra_model.solution.status)
                print(cobra_model.solution.x_dict)
                print("-" * 80)

        # exit()

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


if __name__ == "__main__":

    # Run simulation of the hybrid model
    from settings import comp_full_file

    df = simulate(mixed_sbml=comp_full_file, tend=50.0, step_size=1)
    df.plot(x='time', y=['fba__R1', 'fba__R2', 'fba__R3', 'model__R4'])
    df.plot(x='time', y=['[update__A]',
                         '[update__B1]',
                         '[update__B2]',
                          '[C]',
                          '[model__D]'])

    """
    # Run iterative simulation of the two models
    if 0:
        from settings import fba_file, comp_ode_file
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

        # TODO: save figures as files and csv (results folder)
    """
