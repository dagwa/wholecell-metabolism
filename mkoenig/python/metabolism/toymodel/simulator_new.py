"""
Simulator to run the models consisting of multiple
modeling frameworks.

Simulation of the combined toy model consisting of FBA and kinetic submodels.
Using FBA simulator & kinetic simulator to simulate submodels with
synchronization between the partial simulations.

Simulating the model.

TODO: Fix the zero time point of the simulation, how to handle this correctly
TODO: refactor all the 1 time calculations out of the iteration (replacements)
TODO: write fba ssubmodel results in complete result vector (for plotting & visualization)
"""

from __future__ import print_function, division
import libsbml
import pandas as pd
import roadrunner
import cobra

from pandas import DataFrame
import sbmlutils.comp as comp
import numpy
import warnings

from collections import defaultdict

#################################################
MODEL_FRAMEWORK_FBA = 'fba'
MODEL_FRAMEWORK_ODE = 'ode'
MODEL_FRAMEWORK_STOCHASTIC = 'stochastic'
MODEL_FRAMEWORK_LOGICAL = 'logical'

MODEL_FRAMEWORKS = [
    MODEL_FRAMEWORK_FBA,
    MODEL_FRAMEWORK_ODE,
    MODEL_FRAMEWORK_STOCHASTIC,
    MODEL_FRAMEWORK_LOGICAL,
]
#################################################

class FBAModel(object):
    """ Handling FBA models

    """
    def __init__(self, fba_doc):
        self.fba_doc = fba_doc
        self.fba_model = fba_doc.getModel()

        # bounds are mappings from parameters to reactions
        #       parameter_id -> [rid1, rid2, ...]
        self.ub_parameters = defaultdict(list)
        self.lb_parameters = defaultdict(list)
        self.process_bounds()

    def process_bounds(self):
        """  Determine which parameters are upper or lower bounds.
        :return:
        :rtype:
        """
        for r in self.fba_model.getListOfReactions():
            mr = r.getPlugin("fbc")
            rid = r.getId()
            if mr.isSetUpperFluxBound():
                self.ub_parameters[mr.getUpperFluxBound()].append(rid)
            if mr.isSetLowerFluxBound():
                self.lb_parameters[mr.getLowerFluxBound()].append(rid)


class Simulator(object):
    """ Simulator class to simulate hybrid models.

    The simulator is initialized with the top level sbml file.
    """

    def __init__(self, top_level_file):
        """ Create the simulator with the top level SBML file.

        The models are resolved to their respective simulation framework.
        The top level network must be an ode network.
        """
        self.sbml_top = top_level_file
        # read top level model
        self.doc_top = libsbml.readSBMLFromFile(self.sbml_top)
        self.model_top = self.doc_top.getModel()
        self.framework_top = self.get_framework(self.model_top)
        if self.framework_top is not MODEL_FRAMEWORK_ODE:
            warnings.warn("The top level model framework is not ode: {}".format(self.framework_top))

        self.submodels = defaultdict(list)

        # memory for synchronization between models
        self.global_data = {}
        self._process_top_level()

    @staticmethod
    def get_framework(model):
        """ Get the framework for the given model/submodel object.
        Terms from the SBO modelling framework.

        This is the sbo which is set on the respective model/submodel element

        :param model:
        :return:
        """
        framework = None
        if model.isSetSBOTerm():
            sbo = model.getSBOTerm()
            if sbo == 624:
                framework = MODEL_FRAMEWORK_FBA
            elif sbo == 62:
                framework = MODEL_FRAMEWORK_ODE
            elif sbo == 63:
                framework = MODEL_FRAMEWORK_STOCHASTIC
            elif sbo in [234, 547]:
                framework = MODEL_FRAMEWORK_LOGICAL
            else:
                warnings.warn("Modelling Framework not supported: {}".format(sbo))
        return framework

    def __str__(self):
        """ Information string. """
        # top level
        print('-' * 80)
        print(self.doc_top, self.framework_top)
        print('-' * 80)
        # submodels
        for framework in MODEL_FRAMEWORKS:
            print("{:<10} : {}".format(framework, self.submodels[framework]))
        print('-' * 80)

    def _process_top_level(self):
        """ Process the top level information.

        Reads all the submodels, creates the global data structure.
        Order for executtion

        :return:
        """
        # get list of submodels
        model = self.doc_top.getModel()
        if model is None:
            warnings.warn("No top level model found.")

        # Get submodel frameworks & store in respective list
        top_plugin = self.model_top.getPlugin("comp")
        for submodel in top_plugin.getListOfSubmodels():
            # models are processed in the order they are listed in the listOfSubmodels
            framework = Simulator.get_framework(submodel)
            self.submodels[framework].append(submodel)

        print(self.__str__())


    def _prepare_models(self):
        """ Prepare the models for simulation.

        :return:
        :rtype:
        """
        pass
        '''
        # assign submodels to FBA
        fba_models = {}
        for key, value in model_frameworks.iteritems():
            if value["sbo"] == 624:
                print('FBA model')
                # get SBML file
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
                raise Exception("Modelling framework not supported, or not annotated on submodel: ", sbo)

        # submodels handled by FBA
        print(fba_models)

        # process FBA assignment rules
        # Necessary to find the assignment rules in the top model
        # These cannot be set in an hybrid approach.
        fba_rules = self._process_mixed_assignments(self.model_top)
        print('FBA rules:', fba_rules)

        # remove the FBA assignment rules from the model, so they can be set via the simulator
        for variable in fba_rules.values():
            print(variable)
            self.model_top.removeRuleByVariable(variable)

        import tempfile
        mixed_sbml_cleaned = tempfile.NamedTemporaryFile("w", suffix=".xml")
        libsbml.writeSBMLToFile(doc, mixed_sbml_cleaned.name)
        # mixed_sbml_cleaned.close()


        ###########################
        # Load ODE model
        ###########################
        # the roadrunner ode file is the flattend comp file.
        # FBA subparts do not change change any of the kinetic subparts (only connections via replaced bounds
        # and fluxes).
        # Consequently the ode part can be solved as is, only the iterative update between ode and fba has
        # to be performed
        rr_comp = roadrunner.RoadRunner(mixed_sbml_cleaned.name)
        sel = ['time'] \
              + ["".join(["[", item, "]"]) for item in rr_comp.model.getBoundarySpeciesIds()] \
              + ["".join(["[", item, "]"]) for item in rr_comp.model.getFloatingSpeciesIds()] \
              + rr_comp.model.getReactionIds() + fba_rules.values()
        rr_comp.timeCourseSelections = sel
        rr_comp.reset()

        ###########################
        # Load FBA models
        ###########################
        # via the submodels information it is decided which submodels are belonging to which
        # modeling framework (SBOTerms on submodel)
        doc = libsbml.readSBMLFromFile(mixed_sbml)
        model_frameworks = comp.get_submodel_frameworks(doc)
        model = doc.getModel()

        # get top level reaction ids
        # top_level_rids = frozenset([reaction.getId() for reaction in model.getListOfReactions()])

        # assign submodels to FBA
        fba_models = {}
        for key, value in model_frameworks.iteritems():
            if value["sbo"] == 624:
                print('FBA model')
                # get SBML file
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
                raise Exception("Modelling framework not supported, or not annotated on submodel: ", sbo)

        # submodels handled by FBA
        print(fba_models)
        '''

    def _process_mixed_assignments(model):
        """ Find the assignment rules which assign to a variable a reaction rate.

        This are the assignment rules synchronizing between FBA and ODE models.
        """
        fba_rules = {}

        for rule in model.getListOfRules():
            if not rule.isAssignment():
                continue
            variable = rule.getVariable()
            formula = rule.getFormula()
            parameter = model.getParameter(variable)
            if not parameter:
                continue
            reaction = model.getReaction(formula)
            if not reaction:
                continue
            fba_rules[reaction.getId()] = parameter.getId()
        return fba_rules

    @staticmethod
    def print_flux_bounds(model):
        """ Prints flux bounds for all reactions. """
        info = []
        for r in model.reactions:
            info.append([r.id, r.lower_bound, r.upper_bound])
        df = DataFrame(info, columns=['id', 'lb', 'ub'])
        pd.set_option('display.max_rows', len(df))
        print(df)
        pd.reset_option('display.max_rows')

    def simulate(self, tstart=0.0, tend=10.0, step_size=0.1, debug=False):
        """
        Performs model simulation.

        The simulator figures out based on the SBO terms in the list of submodels, which
        simulation/modelling framework to use.
        The passing of information between FBA and SSA/ODE is based on the list of replacements.
        """

        ###########################
        # Simulation
        ###########################
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
            # all fba submodels are simulated
            for fba_key, fba_info in fba_models.iteritems():
                cobra_model = fba_info['cobra']
                sbml_model = fba_info['doc'].getModel()


                # search in global parameter replacements for replacements which replace the bounds
                # of reactions
                # TODO: calculate the replacements only once
                print("* set bounds *")
                for p in model.getListOfParameters():
                    pid = p.getId()
                    mp = p.getPlugin("comp")
                    for rep_element in mp.getListOfReplacedElements():
                        # the submodel of the replacement belongs to the current model
                        if rep_element.getSubmodelRef() == fba_key:
                            # update upper and lower bounds
                            for rid in ub_parameters.get(pid, []):
                                print(rid, ': (upper) -> ', pid)
                                cobra_reaction = cobra_model.reactions.get_by_id(rid)
                                cobra_reaction.upper_bound = rr_comp[pid]
                            for rid in lb_parameters.get(pid, []):
                                print(rid, ': (lower) -> ', pid)
                                cobra_reaction = cobra_model.reactions.get_by_id(rid)
                                cobra_reaction.lower_bound = rr_comp[pid]

                # [2] optimize
                print("* optimize *")
                # TODO: start with last solution (speed improvement)
                cobra_model.optimize()

                # [3] set fluxes in ODE part
                print("* set fba fluxes in ode *")
                # based on replacements the fluxes are written in the kinetic part
                # set solution fluxes in parameters
                for (rid, flux) in cobra_model.solution.x_dict.iteritems():
                    if rid in fba_rules:
                        rr_comp[fba_rules[rid]] = flux
                        print(rid, ':', fba_rules[rid], "=", flux)
                    else:
                        print(rid, "no boundary flux")

                if debug:
                    print_flux_bounds(cobra_model)
                    print(cobra_model.solution.status)
                    print(cobra_model.solution.x_dict)
                    print("-" * 80)

            # --------------------------------------
            # ODE
            # --------------------------------------
            # With the updated fluxes from the FBA the kinetic part is
            # calculated:

            # simulate (1 step)
            if step_size:
                # constant step size
                result = rr_comp.simulate(0, end=step_size, steps=1)
            else:
                # variable step size
                result = rr_comp.simulate(0, steps=1, variableStep=True)

            # store results
            all_results.append(result[1])
            # set the fba fluxes in the model)
            # TODO: the fba fluxes are not set in the full kinetic result (shown as zero)
            #   these have to be set with the mapping between comp and flattened model
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
        '''

    def plot_species(self, df, path="species.png"):
        """ Plot species.

        :param df:
        :type df:
        :return:
        :rtype:
        """
        # create plots (use ids from flattened model for plotting)
        flat_doc = libsbml.readSBMLFromFile(flattened_file)
        flat_model = flat_doc.getModel()
        species_ids = ["[{}]".format(s.getId()) for s in flat_model.getListOfSpecies()]

        print(species_ids)

        ax_s = df.plot(x='time', y=species_ids)
        fig = ax_s.get_figure()
        fig.savefig(path)

    def plot_reactions(self, df, path="reactions.png"):
        """ Plot reactions.

        :param df:
        :type df:
        :return:
        :rtype:
        """
        flat_doc = libsbml.readSBMLFromFile(flattened_file)
        flat_model = flat_doc.getModel()
        reaction_ids = [r.getId() for r in flat_model.getListOfReactions()]

        # plot reactions
        print(reaction_ids)
        ax_r = df.plot(x='time', y=reaction_ids + ['vR3'])
        fig = ax_r.get_figure()
        fig.savefig(path)

    def save_csv(self, df, path="simulation.csv"):
        """ Save results to csv.
        """
        df.to_csv(path, sep="\t")


########################################################################################################################
if __name__ == "__main__":
    # Run simulation of the hybrid model
    from simsettings import top_level_file, flattened_file, out_dir
    import os
    os.chdir(out_dir)

    # Create simulator instance
    simulator = Simulator(top_level_file=top_level_file)
    df = simulator.simulate(tstart=0.0, tend=50.0, step_size=0.1, debug=True)
    simulator.plot_reactions(df)
    simulator.plot_species(df)
    simulator.save_csv(df)
