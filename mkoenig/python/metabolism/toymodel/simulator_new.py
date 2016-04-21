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

import logging
logging.basicConfig(format='%(message)s', level=logging.DEBUG)
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
    """ Handling FBA submodels models

    """
    # TODO: handle submodels directly defined in model

    def __init__(self, submodel, source, fba_rules):
        self.source = source
        self.submodel = submodel

        # read the model
        self.fba_doc = libsbml.readSBMLFromFile(source)
        self.fba_model = self.fba_doc.getModel()
        self.cobra_model = cobra.io.read_sbml_model(source)

        # parameters to replace in top model
        self.fba_rules = self.process_fba_rules(fba_rules)

        # bounds are mappings from parameters to reactions
        #       parameter_id -> [rid1, rid2, ...]
        self.ub_parameters = defaultdict(list)
        self.lb_parameters = defaultdict(list)
        self.ub_replacements = []
        self.lb_replacements = []
        self.process_bounds()

    def load_cobra_model(self):
        pass

    def process_fba_rules(self, fba_rules):
        """ Returns subset of fba_rules relevant for the FBA model.

        :param fba_rules:
        :type fba_rules:
        :return:
        :rtype:
        """
        rules = {}
        for rid, pid in fba_rules.iteritems():
            if self.fba_model.getReaction(rid) is not None:
                rules[rid] = pid
        return rules

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

    def process_replacements(self, top_model):
        """ Process the global replacements once. """
        for p in top_model.getListOfParameters():
            pid = p.getId()
            mp = p.getPlugin("comp")
            for rep_element in mp.getListOfReplacedElements():
                # the submodel of the replacement belongs to the current model
                if rep_element.getSubmodelRef() == self.submodel.getId():
                    # and parameter is part of the bounds
                    if pid in self.ub_parameters:
                        self.up_replacements.append(pid)
                    if pid in self.lb_parameters:
                        self.lb_replacements.append(pid)

    def update_fba_bounds(self, rr_comp):
        """
        Uses the global parameter replacements for replacements which replace the bounds
        of reactions.

        :param model:
        :type model:
        :return:
        :rtype:
        """
        logging.debug('* update_fba_bounds *')
        for pid in self.ub_replacements:
            for rid in self.ub_parameters.get(pid):
                logging.debug(rid, ': (upper) -> ', pid)
                cobra_reaction = self.cobra_model.reactions.get_by_id(rid)
                cobra_reaction.upper_bound = rr_comp[pid]

        for pid in self.lb_replacements:
            for rid in self.lb_parameters.get(pid):
                logging.debug(rid, ': (lower) -> ', pid)
                cobra_reaction = self.cobra_model.reactions.get_by_id(rid)
                cobra_reaction.lower_bound = rr_comp[pid]

    def optimize(self):
        """ Optimize the model """
        # TODO: start with last solution (speed improvement)
        logging.debug("* optimize *")
        self.cobra_model.optimize()

        # TODO: log
        self.log_flux_bounds()
        logging.debug('Solution status: {}'.format(self.cobra_model.solution.status))
        logging.debug('Solution fluxes: {}'.format(self.cobra_model.solution.x_dict))

    def set_ode_fluxes(self, rr_comp):
        """ Set fluxes in ODE part.

        Based on replacements the fluxes are written in the kinetic part
        :param rr_comp:
        :type rr_comp:
        :return:
        :rtype:
        """
        logging.debug("* set_ode_fluxes *")
        for rid, pid in self.fba_rules.iteritems():
            flux = self.cobra_model.solution.x_dict[rid]
            rr_comp[pid] = flux
            logging.debug('{}: {} = {}'.format(rid, pid, flux))

    def log_flux_bounds(self):
        """ Prints flux bounds for all reactions. """
        info = []
        for r in self.cobra_model.reactions:
            info.append([r.id, r.lower_bound, r.upper_bound])
        df = DataFrame(info, columns=['id', 'lb', 'ub'])
        pd.set_option('display.max_rows', len(df))
        logging.debug(df)
        pd.reset_option('display.max_rows')


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
        self.rr_comp = None
        self.fba_models = []

        self._process_top_level()
        self._prepare_models()

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

        Resolves the replacements and model couplings between the
        different frameworks and creates the respective simulatable
        models for the different frameworks.

        :return:
        :rtype:
        """
        logging.debug('_prepare_models')
        ###########################
        # find FBA rules
        ###########################
        # process FBA assignment rules of the top model
        self.fba_rules = self.find_fba_rules(self.model_top)
        logging.debug('FBA rules:', self.fba_rules)

        ###########################
        # prepare FBA models
        ###########################
        mdoc = self.doc_top.getPlugin("comp")
        for submodel in self.submodels[MODEL_FRAMEWORK_FBA]:
            mref = submodel.getModelRef()
            emd = mdoc.getExternalModelDefinition(mref)
            source = emd.getSource()
            fba_model = FBAModel(submodel=submodel, source=source, fba_rules=self.fba_rules)
            self.fba_models.append(fba_model)

        ###########################
        # prepare ODE model
        ###########################
        # the roadrunner ode file is the flattend comp file.
        # FBA subparts do not change change any of the kinetic subparts (only connections via replaced bounds
        # and fluxes).
        # Consequently the ode part can be solved as is, only the iterative update between ode and fba has
        # to be performed

        # remove FBA assignment rules from the model, so they can be set via the simulator
        for variable in self.fba_rules.values():
            self.model_top.removeRuleByVariable(variable)

        import tempfile
        mixed_sbml_cleaned = tempfile.NamedTemporaryFile("w", suffix=".xml")
        libsbml.writeSBMLToFile(self.doc_top, mixed_sbml_cleaned.name)

        rr_comp = roadrunner.RoadRunner(mixed_sbml_cleaned.name)
        sel = ['time'] \
              + ["".join(["[", item, "]"]) for item in rr_comp.model.getBoundarySpeciesIds()] \
              + ["".join(["[", item, "]"]) for item in rr_comp.model.getFloatingSpeciesIds()] \
              + rr_comp.model.getReactionIds() + self.fba_rules.values()
        rr_comp.timeCourseSelections = sel
        rr_comp.reset()
        self.rr_comp = rr_comp

    def find_fba_rules(self, top_model):
        """ Finds FBA rules in top model.

        Find the assignment rules which assign a reaction rate to a parameter.
        This are the assignment rules synchronizing between FBA and ODE models.

        These are Assignment rules of the form
            pid = rid
        i.e. a reaction rate is assigned to a parameter.
        """
        fba_rules = {}

        for rule in top_model.getListOfRules():
            if not rule.isAssignment():
                continue
            variable = rule.getVariable()
            formula = rule.getFormula()
            parameter = top_model.getParameter(variable)
            if not parameter:
                continue
            reaction = top_model.getReaction(formula)
            if not reaction:
                continue
            fba_rules[reaction.getId()] = parameter.getId()
        return fba_rules

    def simulate(self, tstart=0.0, tend=10.0, step_size=0.1):
        """
        Performs model simulation.

        The simulator figures out based on the SBO terms in the list of submodels, which
        simulation/modelling framework to use.
        The passing of information between FBA and SSA/ODE is based on the list of replacements.
        """

        logging.debug('###########################')
        logging.debug('# Simulation')
        logging.debug('###########################')
        all_results = []
        all_time = []
        result = None
        time = 0.0

        # variable step size integration
        if not step_size:
            self.rr_comp.integrator.setValue('variable_step_size', True)

        while time <= tend:
            logging.debug("-" * 80)
            logging.debug("Time: {}".format(time))

            # --------------------------------------
            # FBA
            # --------------------------------------
            for fba_model in self.fba_models:
                # update fba bounds from ode
                fba_model.update_fba_bounds(self.rr_comp)
                # optimize fba
                fba_model.optimize()
                # set ode fluxes from fba
                fba_model.set_ode_fluxes(self.rr_comp)

            # --------------------------------------
            # ODE
            # --------------------------------------
            if step_size:
                # constant step size
                logging.debug("rr_comp vR3 = {}".format(self.rr_comp['vR3']))
                result = self.rr_comp.simulate(start=0, end=step_size, steps=1)
            else:
                # variable step size
                result = self.rr_comp.simulate(start=0, steps=1)

            # store results
            # TODO: the fba fluxes are not set in the full kinetic result (shown as zero)
            #   these have to be set with the mapping between comp and flattened model
            # TODO: preinit array
            all_results.append(result[1])

            all_time.append(time)

            # store simulation values & get time step
            delta_time = result['time'][1]
            time = time + delta_time

            logging.debug(result)

        # create result matrix
        df = pd.DataFrame(data=all_results, columns=result.colnames)
        df.time = all_time
        return df

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
    logging.getLogger().setLevel(logging.DEBUG)
    from simsettings import top_level_file, flattened_file, out_dir
    import os
    os.chdir(out_dir)
    import timeit

    # Create simulator instance
    simulator = Simulator(top_level_file=top_level_file)

    start_time = timeit.default_timer()
    df = simulator.simulate(tstart=0.0, tend=50.0, step_size=0.1)
    elapsed = timeit.default_timer() - start_time
    logging.info("Simulation time: {}".format(elapsed))
    simulator.plot_reactions(df)
    simulator.plot_species(df)
    simulator.save_csv(df)
