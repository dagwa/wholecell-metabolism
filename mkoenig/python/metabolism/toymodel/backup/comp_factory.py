"""
Create a comp model.
Test script for working with the comp extension in SBML.

One model composition combines all the kinetic models,
in addition the higher level comp model is created which combines everything (i.e. the FBA & ODE models).
For the simulation of the full combined model the tools have to figure out the subparts which are
simulated with which simulation environment.
"""
from __future__ import print_function
import warnings
from libsbml import *
from settings import *

import model_factory
import multiscale.sbmlutils.comp as comp
from multiscale.sbmlutils.factory import *
import multiscale.sbmlutils.sbmlio as sbml_io


def create_comp_ode_model(sbml_file):
    """
    Creates the ODE/SSA comp model.
    These are all the ode submodels combined without the FBA part.
    """
    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    mdoc = doc.getPlugin("comp")

    # create listOfExternalModelDefinitions
    # absolute file path for model resolving essential !
    emd_bounds = comp.create_ExternalModelDefinition(mdoc, "toy_ode_bounds", sbml_file=ode_bounds_file)
    emd_update = comp.create_ExternalModelDefinition(mdoc, "toy_ode_update", sbml_file=ode_update_file)
    emd_model = comp.create_ExternalModelDefinition(mdoc, "toy_ode_model", sbml_file=ode_model_file)

    # create models and submodels
    model = doc.createModel()
    model.setId("toy_comp_ode")
    model.setName("Combined ODE/SSA model")
    model_factory.add_generic_info(model)
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)
    mplugin = model.getPlugin("comp")

    # add listOfSubmodels which reference the External models
    comp.add_submodel_from_emd(mplugin, submodel_sid="bounds", emd=emd_bounds)
    comp.add_submodel_from_emd(mplugin, submodel_sid="update", emd=emd_update)
    comp.add_submodel_from_emd(mplugin, submodel_sid="model", emd=emd_model)

    # replace compartment
    create_compartments(model, [{
        A_ID: "extern", A_NAME: "external compartment", A_SPATIAL_DIMENSION: 3, A_VALUE: 1.0,
        A_UNIT: model_factory.UNIT_VOLUME
    }])
    comp.replace_compartment(model, 'extern', ['model', 'update'])

    # replace species
    create_species(model, [{
        A_ID: "C", A_NAME: "C", A_VALUE: 0.0, A_COMPARTMENT: "extern",
        A_CONSTANT: False, A_BOUNDARY_CONDITION: False, A_UNIT: model_factory.UNIT_AMOUNT,
        A_HAS_ONLY_SUBSTANCE_UNITS: True
    }])
    comp.replace_species(model, 'C', ['model', 'update'])

    # write SBML file
    sbml_io.write_and_check(doc, sbml_file)


def create_comp_full_model(sbml_file):
    """
    Creates the full comp model as combination of FBA and comp models.

    connections
    ------------
    [1] flux bounds
    kinetic reaction bounds => replace the FBA bounds
    comp_ode.submodel_bounds__ub_R1 => fba_R1.upper_bound

    [2] reaction rates (fluxes)
    FBA flux (all reactions) => replaces Reaction flux in kinetic model

    ** FBA **
    <listOfParameters>
        <parameter id="ub_R1" name="ub R1" value="1" units="item_per_s" constant="false"/>
        <parameter id="lb" name="lower bound" value="0" units="item_per_s" constant="true"/>
        <parameter id="ub" name="upper bound" value="1000" units="item_per_s" constant="true"/>
        <parameter id="v_R1" name="R1 flux" value="0" units="item_per_s" constant="false"/>
        <parameter id="v_R2" name="R2 flux" value="0" units="item_per_s" constant="false"/>
        <parameter id="v_R3" name="R3 flux" value="0" units="item_per_s" constant="false"/>
    </listOfParameters>
    """
    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    mdoc = doc.getPlugin("comp")

    # create listOfExternalModelDefinitions
    emd_bounds = comp.create_ExternalModelDefinition(mdoc, "toy_ode_bounds", sbml_file=ode_bounds_file)
    emd_update = comp.create_ExternalModelDefinition(mdoc, "toy_ode_update", sbml_file=ode_update_file)
    emd_model = comp.create_ExternalModelDefinition(mdoc, "toy_ode_model", sbml_file=ode_model_file)
    emd_fba = comp.create_ExternalModelDefinition(mdoc, "toy_fba", sbml_file=fba_file)

    # create models and submodels
    model = doc.createModel()
    model.setId("toy_comp_full")
    model.setName("Combined SSA with FBA model")
    model_factory.add_generic_info(model)
    mplugin = model.getPlugin("comp")
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)

    # add submodel which references the external model definition
    comp.add_submodel_from_emd(mplugin, submodel_sid="bounds", emd=emd_bounds)
    comp.add_submodel_from_emd(mplugin, submodel_sid="update", emd=emd_update)
    comp.add_submodel_from_emd(mplugin, submodel_sid="model", emd=emd_model)
    comp.add_submodel_from_emd(mplugin, submodel_sid="fba", emd=emd_fba)

    #############################
    # compartments
    #############################
    create_compartments(model, [
        {A_ID: "extern", A_NAME: "external compartment", A_VALUE: 1.0, A_SPATIAL_DIMENSION: 3, A_UNIT: model_factory.UNIT_VOLUME},
        {A_ID: 'cell', A_NAME: 'cell', A_VALUE: 1.0, A_SPATIAL_DIMENSION: 3, A_UNIT: model_factory.UNIT_VOLUME}
    ])
    # replaced compartments
    comp.replace_compartment(model, 'extern', ['model', 'update', 'fba'])
    comp.replace_compartment(model, 'cell', ['update', 'fba'])

    #############################
    # species
    #############################
    # replaced species
    # (fba species are not replaced, because they need their boundaryConditions for the FBA,
    #    and do not depend on the actual concentrations)
    create_species(model, [
        {A_ID: 'C', A_NAME: "C", A_VALUE: 0, A_UNIT: model_factory.UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "extern"},
    ])
    comp.replace_species(model, 'C', ['model', 'update'])

    #############################
    # parameters
    #############################
    create_parameters(model, [
        # bounds
        {A_ID: 'ub_R1', A_VALUE: 1.0, A_UNIT: model_factory.UNIT_FLUX,
         A_NAME: 'ub_r1', A_CONSTANT: False},
        # parameters (fluxes)
        {A_ID: "v_R1", A_NAME: "R1 flux", A_VALUE: 0.0, A_UNIT: model_factory.UNIT_FLUX, A_CONSTANT: False},
        {A_ID: "v_R2", A_NAME: "R2 flux", A_VALUE: 0.0, A_UNIT: model_factory.UNIT_FLUX, A_CONSTANT: False},
        {A_ID: "v_R3", A_NAME: "R3 flux", A_VALUE: 0.0, A_UNIT: model_factory.UNIT_FLUX, A_CONSTANT: False},
    ])
    # bounds -> fba
    comp.replace_parameters(model, 'ub_R1', ['fba', 'bounds'])
    # fba -> update
    comp.replace_parameters(model, 'v_R1', ['fba', 'update'])
    comp.replace_parameters(model, 'v_R2', ['fba', 'update'])
    comp.replace_parameters(model, 'v_R3', ['fba', 'update'])

    # write SBML file
    sbml_io.write_and_check(doc, sbml_file)


if __name__ == "__main__":
    create_comp_ode_model(comp_ode_file)
    create_comp_full_model(comp_full_file)

    # flatten the model
    # TODO: how to flatten with libsbml?
    # Use the flattened model for simulation
