"""
Create a comp model.
Test script for working with the comp extension in SBML.

One model composition combines all the kinetic models,
in addition the higher level comp model is created which combines everything (i.e. the FBA & ODE models).
For the simulation of the full combined model the tools have to figure out the subparts which are
simulated with which simulation environment.
"""

from libsbml import *
from settings import comp_file

import model_factory
from multiscale.sbmlutils import comp
from multiscale.sbmlutils.factory import *
import multiscale.sbmlutils.io as sbml_io


def create_comp_model(sbml_file):
    """
    Creates the ODE/SSA comp model.
    """
    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    mdoc = doc.getPlugin("comp")

    # create listOfExternalModelDefinitions
    emd_bounds = comp.create_ExternalModelDefinition(mdoc, "toy_ode_bounds", sbml_file="toy_ode_bounds.xml")
    emd_update = comp.create_ExternalModelDefinition(mdoc, "toy_ode_update", sbml_file="toy_ode_update.xml")
    emd_model = comp.create_ExternalModelDefinition(mdoc, "toy_ode_model", sbml_file="toy_ode_model.xml")

    # create models and submodels
    model = doc.createModel()
    model.setId("toy_ode_comp")
    model.setName("Combined ODE/SSA model")
    model_factory.add_generic_info(model)
    mplugin = model.getPlugin("comp")

    # add listOfSubmodels which reference the External models
    submodel_bounds = mplugin.createSubmodel()
    submodel_bounds.setId("submodel_bounds")
    submodel_bounds.setModelRef(emd_bounds.getModelRef())

    submodel_update = mplugin.createSubmodel()
    submodel_update.setId("submodel_update")
    submodel_update.setModelRef(emd_update.getModelRef())

    submodel_model = mplugin.createSubmodel()
    submodel_model.setId("submodel_model")
    submodel_model.setModelRef(emd_model.getModelRef())

    # shared compartment
    create_compartments(model, [{
        A_ID: "extern", A_NAME: "external compartment", A_SPATIAL_DIMENSION: 3, A_VALUE: 1.0,
        A_UNIT: model_factory.UNIT_VOLUME
    }])

    c_ext = model.getCompartment("extern")
    cplugin = c_ext.getPlugin("comp")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_model")
    replaced_element.setIdRef("extern")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_update")
    replaced_element.setIdRef("extern")

    # shared species
    create_species(model, [{
        A_ID: "C", A_NAME: "C", A_VALUE: 0.0, A_COMPARTMENT: "extern",
        A_CONSTANT: False, A_BOUNDARY_CONDITION: False, A_UNIT: model_factory.UNIT_AMOUNT,
        A_HAS_ONLY_SUBSTANCE_UNITS: True
    }])

    s_C = model.getSpecies("C")
    cplugin = s_C.getPlugin("comp")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_model")
    replaced_element.setIdRef("C")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_update")
    replaced_element.setIdRef("C")

    # write SBML file
    sbml_io.write_and_check(doc, sbml_file)

    # flatten the model
    # TODO: how to flatten with libsbml?
    # Use the flattened model for simulation


if __name__ == "__main__":
    create_comp_model(comp_file)
