"""
Create a comp model.
Test script for working with the comp extension in SBML.

@author: Matthias Koenig
"""
from libsbml import *

from toymodel_settings import toy_comp_file
from toymodel_factory import *


def create_ExternalModelDefinition(mdoc, cid, sbml_file):
    extdef = mdoc.createExternalModelDefinition()

    extdef.setId(cid)
    extdef.setName(cid)
    extdef.setModelRef(cid)
    extdef.setSource(sbml_file)
    return extdef


def create_comp_model(sbml_file):

    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    mdoc = doc.getPlugin("comp")

    # create listOfExternalModelDefinitions
    emd_bounds = create_ExternalModelDefinition(mdoc, "toy_ode_bounds", sbml_file="toy_ode_bounds.xml")
    emd_update = create_ExternalModelDefinition(mdoc, "toy_ode_update", sbml_file="toy_ode_update.xml")
    emd_model = create_ExternalModelDefinition(mdoc, "toy_ode_model", sbml_file="toy_ode_model.xml")

    # create models and submodels
    model = doc.createModel()
    model.setId("toy_ode_comp")
    model.setName("Combined ODE/SSA model")
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


    # write SBML file
    write_and_check(doc, sbml_file)

    # flatten the model
    # Use the flattened model for simulation


if __name__ == "__main__":
    create_comp_model(toy_comp_file)
