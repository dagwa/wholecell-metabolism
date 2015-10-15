"""
Create a comp model.
Test script for working with the comp extension in SBML.

@author: Matthias Koenig
"""
from libsbml import *

from toymodel_settings import comp_file
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

    # shared compartment
    c_ext = create_compartment(model, cid="extern", name="external compartment")

    cplugin = c_ext.getPlugin("comp")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_model")
    replaced_element.setIdRef("extern")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_update")
    replaced_element.setIdRef("extern")

    # shared species
    s_C = create_species(model, sid="C", name="C", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())
    cplugin = s_C.getPlugin("comp")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_model")
    replaced_element.setIdRef("C")
    replaced_element = cplugin.createReplacedElement()
    replaced_element.setSubmodelRef("submodel_update")
    replaced_element.setIdRef("C")

    # write SBML file
    write_and_check(doc, sbml_file)

    # flatten the model
    # Use the flattened model for simulation


if __name__ == "__main__":
    create_comp_model(comp_file)
