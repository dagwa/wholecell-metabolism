"""
Create a comp model.
Test script for working with the comp extension in SBML.

@author: Matthias Koenig
"""
from libsbml import *

from toymodel_settings import toy_comp_file
from toymodel_factory import *

def create_comp_model(sbml_file):

    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", False)
    mdoc = doc.getPlugin("comp")

    model = doc.createModel()
    mplugin = model.getPlugin("comp")

    # create listOfExternalModelDefinitions
    extdef = CompSBMLDocumentPlugin.createExternalModelDefinition(mdoc)
    extdef.setId("toy_ode_bounds")
    extdef.setName("toy_ode_bounds")
    extdef.setModelRef("toy_ode_bounds")
    extdef.setSource("toy_ode_bounds.xml")


    # add listOfSubmodels which reference the External models

    # add list of Ports


    # model
    model = doc.createModel()
    model.setId("toy_ode_comp")
    model.setName("Combined ODE/SSA model")

    # write SBML file
    write_and_check(doc, sbml_file)

    # flatten the model
    # Use the flattened model for simulation


if __name__ == "__main__":
    create_comp_model(toy_comp_file)
