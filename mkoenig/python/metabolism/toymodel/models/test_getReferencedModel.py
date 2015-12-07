from __future__ import print_function
from libsbml import *


def create_ExternalModelDefinition(mdoc, cid, sbml_file):
    extdef = mdoc.createExternalModelDefinition()
    extdef.setId(cid)
    extdef.setName(cid)
    extdef.setModelRef(cid)
    extdef.setSource(sbml_file)
    return extdef


sbmlns = SBMLNamespaces(3, 1, "comp", 1)
doc = SBMLDocument(sbmlns)
doc.setPackageRequired("comp", True)
mdoc = doc.getPlugin("comp")

# create listOfExternalModelDefinitions

from toymodel.settings import fba_file
print(fba_file)

emd = create_ExternalModelDefinition(mdoc, "toy_fba", sbml_file="toy_fba.xml")
model = emd.getReferencedModel()
print(model)

emd = create_ExternalModelDefinition(mdoc, "toy_fba", sbml_file=fba_file)
model = emd.getReferencedModel()
print(model)
