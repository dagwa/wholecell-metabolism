"""
Create the Toy sub models.
"""

from libsbml import *

# Create the FBA submodel and run an optimization
sbmlns = SBMLNamespaces(3, 1, "fbc", 1)
doc_fba = SBMLDocument(sbmlns)
doc_fba.setPackageRequired("fbc", False)
model = doc_fba.createModel()
mplugin = model.getPlugin("fbc")

# model name and id for process
model.setId("fba_toy")
model.setName("FBA submodel")


# create compartments
def createCompartment(id, name, size, dims, constant):
    """ Create compartments. """
    c = model.createCompartment()
    c.setId(id)
    c.setName(name)
    c.setSize(size)
    c.setSpatialDimensions(dims)
    c.setConstant(constant)
    return c

c_ext = createCompartment(id="ext", name="", size=1.0, dims=3, constant=True)
c_int = createCompartment(id="int", name="", size=1.0, dims=3, constant=True)

# create species
def createSpecies(id, name, initialAmount, constant, boundaryCondition, compartment,
                  hasOnlySubstanceUnits=False):
    s = model.createSpecies()
    s.setId(id)
    s.setName(name)
    s.setInitialAmount(initialAmount)
    s.setConstant(constant)
    s.setBoundaryCondition(boundaryCondition)
    s.setCompartment(compartment)
    s.setHasOnlySubstanceUnits(hasOnlySubstanceUnits)
    return s

s_A = createSpecies(id="A", name="A", initialAmount=10, constant=False,
                    boundaryCondition=False, compartment=c_ext.getId())
s_B1 = createSpecies(id="B1", name="B1", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_int.getId())
s_B2 = createSpecies(id="B2", name="B2", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_int.getId())
s_C = createSpecies(id="C", name="C", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_ext.getId())

# Create paramters

# create reactions

# create objective function

