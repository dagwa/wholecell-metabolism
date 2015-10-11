"""
Create the Toy sub models.
"""

from libsbml import *


# Create the FBA submodel and run an optimization
sbmlns = SBMLNamespaces(3, 1, "fbc", 2)
doc_fba = SBMLDocument(sbmlns)
doc_fba.setPackageRequired("fbc", False)
model = doc_fba.createModel()
mplugin = model.getPlugin("fbc")
mplugin.setStrict(True)

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
c_int = createCompartment(id="cell", name="", size=1.0, dims=3, constant=True)

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
                    boundaryCondition=True, compartment=c_ext.getId())
s_B1 = createSpecies(id="B1", name="B1", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_int.getId())
s_B2 = createSpecies(id="B2", name="B2", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_int.getId())
s_C = createSpecies(id="C", name="C", initialAmount=0, constant=False,
                    boundaryCondition=True, compartment=c_ext.getId())
s_D = createSpecies(id="D", name="D", initialAmount=0, constant=False,
                    boundaryCondition=False, compartment=c_ext.getId())

# Create parameters
def createParameter(id, name, constant, value):
    p = model.createParameter()
    p.setId(id)
    p.setName(name)
    p.setConstant(constant)
    p.setValue(value)
    return p

p_r1 = createParameter(id="r1", name="r1", constant=True, value=1.0)
p_r3 = createParameter(id="r3", name="0.0", constant=True, value=0.0)
p_lb = createParameter(id="lb", name="lower bound", constant=True, value=0.0)
p_ub = createParameter(id="ub", name="upper bound", constant=True, value=1000.0)

# create reactions with kinetic laws
def createReaction(id, name, fast=False, reversible=True,
                     reactants={}, products={}):
    r = model.createReaction()
    r.setId(id)
    r.setName(name)
    r.setFast(fast)
    r.setReversible(reversible)

    for sid, stoichiometry in reactants.iteritems():
        rt = r.createReactant()
        rt.setSpecies(sid)
        rt.setStoichiometry(abs(stoichiometry))
        rt.setConstant(True)

    for sid, stoichiometry in products.iteritems():
        rt = r.createProduct()
        rt.setSpecies(sid)
        rt.setStoichiometry(abs(stoichiometry))
        rt.setConstant(True)

    return r

r_R1 = createReaction(id="R1", name="R1", fast=False, reversible=True,
               reactants={"A": 1}, products={"B1": 1})
print r_R1
rplugin = r_R1.getPlugin("fbc")
rplugin.setLowerFluxBound("lb")
rplugin.setUpperFluxBound("ub")

r_R2 = createReaction(id="R2", name="R2", fast=False, reversible=True,
               reactants={"B1": 1}, products={"B2": 1})
r_R3 = createReaction(id="R3", name="R3", fast=False, reversible=True,
               reactants={"B2": 1}, products={"C": 1})
r_R4 = createReaction(id="R4", name="R4", fast=False, reversible=True,
               reactants={"C": 1}, products={"D": 1})

# create objective function
objective = mplugin.createObjective()
objective.setId("R3_maximize")
objective.setType("maximize")
mplugin.setActiveObjectiveId("R3_maximize")
fluxObjective = objective.createFluxObjectie()
fluxObjective.setReaction("R3")
fluxObjective.setCoefficient(1.0)

# create flux bounds
def createFluxBounds(r, lb="lb", ub="ub"):
    rplugin = r.getPlugin("fbc")
    rplugin.setLowerFluxBound(lb)
    rplugin.setUpperFluxBound(ub)

# special
createFluxBounds(r_R1, lb="lb", ub="r1")
createFluxBounds(r_R2, lb="lb", ub="ub")
createFluxBounds(r_R3, lb="lb", ub="ub")
createFluxBounds(r_R4, lb="lb", ub="ub")

# write sbml
fba_out = "fba_toy.xml"
writer = SBMLWriter()
writer.writeSBML(doc_fba, fba_out)

from sbml_tools.checks import check_sbml
check_sbml(fba_out)

"""
#########################################################################################3
# simulate the kinetic model




"""
#########################################################################################3
# simulate the FBA model

import cobra
cobra_model = cobra.io.read_sbml_model(fba_out)

# [1] Simple FBA
cobra_model.optimize()
print cobra_model.solution.status
# Output:
# 'optimal'
print cobra_model.solution.f
{reaction: reaction.objective_coefficient for reaction in cobra_model.reactions
 if reaction.objective_coefficient > 0}


from fba.cobra.cobra_tools import print_flux_bounds
print_flux_bounds(cobra_model)


# [2] Dynamic FBA
# Change the boundaries dynamically, i.e. changing a parameter/concentration which is used
# in the calculation of the FBA boundaries.
# => recalculate the boundaries

