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

createReaction(id="R1", name="R1", fast=False, reversible=True,
               reactants={"A": 1}, products={"B1": 1})
createReaction(id="R2", name="R2", fast=False, reversible=True,
               reactants={"B1": 1}, products={"B2": 1})
createReaction(id="R3", name="R3", fast=False, reversible=True,
               reactants={"B2": 1}, products={"C": 1})

# create objective function
objective = mplugin.createObjective()
objective.setId("R3_maximize")
objective.setType("maximize");
mplugin.setActiveObjectiveId("maximize");
fluxObjective = objective.createFluxObjective()
fluxObjective.setReaction("R3")
fluxObjective.setCoefficient(1.0)

"""
        # <fluxbounds>
        # parameters for dynamical calculation of flux bounds
        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds', 'lb_fbaEnzymeBounds', 'ub_fbaEnzymeBounds'):
            par = model.createParameter()
            par.setId('{}__{}'.format(index, p_name))
            par.setValue(r_fba_df[p_name][index])
            par.setConstant(True)
        # The reaction flux bounds are set as hard upper and lower flux bounds
        # These are NOT the dynamical flux bounds.

        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds'):
            # "lessEqual", "greaterEqual", "equal"
            bound = mplugin.createFluxBound();
            bound.setReaction(index);
            bound.setValue(r_fba_df[p_name][index])
            if p_name.startswith('lb'):
                bound.setId('lb__{}'.format(index))
                bound.setOperation("greaterEqual")
            if p_name.startswith('ub'):
                bound.setId('ub__{}'.format(index))
                bound.setOperation("lessEqual")

    # <objective function>

"""
# write sbml
fba_out = "fba_toy.xml"
writer = SBMLWriter()
writer.writeSBML(doc_fba, fba_out)

from sbml_tools.checks import check_sbml
check_sbml(fba_out)

#########################################################################################3
# simulate the FBA model

import cobra
cobra_model = cobra.io.read_sbml_model(fba_out)


### Running FBA
cobra_model.optimize()
print cobra_model.solution.status
# Output:
# 'optimal'
print cobra_model.solution.f
{reaction: reaction.objective_coefficient for reaction in cobra_model.reactions
 if reaction.objective_coefficient > 0}
