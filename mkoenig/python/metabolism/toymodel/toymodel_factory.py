"""
Create the Toy sub models for coupling ODE models with FBA models.
Using hierarchical model composition to create the main model and the kinetic model
subparts.

@author: Matthias Koenig
"""
from libsbml import *

def create_compartment(model, cid, name, size=1.0, dims=3, constant=True):
    """ Create compartments. """
    c = model.createCompartment()
    c.setId(cid)
    c.setName(name)
    c.setSize(size)
    c.setSpatialDimensions(dims)
    c.setConstant(constant)
    return c

def create_species(model, sid, name, initialAmount, constant, boundaryCondition, compartment,
                  hasOnlySubstanceUnits=False):
    """ Create species. """
    s = model.createSpecies()
    s.setId(sid)
    s.setName(name)
    s.setInitialAmount(initialAmount)
    s.setConstant(constant)
    s.setBoundaryCondition(boundaryCondition)
    s.setCompartment(compartment)
    s.setHasOnlySubstanceUnits(hasOnlySubstanceUnits)
    return s

def create_parameter(model, pid, name, constant, value):
    """ Create paramters. """
    p = model.createParameter()
    p.setId(pid)
    p.setName(name)
    p.setConstant(constant)
    p.setValue(value)
    return p

def create_reaction(model, rid, name, fast=False, reversible=True,
                     reactants={}, products={}):
    """ Create basic reaction structure. """
    r = model.createReaction()
    r.setId(rid)
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


def set_flux_bounds(r, lb="lb", ub="ub"):
    """ Set flux bounds on given reaction. """
    rplugin = r.getPlugin("fbc")
    rplugin.setLowerFluxBound(lb)
    rplugin.setUpperFluxBound(ub)

def create_objective(mplugin, oid, otype, fluxObjectives, active=True):
    objective = mplugin.createObjective()
    objective.setId(oid)
    objective.setType(otype)
    if active:
        mplugin.setActiveObjectiveId("R3_maximize")
    for rid, coefficient in fluxObjectives.iteritems():
        fluxObjective = objective.createFluxObjective()
        fluxObjective.setReaction(rid)
        fluxObjective.setCoefficient(coefficient)
    return objective


def write_and_check(doc, sbml_file):
    # write and check the SBML file
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_file)
    from sbml_tools.checks import check_sbml
    check_sbml(sbml_file)

####################################################
# FBA submodel
####################################################
# FBA submodel in FBC v2 which uses parameters as flux bounds.
def create_fba(sbml_file):
    """
    Create the fba model.
    TODO: add units for species and fluxes, important for the correct update and synchronization.
    TODO: load into multiscale model data base
    TODO: SBO term for fluxBounds, species, ...
    """
    # Create the FBA submodel and run an optimization
    sbmlns = SBMLNamespaces(3, 1, "fbc", 2)
    doc_fba = SBMLDocument(sbmlns)
    doc_fba.setPackageRequired("fbc", False)
    model = doc_fba.createModel()
    mplugin = model.getPlugin("fbc")
    mplugin.setStrict(False)

    # model
    model.setId("fba_toy")
    model.setName("FBA submodel")

    # compartments
    c_ext = create_compartment(model, cid="ext", name="external compartment")
    c_int = create_compartment(model, cid="cell", name="cell")
    
    # boundary species
    s_A = create_species(model, sid="A", name="A", initialAmount=10, constant=False,
                        boundaryCondition=True, compartment=c_ext.getId())
    s_C = create_species(model, sid="C", name="C", initialAmount=0, constant=False,
                        boundaryCondition=True, compartment=c_ext.getId())
    # internal species
    s_B1 = create_species(model, sid="B1", name="B1", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_int.getId())
    s_B2 = create_species(model, sid="B2", name="B2", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_int.getId())
    s_D = create_species(model, sid="D", name="D", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())
    # parameters
    p_r1 = create_parameter(model, pid="r1", name="r1", constant=False, value=1.0)
    p_r3 = create_parameter(model, pid="r3", name="r3", constant=True, value=0.0)
    p_lb = create_parameter(model, pid="lb", name="lower bound", constant=True, value=0.0)
    p_ub = create_parameter(model, pid="ub", name="upper bound", constant=True, value=1000.0)

    # reactions
    r_R1 = create_reaction(model, rid="R1", name="A import (R1)", fast=False, reversible=True,
               reactants={"A": 1}, products={"B1": 1})
    r_R2 = create_reaction(model, rid="R2", name="B1 <-> B2 (R2)", fast=False, reversible=True,
               reactants={"B1": 1}, products={"B2": 1})
    r_R3 = create_reaction(model, rid="R3", name="B2 export (R3)", fast=False, reversible=True,
               reactants={"B2": 1}, products={"C": 1})

    # flux bounds
    set_flux_bounds(r_R1, lb="lb", ub="r1")
    set_flux_bounds(r_R2, lb="lb", ub="ub")
    set_flux_bounds(r_R3, lb="lb", ub="ub")

    # objective function
    create_objective(mplugin, oid="R3_maximize", otype="maximize", fluxObjectives={"R3": 1.0})

    # write SBML file
    write_and_check(doc_fba, sbml_file)

####################################################
# ODE submodels
####################################################

def create_assignment_rule(model, sid, formula):
    rule = model.createAssignmentRule()
    rule.setVariable(sid)
    astnode = parseL3FormulaWithModel(formula, model)
    if not astnode:
        print('Formula could not be parsed:', formula)
        print(getLastParseL3Error())
    rule.setMath(astnode)
    return rule

def create_rate_rule(model, sid, formula):
    rule = model.createRateRule()
    rule.setVariable(sid)
    astnode = parseL3FormulaWithModel(formula, model)
    if not astnode:
        print('Formula could not be parsed:', formula)
        print(getLastParseL3Error())
    rule.setMath(astnode)
    return rule



# [1] dynamical flux bounds calculation
def create_ode_bounds(sbml_file):
    """" Submodel for dynamically calculating the flux bounds. """
    sbmlns = SBMLNamespaces(3, 1)
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()

    # model
    model.setId("ode_bounds_toy")
    model.setName("ODE bound calculation submodel")

    # parameters to update
    p_r1 = create_parameter(model, pid="r1", name="r1", constant=False, value=1.0)
    p_k1 = create_parameter(model, pid="k1", name="k1", constant=False, value=-0.1)

    # assignment rules
    create_rate_rule(model, sid="r1", formula="k1*r1")

    # write SBML file
    write_and_check(doc, sbml_file)


if __name__ == "__main__":
    # write & check sbml
    from toymodel_settings import fba_file
    from toymodel_settings import ode_bounds_file
    create_fba(fba_file)
    create_ode_bounds(ode_bounds_file)