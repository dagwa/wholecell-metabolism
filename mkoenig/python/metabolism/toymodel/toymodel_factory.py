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
                     reactants={}, products={}, formula=None):
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

    if formula:
        # set formula in reaction
        astnode = parseL3FormulaWithModel(formula, model)
        if not astnode:
            print('Formula could not be parsed:', formula)
            print(getLastParseL3Error())
        law = r.createKineticLaw()
        law.setMath(astnode)

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

def write_and_check(doc, sbml_file):
    # write and check the SBML file
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_file)
    from sbml_tools.checks import check_sbml
    check_sbml(sbml_file)

####################################################
# ODE flux bounds
####################################################
# ODE model for dynamical flux bound calculation
def create_ode_bounds(sbml_file):
    """" Submodel for dynamically calculating the flux bounds. """
    sbmlns = SBMLNamespaces(3, 1)
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()

    # model
    model.setId("toy_ode_bounds")
    model.setName("ODE bound calculation submodel")

    # parameters to update
    create_parameter(model, pid="ub_R1", name="ub r1", constant=False, value=1.0)

    # assignment rules
    create_parameter(model, pid="k1", name="k1", constant=False, value=-0.2)
    create_rate_rule(model, sid="ub_R1", formula="k1*ub_R1")

    # write SBML file
    write_and_check(doc, sbml_file)

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
    model.setId("toy_fba")
    model.setName("FBA submodel")

    # compartments
    c_ext = create_compartment(model, cid="extern", name="external compartment")
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
    # parameters (bounds)
    create_parameter(model, pid="ub_R1", name="ub R1", constant=False, value=1.0)
    create_parameter(model, pid="lb", name="lower bound", constant=True, value=0.0)
    create_parameter(model, pid="ub", name="upper bound", constant=True, value=1000.0)
    # parameters (fluxes)
    create_parameter(model, pid="v_R1", name="R1 flux", constant=False, value=0.0)
    create_parameter(model, pid="v_R2", name="R2 flux", constant=False, value=0.0)
    create_parameter(model, pid="v_R3", name="R3 flux", constant=False, value=0.0)

    # reactions with constant flux
    r_R1 = create_reaction(model, rid="R1", name="A import (R1)", fast=False, reversible=True,
                           reactants={"A": 1}, products={"B1": 1},
                           formula="v_R1")
    r_R2 = create_reaction(model, rid="R2", name="B1 <-> B2 (R2)", fast=False, reversible=True,
                           reactants={"B1": 1}, products={"B2": 1},
                           formula="v_R2")
    r_R3 = create_reaction(model, rid="R3", name="B2 export (R3)", fast=False, reversible=True,
                           reactants={"B2": 1}, products={"C": 1},
                           formula="v_R3")

    # flux bounds
    set_flux_bounds(r_R1, lb="lb", ub="ub_R1")
    set_flux_bounds(r_R2, lb="lb", ub="ub")
    set_flux_bounds(r_R3, lb="lb", ub="ub")

    # objective function
    create_objective(mplugin, oid="R3_maximize", otype="maximize", fluxObjectives={"R3": 1.0})

    # write SBML file
    write_and_check(doc_fba, sbml_file)

####################################################
# ODE species update
####################################################
# model for update of species count


def create_ode_update(sbml_file, fba_file):
    """ Submodel for dynamically updating the metabolite count.
        Very similar model to the FBA model.
    """
    # Create the update model without fbc content
    sbmlns = SBMLNamespaces(3, 1)
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()

    # model
    model.setId("toy_ode_update")
    model.setName("ODE metabolite update submodel")

    # compartments
    c_ext = create_compartment(model, cid="extern", name="external compartment")
    c_int = create_compartment(model, cid="cell", name="cell")

    # boundary species released
    s_A = create_species(model, sid="A", name="A", initialAmount=10, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())
    s_C = create_species(model, sid="C", name="C", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())

    # internal species
    s_B1 = create_species(model, sid="B1", name="B1", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_int.getId())
    s_B2 = create_species(model, sid="B2", name="B2", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_int.getId())

    # parameters (fluxes)
    create_parameter(model, pid="v_R1", name="R1 flux", constant=False, value=0.0)
    create_parameter(model, pid="v_R2", name="R2 flux", constant=False, value=0.0)
    create_parameter(model, pid="v_R3", name="R3 flux", constant=False, value=0.0)

    # reactions with constant flux
    r_R1 = create_reaction(model, rid="R1", name="A import (R1)", fast=False, reversible=True,
                           reactants={"A": 1}, products={"B1": 1},
                           formula="v_R1")
    r_R2 = create_reaction(model, rid="R2", name="B1 <-> B2 (R2)", fast=False, reversible=True,
                           reactants={"B1": 1}, products={"B2": 1},
                           formula="v_R2")
    r_R3 = create_reaction(model, rid="R3", name="B2 export (R3)", fast=False, reversible=True,
                           reactants={"B2": 1}, products={"C": 1},
                           formula="v_R3")
    # write SBML file
    write_and_check(doc, sbml_file)

    """
    # Should be possible to just reuse the FBA model for the update,
    # but comp has problems with the fbc part.

    reader = SBMLReader()
    doc = reader.readSBMLFromFile(fba_file)

    # model
    model = doc.getModel()
    model.setId("toy_ode_update")
    model.setName("ODE metabolite update submodel")

    # boundary conditions of FBA have to be released
    # for the update of the species
    s_A = model.getSpecies("A")
    s_A.setBoundaryCondition(False)
    s_C = model.getSpecies("C")
    s_C.setBoundaryCondition(False)

    # write SBML file
    write_and_check(doc, sbml_file)
    """

####################################################
# ODE/SSA model
####################################################
def create_ode_model(sbml_file):
    """" Kinetic submodel (coupled model to FBA). """
    sbmlns = SBMLNamespaces(3, 1)
    doc = SBMLDocument(sbmlns)

    # model
    model = doc.createModel()
    model.setId("toy_ode_model")
    model.setName("ODE/SSA submodel")

    # compartments
    c_ext = create_compartment(model, cid="extern", name="external compartment")

    # boundary species
    s_C = create_species(model, sid="C", name="C", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())
    s_D = create_species(model, sid="D", name="D", initialAmount=0, constant=False,
                        boundaryCondition=False, compartment=c_ext.getId())

    # parameters
    create_parameter(model, pid="k_R4", name="k R4", constant=True, value=0.1)

    # kinetic reaction (MMK)
    r_R4 = create_reaction(model, rid="R4", name="C -> D", fast=False, reversible=False,
                           reactants={"C": 1}, products={"D": 1},
                           formula="k_R4*C")

    # write SBML file
    write_and_check(doc, sbml_file)

####################################################
# Comp model
####################################################
# Combined comp model of all the kinetic parts.
# - bounds calculation
# - metabolite updates
# - kinetic submodel
def create_ode_comp(sbml_file):
    """" Kinetic comp model """
    pass


if __name__ == "__main__":
    # write & check sbml
    from toymodel_settings import fba_file, ode_bounds_file, ode_update_file
    from toymodel_settings import ode_model_file

    create_fba(fba_file)
    create_ode_bounds(ode_bounds_file)
    create_ode_update(ode_update_file, fba_file)
    create_ode_model(ode_model_file)
