"""
Module creates the sub models and combined comp model for the toy model case.
The toy model consists of FBA submodels, deterministic ODE models and stochastic
ODE models.
The SBML comp extension is used for hierarchical model composition, i.e. to create
the main model and the kinetic model subparts.
"""
from libsbml import *
import multiscale.sbmlutils.io as sbml_io
import multiscale.sbmlutils.annotation as sbml_annotation

from multiscale.sbmlutils.factory import *

########################################################################
# General model information
########################################################################
version = 2
notes = XMLNode.convertStringToXMLNode("""
    <body xmlns='http://www.w3.org/1999/xhtml'>
    <h1>Wholecell Toy Model</h1>
    <h2>Description</h2>
    <p>This is a toy model for coupling models with different modeling frameworks via comp.</p>

    <div class="dc:publisher">This file has been produced by
      <a href="https://livermetabolism.com/contact.html" title="Matthias Koenig" target="_blank">Matthias Koenig</a>.
      </div>

    <h2>Terms of use</h2>
      <div class="dc:rightsHolder">Copyright © 2015 Wholecell Consortium.</div>
      <div class="dc:license">
      <p>Redistribution and use of any part of this model, with or without modification, are permitted provided that
      the following conditions are met:
        <ol>
          <li>Redistributions of this SBML file must retain the above copyright notice, this list of conditions
              and the following disclaimer.</li>
          <li>Redistributions in a different form must reproduce the above copyright notice, this list of
              conditions and the following disclaimer in the documentation and/or other materials provided
          with the distribution.</li>
        </ol>
        This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
             the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.</p>
      </div>
    </body>
""")
creators = {
    # id : ('FamilyName', 'GivenName', 'Email', 'Organization')
    'mk': ('Koenig', 'Matthias', 'konigmatt@googlemail.com', 'Charite Berlin'),
}
main_units = {
    'time': 's',
    'extent': UNIT_KIND_ITEM,
    'substance': UNIT_KIND_ITEM,
    'length': 'm',
    'area': 'm2',
    'volume': 'm3',
}
units = {
    's': [(UNIT_KIND_SECOND, 1.0)],
    'kg': [(UNIT_KIND_KILOGRAM, 1.0)],
    'm': [(UNIT_KIND_METRE, 1.0)],
    'm2': [(UNIT_KIND_METRE, 2.0)],
    'm3': [(UNIT_KIND_METRE, 3.0)],
    'mM': [(UNIT_KIND_MOLE, 1.0, 0),
           (UNIT_KIND_METRE, -3.0)],
    'item_per_s': [(UNIT_KIND_ITEM, 1.0),
                   (UNIT_KIND_SECOND, -1.0)],
}
########################################################################

def add_generic_info(model):
    sbml_annotation.set_model_history(model, creators)
    create_unit_definitions(model, units)
    set_main_units(model, units)
    model.setNotes(notes)


# TODO: put all the code in the factory

# def create_assignment_rule(model, sid, formula):
#     rule = model.createAssignmentRule()
#     rule.setVariable(sid)
#     astnode = parseL3FormulaWithModel(formula, model)
#     if not astnode:
#         print('Formula could not be parsed:', formula)
#         print(getLastParseL3Error())
#     rule.setMath(astnode)
#     return rule
#
#
# def create_rate_rule(model, sid, formula):
#     rule = model.createRateRule()
#     rule.setVariable(sid)
#     astnode = parseL3FormulaWithModel(formula, model)
#     if not astnode:
#         print('Formula could not be parsed:', formula)
#         print(getLastParseL3Error())
#     rule.setMath(astnode)
#     return rule



# def create_reaction(model, rid, name, fast=False, reversible=True,
#                      reactants={}, products={}, formula=None):
#     """ Create basic reaction structure. """
#     r = model.createReaction()
#     r.setId(rid)
#     r.setName(name)
#     r.setFast(fast)
#     r.setReversible(reversible)
#
#     for sid, stoichiometry in reactants.iteritems():
#         rt = r.createReactant()
#         rt.setSpecies(sid)
#         rt.setStoichiometry(abs(stoichiometry))
#         rt.setConstant(True)
#
#     for sid, stoichiometry in products.iteritems():
#         rt = r.createProduct()
#         rt.setSpecies(sid)
#         rt.setStoichiometry(abs(stoichiometry))
#         rt.setConstant(True)
#
#     if formula:
#         # set formula in reaction
#         astnode = parseL3FormulaWithModel(formula, model)
#         if not astnode:
#             print('Formula could not be parsed:', formula)
#             print(getLastParseL3Error())
#         law = r.createKineticLaw()
#         law.setMath(astnode)
#
#     return r
#
#
# def set_flux_bounds(r, lb="lb", ub="ub"):
#     """ Set flux bounds on given reaction. """
#     rplugin = r.getPlugin("fbc")
#     rplugin.setLowerFluxBound(lb)
#     rplugin.setUpperFluxBound(ub)
#
#
# def create_objective(mplugin, oid, otype, fluxObjectives, active=True):
#     objective = mplugin.createObjective()
#     objective.setId(oid)
#     objective.setType(otype)
#     if active:
#         mplugin.setActiveObjectiveId("R3_maximize")
#     for rid, coefficient in fluxObjectives.iteritems():
#         fluxObjective = objective.createFluxObjective()
#         fluxObjective.setReaction(rid)
#         fluxObjective.setCoefficient(coefficient)
#     return objective



####################################################
# ODE flux bounds
####################################################
def create_ode_bounds(sbml_file):
    """"
    Submodel for dynamically calculating the flux bounds.
    """
    sbmlns = SBMLNamespaces(3, 1)
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()
    model.setId("toy_ode_bounds")
    model.setName("ODE bound calculation submodel")

    parameters = [
        {A_ID: 'ub_R1', A_VALUE: 1.0, A_UNIT: None, A_NAME: 'ub_r1', A_CONSTANT: False},
        {A_ID: 'k1', A_VALUE: -0.2, A_UNIT: None, A_NAME: "k1", A_CONSTANT: False},
    ]

    rate_rules = [
        {A_ID: "ub_R1", A_VALUE: "k1*ub_R1"}
    ]

    # create model
    create_parameters(model, parameters)
    create_rate_rules(model, rate_rules)

    sbml_io.write_and_check(doc, sbml_file)


####################################################
# FBA submodel
####################################################
def create_fba(sbml_file):
    """
    Create the fba model.
    FBA submodel in FBC v2 which uses parameters as flux bounds.

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
    model.setId('toy_fba')
    model.setName('FBA submodel')

    compartments = [
        {A_ID: 'extern', A_VALUE: 1.0, A_UNIT:"m3", A_NAME: 'external compartment', A_SPATIAL_DIMENSION: 3},
        {A_ID: 'cell', A_VALUE: 1.0, A_UNIT:"m3", A_NAME: 'cell', A_SPATIAL_DIMENSION: 3}
    ]
    create_compartments(model, compartments)
    
    # boundary species
    species = [
        {A_ID: "A", A_NAME: "A", A_VALUE:10.0, constant=False,
                        boundaryCondition=True, compartment=c_ext.getId())
    s_C = create_species(model, sid="C", name="C", initialAmount=0, constant=False,
                        boundaryCondition=True, compartment=c_ext.getId()}
    ]

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
    sbml_io.write_and_check(doc_fba, sbml_file)

####################################################
# ODE species update
####################################################
# model for update of species count
def create_ode_update(sbml_file, fba_file):
    """ Submodel for dynamically updating the metabolite count.
        Very similar model to the FBA model.
    """
    # read FBA model
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
    sbml_io.write_and_check(doc, sbml_file)

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
    sbml_io.write_and_check(doc, sbml_file)

####################################################
# Comp model
####################################################
# Combined comp model of all the kinetic parts.
# - bounds calculation
# - metabolite updates
# - kinetic submodel


def create_ode_comp(sbmlfile):
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
