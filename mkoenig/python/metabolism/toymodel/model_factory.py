# -*- coding=utf-8 -*-
"""
Module creates the sub models and combined comp model for the toy model case.
The toy model consists of FBA submodels, deterministic ODE models and stochastic
ODE models.
The SBML comp extension is used for hierarchical model composition, i.e. to create
the main model and the kinetic model subparts.
"""
from libsbml import *
import multiscale.sbmlutils.sbmlio as sbml_io
import multiscale.sbmlutils.annotation as sbml_annotation
import multiscale.sbmlutils.comp as comp

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
      <div class="dc:rightsHolder">Copyright © 2016 Wholecell Consortium.</div>
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
    'mk': {'FamilyName':'Koenig', 'GivenName': 'Matthias',
           'Email': 'konigmatt@googlemail.com', 'Organization': 'Humboldt University Berlin'},
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
    'per_s': [(UNIT_KIND_SECOND, -1.0)],
    'item_per_s': [(UNIT_KIND_ITEM, 1.0),
                   (UNIT_KIND_SECOND, -1.0)],
    'item_per_m3': [(UNIT_KIND_ITEM, 1.0),
                    (UNIT_KIND_METRE, -3.0)],
}

UNIT_AMOUNT = UNIT_KIND_ITEM
UNIT_VOLUME = 'm3'
UNIT_CONCENTRATION = 'item_per_m3'
UNIT_FLUX = 'item_per_s'


########################################################################
def add_generic_info(model):
    sbml_annotation.set_model_history(model, creators)
    create_unit_definitions(model, units)
    set_main_units(model, main_units)
    model.setNotes(notes)


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
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)
    add_generic_info(model)

    parameters = [
        {A_ID: 'ub_R1', A_VALUE: 1.0, A_UNIT: UNIT_FLUX, A_NAME: 'ub_r1', A_CONSTANT: False},
        {A_ID: 'k1', A_VALUE: -0.2, A_UNIT: "per_s", A_NAME: "k1", A_CONSTANT: False},
    ]
    create_parameters(model, parameters)

    rate_rules = [
        {A_ID: "ub_R1", A_VALUE: "k1*ub_R1"}
    ]
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
    sbmlns = SBMLNamespaces(3, 1, "fbc", 2)
    sbmlns.addPackageNamespace("comp", 1)  # comp for port definition
    doc_fba = SBMLDocument(sbmlns)
    doc_fba.setPackageRequired("comp", True)
    mdoc = doc_fba.getPlugin("comp")

    doc_fba.setPackageRequired("fbc", False)
    model = doc_fba.createModel()
    mplugin = model.getPlugin("fbc")
    mplugin.setStrict(False)

    # model
    model.setId('toy_fba')
    model.setName('FBA submodel')
    model.setSBOTerm(comp.SBO_FLUX_BALANCE_FRAMEWORK)
    add_generic_info(model)

    compartments = [
        {A_ID: 'extern', A_VALUE: 1.0, A_UNIT: UNIT_VOLUME, A_NAME: 'external compartment', A_SPATIAL_DIMENSION: 3},
        {A_ID: 'cell', A_VALUE: 1.0, A_UNIT: UNIT_VOLUME, A_NAME: 'cell', A_SPATIAL_DIMENSION: 3}
    ]
    create_compartments(model, compartments)

    species = [
        # external
        {A_ID: 'A', A_NAME: "A", A_VALUE: 10, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "extern", A_BOUNDARY_CONDITION: True},
        {A_ID: 'C', A_NAME: "C", A_VALUE: 0, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "extern", A_BOUNDARY_CONDITION: True},
        # internal
        {A_ID: 'B1', A_NAME: "B1", A_VALUE: 0, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "cell"},
        {A_ID: 'B2', A_NAME: "B2", A_VALUE: 0, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "cell"},
    ]
    create_species(model, species)

    parameters = [
        # bounds
        {A_ID: "ub_R1", A_NAME: "ub R1", A_VALUE: 1.0, A_UNIT: UNIT_FLUX, A_CONSTANT: False},
        {A_ID: "lb", A_NAME: "lower bound", A_VALUE: 0.0, A_UNIT: UNIT_FLUX, A_CONSTANT: True},
        {A_ID: "ub", A_NAME: "upper bound", A_VALUE: 1000.0, A_UNIT: UNIT_FLUX, A_CONSTANT: True},
    ]
    create_parameters(model, parameters)

    # reactions with constant flux
    r1 = create_reaction(model, rid="R1", name="A import (R1)", fast=False, reversible=True,
                           reactants={"A": 1}, products={"B1": 1})
    r2 = create_reaction(model, rid="R2", name="B1 <-> B2 (R2)", fast=False, reversible=True,
                           reactants={"B1": 1}, products={"B2": 1})
    r3 = create_reaction(model, rid="R3", name="B2 export (R3)", fast=False, reversible=True,
                           reactants={"B2": 1}, products={"C": 1})

    # flux bounds
    set_flux_bounds(r1, lb="lb", ub="ub_R1")
    set_flux_bounds(r2, lb="lb", ub="ub")
    set_flux_bounds(r3, lb="lb", ub="ub")

    # objective function
    create_objective(mplugin, oid="R3_maximize", otype="maximize", fluxObjectives={"R3": 1.0})

    # create ports
    # <comp:listOfPorts>
    #   <comp:port comp:id="R3_port" comp:idRef="R3" sboTerm="SBO:0000599"/>
    #  <comp:port comp:id="r1_port_1" comp:idRef="r1" sboTerm="SBO:0000599"/>
    # </comp:listOfPorts>
    comp._create_port(model, pid="R3_port", idRef="R3")
    comp._create_port(model, pid="r1_port_1", idRef="ub_R1")


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
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)

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
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)
    add_generic_info(model)

    compartments = [
        {A_ID: 'extern', A_VALUE: 1.0, A_UNIT: "m3", A_NAME: 'external compartment', A_SPATIAL_DIMENSION: 3},
    ]
    create_compartments(model, compartments)

    species = [
        # external
        {A_ID: 'C', A_NAME: "C", A_VALUE: 0, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS:True,
         A_COMPARTMENT: "extern", A_BOUNDARY_CONDITION: False},
        {A_ID: 'D', A_NAME: "D", A_VALUE: 0, A_UNIT: UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS:True,
         A_COMPARTMENT: "extern", A_BOUNDARY_CONDITION: False},
    ]
    create_species(model, species)

    parameters = [
        {A_ID: "k_R4", A_NAME: "k R4", A_VALUE: 0.1, A_CONSTANT: True, A_UNIT: "per_s"}
    ]
    create_parameters(model, parameters)

    # kinetic reaction (MMK)
    r4 = create_reaction(model, rid="R4", name="C -> D", fast=False, reversible=False,
                         reactants={"C": 1}, products={"D": 1}, formula="k_R4*C")

    # write SBML file
    sbml_io.write_and_check(doc, sbml_file)


########################################################################################################################
if __name__ == "__main__":
    # write & check sbml
    from settings import fba_file, ode_bounds_file, ode_update_file
    from settings import ode_model_file

    create_fba(fba_file)
    create_ode_bounds(ode_bounds_file)
    create_ode_update(ode_update_file, fba_file)
    create_ode_model(ode_model_file)

    import multiscale.multiscalesite.simapp.db.api as db_api
    db_api.create_model(fba_file, model_format=db_api.CompModelFormat.SBML)
    db_api.create_model(ode_bounds_file, model_format=db_api.CompModelFormat.SBML)
    db_api.create_model(ode_update_file, model_format=db_api.CompModelFormat.SBML)
    db_api.create_model(ode_model_file, model_format=db_api.CompModelFormat.SBML)
