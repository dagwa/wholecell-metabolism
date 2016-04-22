# -*- coding=utf-8 -*-
"""
Creates the SBML for the metabolism submodul
from matlab matrices and database information.

@author mkoenig
"""
# TODO: handle exchange reactions correctly
# TODO: run example simulations (in cobra)


# TODO: update the annotations/integrate


# TODO: create the kinetic ode model/comp model
# TODO: glycolysis subnetwork




from __future__ import print_function, division
from libsbml import *
import libsbml
from sbmlutils.validation import validate_sbml, check
import sbmlutils.sbmlio as sbml_io
import sbmlutils.annotation as sbml_annotation
import sbmlutils.comp as comp

from sbmlutils.factory import *

from metabolism_settings import RESULTS_DIR, VERSION        

import re

import pandas as pd
from pandas import DataFrame
import warnings

##########################################################################
# Model information
##########################################################################
version = VERSION
model_id = 'WCM_3_10'
model_name = 'Whole Cell 2015 - Metabolism'

notes = XMLNode.convertStringToXMLNode("""
    <body xmlns='http://www.w3.org/1999/xhtml'>
    <h1>Wholecell Model</h1>
    <h2>Description</h2>
    <p>This is the metabolic submodel of the Karr model.</p>

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

creators = [
    ['Bergmann', 'Frank', 'fbergman@caltech.edu', 'Caltech'],
    ['Koenig', 'Matthias', 'konigmatt@googlemail.com', 'Charite Berlin'],
    ['Smallbone', 'Kieran', 'kieran.smallbone@manchester.ac.uk', 'University of Manchester'],
    ['Tokic', 'Milenko', 'milenko.tokic@epfl.ch', 'EPFL'],
    ['Costa', 'Rafael', 'rcosta@kdbio.inesc-id.pt', 'University of Lisbon'],
    ['Baghalian', 'Kambiz', 'kambiz.baghalian@plants.ox.ac.uk', 'University of Oxford'],
]
creators = [dict(zip(['FamilyName', 'GivenName', 'Email', 'Organization'], entry)) for entry in creators]

cvterms = [
    [MODEL_QUALIFIER, BQM_IS_DERIVED_FROM, "http://identifiers.org/doi/10.1016/j.cell.2012.05.044"],
    [MODEL_QUALIFIER, BQM_IS_DESCRIBED_BY, "http://identifiers.org/pubmed/22817898"],
    [MODEL_QUALIFIER, BQM_IS, 'http://identifiers.org/mamo/MAMO_0000040'],  # metabolic network
    [MODEL_QUALIFIER, BQM_IS, 'http://identifiers.org/mamo/MAMO_0000009'],  # constraint-based model
    [BIOLOGICAL_QUALIFIER, BQB_HAS_TAXON, "http://identifiers.org/taxonomy/243273"],  # Mycoplasma genitalium G37
]
cvterms_df = DataFrame(data=cvterms, columns=['Qualifier', 'QualifierType', 'Resource'])

##########################################################################
# Units
##########################################################################
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
UNIT_AREA = 'm2'
UNIT_VOLUME = 'm3'
UNIT_CONCENTRATION = 'item_per_m3'
UNIT_FLUX = 'item_per_s'
########################################################################


def set_model_information(model):
    """ Writes the core information in the Karr model. """
    model.setId(model_id)
    model.setName(model_name)
    sbml_annotation.set_model_history(model, creators)
    create_unit_definitions(model, units)
    set_main_units(model, main_units)
    model.setNotes(notes)
    set_cv_terms(model, cvterms_df)


def set_cv_terms(model, cvterms_df):
    """ Set model cv terms from DataFrame. """
    from annotation import create_meta_id
    if not model.isSetMetaId():
        model.setMetaId(create_meta_id(model.getId()))

    # write all the annotations
    for index, row in cvterms_df.iterrows():
        qualifier = row.Qualifier
        qualifier_type = row.QualifierType
        resource = row.Resource

        cv = CVTerm()
        cv.setQualifierType(qualifier)
        if row.Qualifier == MODEL_QUALIFIER:
            cv.setModelQualifierType(qualifier_type)
        elif row.Qualifier == BIOLOGICAL_QUALIFIER:
            cv.setBiologicalQualifierType(qualifier_type)
        cv.addResource(resource)
        check(model.addCVTerm(cv), 'add cv term')


##########################################################################
# Compartments
##########################################################################
compartments = [
    # Compartment 'none' is added to account for pseudo-metabolites used in the FBA.
    {A_ID: 'c', A_NAME: 'cytosol', A_VALUE: 1.0, A_UNIT: UNIT_VOLUME, A_SPATIAL_DIMENSION: 3, A_CONSTANT: False},
    {A_ID: 'm', A_NAME: 'membrane', A_VALUE: 1.0, A_UNIT: UNIT_AREA, A_SPATIAL_DIMENSION: 2, A_CONSTANT: False},
    {A_ID: 'e', A_NAME: 'extracellular', A_VALUE: 1.0, A_UNIT: UNIT_VOLUME, A_SPATIAL_DIMENSION: 3, A_CONSTANT: False},
    {A_ID: 'n', A_NAME: 'none', A_VALUE: 1.0, A_UNIT: UNIT_VOLUME, A_SPATIAL_DIMENSION: 3, A_CONSTANT: False},
]


def create_metabolism_sbml():
    """ Creates the metabolic Karr model.

    :return:
    :rtype:
    """
    tol = 1E-12  # within this tolerance matrix elements are considered != 0

    ###########################################################################
    # Load matrix information
    ###########################################################################
    matrix_dir = os.path.join(RESULTS_DIR, 'fba_matrices')

    # handle the sodium NA id (not parsing as NaN)
    s_fba_df = pd.read_csv(os.path.join(matrix_dir, 's_fba.csv'), sep="\t",
                           keep_default_na=False, na_values='nan')
    s_fba_df.set_index('sid', inplace=True)

    r_fba_df = pd.read_csv(os.path.join(matrix_dir, 'r_fba.csv'), sep="\t")
    r_fba_df.set_index('rid', inplace=True)

    e_df = pd.read_csv(os.path.join(matrix_dir, 'e_fba.csv'), sep="\t")
    e_df.set_index('eid', inplace=True)

    mat_stoichiometry = pd.read_csv(os.path.join(matrix_dir, 'fbaReactionStoichiometryMatrix.csv'), sep="\t")
    mat_stoichiometry.set_index(s_fba_df.index, inplace=True)

    mat_catalysis = pd.read_csv(os.path.join(matrix_dir, 'fbaReactionCatalysisMatrix.csv'), sep="\t")
    mat_catalysis.set_index(r_fba_df.index, inplace=True)

    ###########################################################################
    # SBML Creation
    ###########################################################################
    # SBML model with FBC support
    sbmlns = SBMLNamespaces(3, 1)
    sbmlns.addPackageNamespace("fbc", 2)
    sbmlns.addPackageNamespace("comp", 1)

    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    mplugin = model.getPlugin("fbc")
    mplugin.setStrict(False)  # +-INF bounds in the Karr model

    # history & creators
    set_model_information(model)

    # compartments
    print('* create compartments *')
    create_compartments(model, compartments)

    # <metabolites>
    print('* create metabolites *')
    metabolite_count = 1000  # to get rid of warnings, not required for FBA simulations
    protein_count = 100

    # Metabolites in the FBA problem (rows) are encoded as species
    for index, row in s_fba_df.iterrows():
        species = [{
            A_ID: index,
            A_NAME: row['name'],
            A_VALUE: metabolite_count,
            A_UNIT: UNIT_AMOUNT,
            A_HAS_ONLY_SUBSTANCE_UNITS: True,
            A_COMPARTMENT: row['compartment'],
            A_BOUNDARY_CONDITION: False,
            A_CONSTANT: False
        }, ]
        sdict = create_species(model, species)
        s = sdict[index]

        # chemical formula and charge => for balance
        splugin = s.getPlugin("fbc")
        formula = row['formula']
        if not pd.isnull(formula):
            splugin.setChemicalFormula(formula)
        charge = row['charge']
        # string to int disaster due to NA name for sodium
        # this is ugly but works
        if not pd.isnull(charge) and len(charge) != 0:
            splugin.setCharge(int(float(charge)))

    # <proteins> [104]
    print('* create proteins *')
    # Proteins (ProteinMonomer & ProteinComplex) are encoded as species.
    # The main reasons are:
    #   1. The protein count is changing during the simulation and input of the dynamical flux bound
    #      calculation.
    #   2. The proteins can be encoded as modifiers of the respective reactions.
    #      This provides clarity for the reaction <- protein <- gene information
    #      And provides important information for possible visualization.

    def create_protein_species(sid, name):
        # TODO: proper way to find location of reactions & proteins
        # Not important for simulation, only for visualization
        species = [{
            A_ID: sid,
            A_NAME: name,
            A_VALUE: protein_count,
            A_UNIT: UNIT_AMOUNT,
            A_HAS_ONLY_SUBSTANCE_UNITS: True,
            A_COMPARTMENT: 'c',  # this is just fix
            A_BOUNDARY_CONDITION: False,
            A_CONSTANT: False
        }, ]
        sdict = create_species(model, species)
        return sdict[sid]

    # Create all Enzymes
    for index, row in e_df.iterrows():
        # check if the protein is already a species (due to involvment in reaction)
        s = model.getSpecies(index)
        if s is not None:
            print(index, 'is already species.')
        else:
            # create species for protein
            name = row['name']
            if pd.isnull(name):
                name = None
            create_protein_species(sid=index, name=name)

    # Handle the special case of ACP (MG_287_MONOMER - acyl carrier protein)
    # Neither substrate nor enzyme (no part of FBA, but part of fluxbound calculation)
    create_protein_species(sid='MG_287_MONOMER', name="acyl carrier protein")

    # <reactions>
    print('* create reactions *')
    # Reactions are all columns in the FBA stoichiometric matrix.
    # This includes some pseudo-reactions (internal & external exchange) which
    # are not represented in the knowledg base.
    for index, row in r_fba_df.iterrows():

        # create reaction
        r = model.createReaction()
        r.setId(index)
        name = row['name']
        if not pd.isnull(name):
            r.setName(name)
        r.setFast(False)

        # <reversibility>
        # The reversibility can be calculated from the reaction bounds. In some
        # cases the reversibility is in backward direction. This would require the
        # change of reactants and products for encoding. The SBML is strictly reproducing
        # the FBA problem, so that no reversibilities are defined in the SBML.
        # Reversibility is functional in the FBA via the actual flux bounds.
        r.setReversible(True)  # some are irreversible via Flux bounds in forward or backward direction

        # exchange reactions
        if index.startswith('metabolite_ex_'):
            r.setSBOTerm(627)  # exchange reaction

        # get FBC plugin
        rplugin = r.getPlugin("fbc")

        # set proteins as modifiers from catalysis matrix
        row = mat_catalysis.ix[index]
        row = row[row > tol]
        for eid, value in row.iteritems():
            # set protein as modifier
            mod = r.createModifier()
            mod.setSpecies(eid)

            # gene associations
            # TODO: create the ga objects !
            gene_str = e_df['genes'][eid]
            genes = [g.strip() for g in gene_str.split(',')]
            genes_formula = ' AND '.join(genes)

            gpa = rplugin.createGeneProductAssociation()
            gpa.setId('ga__{}__{}'.format(index, eid))
            gpa.setAssociation(genes_formula)

        # stoichiometry from stoichiometric matrix  # [376x504]
        # find non-zero elements in the reaction column
        col = mat_stoichiometry[index]
        col = col[abs(col) > tol]

        # exchange reactions are defined in export direction
        if index.startswith('metabolite_ex_'):
            col = -1.0 * col
        # create species references depending on stoichiometry
        for sid, stoichiometry in col.iteritems():
            if stoichiometry < 0:
                rt = r.createReactant()
                rt.setSpecies(sid)
                rt.setStoichiometry(abs(stoichiometry))
                rt.setConstant(True)
                if index.startswith('metabolite_ex_'):
                    r.setId("metabolite_ex_{}".format(sid))
            if stoichiometry > 0:
                pt = r.createProduct()
                pt.setSpecies(sid)
                pt.setStoichiometry(stoichiometry)
                pt.setConstant(True)

        # The reaction flux bounds are set as hard upper and lower flux bounds
        # These are NOT the dynamical flux bounds.
        # The actual flux bounds are calculated based on an outer model which evaluates AssignmentRules for the
        # parameters (taking into account metabolite counts, protein counts, kcat, ...)
        rid = r.getId()
        ub_id = 'ub_{}'.format(rid)
        lb_id = 'lb_{}'.format(rid)

        parameters = [
            # bounds
            {A_ID: lb_id, A_NAME: lb_id, A_VALUE: float("-inf"), A_UNIT: UNIT_FLUX, A_CONSTANT: False},
            {A_ID: ub_id, A_NAME: ub_id, A_VALUE: float("inf"), A_UNIT: UNIT_FLUX, A_CONSTANT: False},
        ]
        create_parameters(model, parameters)
        set_flux_bounds(r, lb=lb_id, ub=ub_id)

    # <objective function>
    print('* create objective function *')
    objective = mplugin.createObjective()
    objective.setId("growth")
    objective.setType("maximize")
    mplugin.setActiveObjectiveId("growth")
    for index, row in r_fba_df.iterrows():
        coeff = row['fbaObjective']
        if not pd.isnull(coeff) and abs(coeff) > tol:
            fluxObjective = objective.createFluxObjective()
            fluxObjective.setReaction(index)
            fluxObjective.setCoefficient(coeff)

    # write sbml
    sbml_out = os.path.join(RESULTS_DIR, "Metabolism_matrices_{}_L3V1.xml".format(VERSION))
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_out)
    print(sbml_out)

    # perform the full validation checks & validation
    print('* validate_sbml *')
    validate_sbml(sbml_out, ucheck=False)

    """
    # Conversion to cobra model
    print('* convert to cobra *')
    conversion_properties = libsbml.ConversionProperties()
    conversion_properties.addOption("convert fbc to cobra", True, "Convert FBC model to Cobra model")
    result = doc.convert(conversion_properties)
    if result != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise (Exception("Conversion of SBML+fbc to COBRA failed"))
    sbml_cobra = os.path.join(RESULTS_DIR, "Metabolism_matrices_cobra_{}_L3V1.xml".format(VERSION))
    print(sbml_cobra)
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_cobra)
    """

if __name__ == "__main__":
    create_metabolism_sbml()
