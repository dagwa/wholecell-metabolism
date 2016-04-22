
from libsbml import *
XMLOutputStream.setWriteTimestamp(False)
import sbmlutils.sbmlio as sbml_io
import sbmlutils.annotation as sbml_annotation
import sbmlutils.comp as comp

from sbmlutils.factory import *

import fba_model


####################################################
# ODE flux bounds
####################################################
def create_bounds_model(sbml_out, sbml_fba):
    """"
    Submodel for dynamically calculating the flux bounds.
    The dynamically changing flux bounds are the input to the
    FBA model.

    The information from the fba model is used to calculate the bounds.
    In additon the kcat & other information is needed.

    """
    # TODO: read kcat+ and kcat-, kex- and kex+ from database and use for flux bound calculation
    doc_fba = libsbml.readSBMLFromFile(sbml_fba)
    model_fba = doc_fba.getModel()

    model_id = 'WCM_3_10_metabolism__bounds'
    model_name = 'WCM_3_10_metabolism bound calculation'

    sbmlns = SBMLNamespaces(3, 1, 'comp', 1)
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)
    fba_model.set_model_information(model, model_id=model_id, model_name=model_name)

    # --- compartments ---
    # add all compartments from fba network (needed for species)
    for c in model_fba.getListOfCompartments():
        model.addCompartment(c)
        cid = c.getId()
        comp._create_port(model, pid="{}_port".format(cid), idRef=cid, portType=comp.PORT_TYPE_PORT)

    # --- species ---
    # add all species from  fba network
    for s in model_fba.getListOfSpecies():
        model.addSpecies(s)
        sid = s.getId()
        comp._create_port(model, pid="{}_port".format(sid), idRef=sid, portType=comp.PORT_TYPE_PORT)

    # --- parameters ---
    # upper and lower bound parameters
    for p in model_fba.getListOfParameters():
        pid = p.getId()
        if pid.startswith('ub_') or pid.startswith('lb_'):
            model.addParameter(p)
            comp._create_port(model, pid="{}_port".format(pid), idRef=pid, portType=comp.PORT_TYPE_PORT)

    # TODO: add kcat and exchange parameters
    for r in model_fba.getListOfReactions():
        rid = r.getId()

    parameters = [
        {A_ID: 'm', A_VALUE: 1.0, A_UNIT: "kg", A_NAME: "cell mass", A_CONSTANT: False},
        {A_ID: 'kcat_fw', A_VALUE: 10, A_UNIT: "per_s", A_NAME: "kcat forward", A_CONSTANT: False},
        {A_ID: 'kcat_bw', A_VALUE: 10, A_UNIT: "per_s", A_NAME: "kcat backward", A_CONSTANT: False},
        {A_ID: 'kex_fw', A_VALUE: 10, A_UNIT: "per_s", A_NAME: "kcat forward", A_CONSTANT: False},
        {A_ID: 'kex_bw', A_VALUE: 10, A_UNIT: "per_s", A_NAME: "kcat backward", A_CONSTANT: False},
    ]
    create_parameters(model, parameters)

    # rate rule for upper and lower bounds
    rate_rules = []
    for r in model_fba.getListOfReactions():
        rid = r.getId()
        enzyme = None
        modifiers = r.getListOfModifiers()
        if len(modifiers) is 1:
            enzyme = modifiers.get(0)
        # catalyzed reaction
        if enzyme:
            ub_formula = '{} * {}'.format('kcat_fw', enzyme)
            lb_formula = '-{} * {}'.format('kcat_bw', enzyme)
        # chemical reaction
        else:
            # TODO: piecewise by enzyme
            ub_formula = '{} * 1 item'.format('kcat_fw')
            lb_formula = '-{} * 1 item'.format('kcat_bw')

        print(ub_formula)
        print(lb_formula)
        rate_rules.append({A_ID: 'ub_{}'.format(rid), A_VALUE: ub_formula})
        rate_rules.append({A_ID: 'lb_{}'.format(rid), A_VALUE: lb_formula})

    create_rate_rules(model, rate_rules)

    sbml_io.write_and_check(doc, sbml_out)

'''
####################################################
# ODE species update
####################################################
def create_ode_update(sbml_out, sbml_fba):
    """
        Submodel for dynamically updating the metabolite count.
        This updates the ode model based on the FBA fluxes.
    """
    model_id = 'WCM_3_10_metabolism__update'
    model_name = 'WCM_3_10_metabolism species update'

    sbmlns = SBMLNamespaces(3, 1, 'comp', 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)

    # model
    model = doc.createModel()
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)
    fba_model.set_model_information(model, model_id=model_id, model_name=model_name)

    # compartments
    create_compartments(model, fba_model.compartments)

    # TODO: add boundary species (? how to find, how done in the toy example)
    # only update the boundarySpecies in the reactions
    species = [
        {A_ID: 'C', A_NAME: "C", A_VALUE: 0, A_UNIT: fba_model.UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "extern", A_BOUNDARY_CONDITION: False},
    ]
    create_species(model, species)

    # TODO
    parameters = [
        {A_ID: "vR3", A_NAME: "vR3 (FBA flux)", A_VALUE: 0.1, A_CONSTANT: True, A_UNIT: "item_per_s"}
    ]
    create_parameters(model, parameters)

    # TODO: kinetic reactions
    # kinetic reaction (MMK)
    r4 = create_reaction(model, rid="R3", name="-> C", fast=False, reversible=False,
                         reactants={}, products={"C": 1}, formula="vR3", compartment="extern")

    # TODO: ports
    comp._create_port(model, pid="vR3_port", idRef="vR3", portType=comp.PORT_TYPE_PORT)
    comp._create_port(model, pid="C_port", idRef="C", portType=comp.PORT_TYPE_PORT)
    comp._create_port(model, pid="extern_port", idRef="extern", portType=comp.PORT_TYPE_PORT)

    # write SBML file
    sbml_io.write_and_check(doc, sbml_out)

####################################################
# Top level
####################################################
# TODO:
def create_top_level_model(sbml_file):
    """
    Creates the full comp model as combination of FBA and comp models.
    """
    sbmlns = SBMLNamespaces(3, 1, "comp", 1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("comp", True)
    doc.setPackageRequired("fbc", False)

    mdoc = doc.getPlugin("comp")

    # create listOfExternalModelDefinitions
    emd_bounds = comp.create_ExternalModelDefinition(mdoc, "toy_ode_bounds", sbml_file=ode_bounds_file)
    emd_fba = comp.create_ExternalModelDefinition(mdoc, "toy_fba", sbml_file=fba_file)
    emd_update = comp.create_ExternalModelDefinition(mdoc, "toy_ode_update", sbml_file=ode_update_file)
    emd_model = comp.create_ExternalModelDefinition(mdoc, "toy_ode_model", sbml_file=ode_model_file)

    # create models and submodels
    model = doc.createModel()
    model.setId("toy_top_level")
    model.setName("Top level model")
    model_factory.add_generic_info(model)
    mplugin = model.getPlugin("comp")
    model.setSBOTerm(comp.SBO_CONTINOUS_FRAMEWORK)

    # add submodel which references the external model definition
    comp.add_submodel_from_emd(mplugin, submodel_sid="bounds", emd=emd_bounds)
    comp.add_submodel_from_emd(mplugin, submodel_sid="fba", emd=emd_fba)
    comp.add_submodel_from_emd(mplugin, submodel_sid="update", emd=emd_update)
    comp.add_submodel_from_emd(mplugin, submodel_sid="model", emd=emd_model)

    # --- compartments ---
    create_compartments(model, [
        {A_ID: "extern", A_NAME: "external compartment", A_VALUE: 1.0, A_SPATIAL_DIMENSION: 3, A_UNIT: model_factory.UNIT_VOLUME},
        {A_ID: 'cell', A_NAME: 'cell', A_VALUE: 1.0, A_SPATIAL_DIMENSION: 3, A_UNIT: model_factory.UNIT_VOLUME}
    ])

    # --- species ---
    # replaced species
    # (fba species are not replaced, because they need their boundaryConditions for the FBA,
    #    and do not depend on the actual concentrations)
    create_species(model, [
        {A_ID: 'C', A_NAME: "C", A_VALUE: 0, A_UNIT: model_factory.UNIT_AMOUNT, A_HAS_ONLY_SUBSTANCE_UNITS: True,
         A_COMPARTMENT: "extern"},
    ])

    # --- parameters ---
    create_parameters(model, [
        # bounds
        {A_ID: 'ub_R1', A_VALUE: 1.0, A_UNIT: model_factory.UNIT_FLUX, A_NAME: 'ub_R1', A_CONSTANT: False},
        {A_ID: "vR3", A_NAME: "vR3 (FBA flux)", A_VALUE: 0.1, A_UNIT: model_factory.UNIT_FLUX, A_CONSTANT: False}
    ])

    # --- reactions ---
    # dummy reaction in top model
    r1 = create_reaction(model, rid="R3", name="R3 dummy", fast=False, reversible=True,
                         reactants={}, products={"C": 1}, compartment="extern")
    # assignment rule
    create_assignment_rules(model, [{A_ID: "vR3", A_VALUE: "R3"}])

    # --- replacements ---
    # replace compartments
    comp.replace_elements(model, 'extern', ref_type=comp.SBASE_REF_TYPE_PORT,
                          replaced_elements={'fba': ['extern_port'], 'update': ['extern_port'], 'model': ['extern_port']})

    comp.replace_elements(model, 'cell', ref_type=comp.SBASE_REF_TYPE_PORT,
                          replaced_elements={'fba': ['cell_port']})

    # replace parameters
    comp.replace_elements(model, 'ub_R1', ref_type=comp.SBASE_REF_TYPE_PORT,
                          replaced_elements={'bounds': ['ub_R1_port'], 'fba': ['ub_R1_port']})
    comp.replace_elements(model, 'vR3', ref_type=comp.SBASE_REF_TYPE_PORT,
                          replaced_elements={'update': ['vR3_port']})

    # replace species
    comp.replace_elements(model, 'C', ref_type=comp.SBASE_REF_TYPE_PORT,
                          replaced_elements={'fba': ['C_port'], 'update': ['C_port'],
                                             'model': ['C_port']})

    # replace reaction by fba reaction
    comp.replaced_by(model, 'R3', ref_type=comp.SBASE_REF_TYPE_PORT,
                     submodel='fba', replaced_by="R3_port")

    # replace units
    for uid in ['s', 'kg', 'm3', 'm2', 'mM', 'item_per_m3', 'm', 'per_s', 'item_per_s']:
        comp.replace_element_in_submodels(model, uid, ref_type=comp.SBASE_REF_TYPE_UNIT,
                                      submodels=['bounds', 'fba', 'update', 'model'])

    # write SBML file
    sbml_io.write_and_check(doc, sbml_file)
'''

######################################################################
if __name__ == "__main__":
    from metabolism_settings import RESULTS_DIR, VERSION
    import fba_model

    sbml_fba = os.path.join(RESULTS_DIR, "wholecell_fba_{}.xml".format(VERSION))
    # fba_model.create_fba_model(sbml_fba)

    sbml_bounds = os.path.join(RESULTS_DIR, "wholecell_bounds_{}.xml".format(VERSION))
    create_bounds_model(sbml_bounds, sbml_fba)

