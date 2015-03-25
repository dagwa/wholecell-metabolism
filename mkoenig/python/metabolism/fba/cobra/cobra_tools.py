'''
Reading the FBC v1 information in cobrapy.

@author: Matthias Koenig
@date: 2015-03-25
'''
import cobra
from libsbml import readSBML

# -----------------------------------------------------------------------------
# Objectives & coefficients
# -----------------------------------------------------------------------------

def _get_active_objective_from_fbc(sbml):
    ''' Read the active objective from fbc model. '''
    obj_coefs = {}
    doc = readSBML(sbml)
    if (doc.getPlugin("fbc") != None):
        model_sbml = doc.getModel()
        mplugin = model_sbml.getPlugin("fbc");
        # read the active objective        
        objective = mplugin.getActiveObjective()
        for fluxobj in objective.getListOfFluxObjectives():
            obj_coefs[fluxobj.getReaction()] = fluxobj.getCoefficient()
    return obj_coefs

def reset_objective_coefficients(model):
    ''' Set all objective coefficents to 0.0. '''
    for reaction in model.reactions:
        reaction.objective_coefficient = 0
    return None

    
def set_objective_coefficients(model, coefs_dict):
    ''' Set given objective coefficients. '''
    reset_objective_coefficients(model)
    for rid, coef in coefs_dict.iteritems():
        reaction = model.reactions.get_by_id(rid)
        reaction.objective_coefficient = coef

    
def get_objective_coefficients(model):
    ''' Returns objective coefficents != 0. '''
    return {reaction: reaction.objective_coefficient for reaction in model.reactions if reaction.objective_coefficient != 0}

# -----------------------------------------------------------------------------
# FluxBounds
# -----------------------------------------------------------------------------
def get_flux_bounds_from_fbc():
    # TODO: implement
    pass


# -----------------------------------------------------------------------------
# GeneAssociations
# -----------------------------------------------------------------------------
def _set_gene_associations_from_fbc(model_cobra, sbml):
    ''' Read GeneAssociations from FBC v1. '''
    doc = readSBML(sbml)
    if (doc.getPlugin("fbc") != None):
        model_sbml = doc.getModel()
        mplugin = model_sbml.getPlugin("fbc");
        for ga in mplugin.getListOfGeneAssociations():        
            # get the rule string for cobrapy        
            ass = ga.getAssociation()
            rule = ass.toInfix()
            rule = rule.replace('(', '')
            rule = rule.replace('(', '')
            # get reaction and set gene rules
            reaction_id = ga.getReaction() 
            reaction = model_cobra.reactions.get_by_id(reaction_id)
            reaction.gene_reaction_rule = rule # property takes care of all the logic
    return model_cobra


# -----------------------------------------------------------------------------
# Reading FBC 1 Models
# -----------------------------------------------------------------------------
def read_sbml_fbc_model(sbml):
    model = cobra.io.read_sbml_model(sbml)
    # gene associations
    _set_gene_associations_from_fbc(model, sbml)
    # objective coefficients
    obj_coefs = _get_active_objective_from_fbc(sbml)
    set_objective_coefficients(model, obj_coefs)
    # flux bounds
    
    return model
    
##############################################################################
if __name__ == '__main__':
    import os
    from metabolism_settings import VERSION, RESULTS_DIR 
    sbml = os.path.join(RESULTS_DIR, "Metabolism_matrices_{}_L3V1.xml".format(VERSION))

    # model = cobra.io.read_sbml_model(sbml)
    model = read_sbml_fbc_model(sbml) # with additional FBC v1 information

    print len(model.reactions)    # 504
    print len(model.metabolites)  # 479  (336 + 104 -1) species + proteins - protein_species
    print len(model.genes)        # 142  (104)
    print model.compartments

    for key, value in get_objective_coefficients(model).iteritems():
        print key, value    