'''
Reading the FBC v1 information in cobrapy.

@author: Matthias Koenig
@date: 2015-03-25
'''
import cobra
from libsbml import readSBML
import pandas as pd
from pandas import DataFrame
import warnings

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
def get_flux_bounds_from_fbc(sbml):
    bound_ids = []
    reaction_ids = []
    operations = []
    values = []
    
    doc = readSBML(sbml)
    if (doc.getPlugin("fbc") != None):
        model_sbml = doc.getModel()
        mplugin = model_sbml.getPlugin("fbc");
    
        # read flux bounds
        for fb in mplugin.getListOfFluxBounds():
            bound_ids.append(fb.getId())
            reaction_ids.append(fb.getReaction())
            operations.append(fb.getOperation())
            values.append(fb.getValue())
                
    # create pandas
    return DataFrame({'reaction': reaction_ids, 
                      'operation': operations, 
                      'value' :values}, index=bound_ids)

def _set_flux_bounds_from_fbc(model, sbml):
    ''' Sets all the lower and upper fluxbounds from the FBC. '''
    
    bounds_df = get_flux_bounds_from_fbc(sbml)
    for index, row in bounds_df.iterrows():
        reaction_id = row.reaction
        reaction = model.reactions.get_by_id(reaction_id)
        if row.operation in ['less', 'lessEqual']:
            reaction.upper_bound = row.value
        elif row.operation in ['greater', 'greaterEqual']:
            reaction.lower_bound = row.value
        else:
            warnings.warn('Operation not supported on FluxBound: {}'.format(row.operation))


def set_flux_bounds(model, bounds_df):
    for reaction_id, row in bounds_df.iterrows():
        reaction = model.reactions.get_by_id(reaction_id)
        reaction.upper_bound = row.upperBounds
        reaction.lower_bound = row.lowerBounds


def set_max_bound(model, max_value):
    for reaction in model.reactions:
        reaction.upper_bound = min(reaction.upper_bound, max_value)
        reaction.lower_bound = max(reaction.lower_bound, -max_value)
    

def print_full_df(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

def print_flux_bounds(model):
    ''' Prints flux bounds for all reactions. '''
    print '*'*80
    info = []
    for r in model.reactions:
        info.append([r.id, r.lower_bound, r.upper_bound])
        df = DataFrame(info, columns=['id', 'lw', 'ub'])
    print_full_df(df)
    print '*'*80


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
            rule = rule.replace(')', '')
            # get reaction and set gene rules
            reaction_id = ga.getReaction() 
            reaction = model_cobra.reactions.get_by_id(reaction_id)
            # print reaction_id, ':', rule
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
    _set_flux_bounds_from_fbc(model, sbml)
    
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
        
    bounds_df = get_flux_bounds_from_fbc(sbml)
    print bounds_df