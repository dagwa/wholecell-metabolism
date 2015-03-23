'''
Write SBML for metabolism process directly from database information.

@author: Matthias Koeng
@date: 2015-03-09
'''
from libsbml import *
from public.models import Metabolite, Reaction, Protein, ProteinComplex
from public.models import Entry, Molecule
from sbml_tools.validator import SBMLValidator
from sbml_tools.checks import check


from metabolism_settings import DATA_DIR, RESULTS_DIR, VERSION        

import pandas as pd
from pandas import DataFrame

from libsbml import UNIT_KIND_SECOND, UNIT_KIND_MOLE,\
    UNIT_KIND_METRE, UNIT_KIND_KILOGRAM, UNIT_KIND_LITRE

#########################################################################
main_units = dict()
units = dict()
names = dict()
pars = []
external = []
assignments = []
rules = []

##########################################################################
# Units
########################################################################## 
main_units['time'] = 's'
main_units['extent'] = UNIT_KIND_MOLE
main_units['substance'] = UNIT_KIND_MOLE
main_units['length'] = 'm'
main_units['area'] = 'm2'
main_units['volume'] = 'm3'    
    
        
units['s'] = [(UNIT_KIND_SECOND, 1.0, 0)]
units['kg'] = [(UNIT_KIND_KILOGRAM, 1.0, 0)]
units['m'] = [(UNIT_KIND_METRE, 1.0, 0)]
units['m2'] = [(UNIT_KIND_METRE, 2.0, 0)]
units['m3'] = [(UNIT_KIND_METRE, 3.0, 0)]
units['per_s'] = [(UNIT_KIND_SECOND, -1.0, 0)]
units['mole_per_s'] = [(UNIT_KIND_MOLE, 1.0, 0), 
                       (UNIT_KIND_SECOND, -1.0, 0)]
units['mole_per_s_per_mM'] = [(UNIT_KIND_METRE, 3.0, 0),
                       (UNIT_KIND_SECOND, -1.0, 0) ]
units['mole_per_s_per_mM2'] = [(UNIT_KIND_MOLE, -1.0, 0), (UNIT_KIND_METRE, 6.0, 0), 
                       (UNIT_KIND_SECOND, -1.0, 0) ]

units['m_per_s'] = [(UNIT_KIND_METRE, 1.0, 0), 
                    (UNIT_KIND_SECOND, -1.0, 0)]
units['m2_per_s'] = [(UNIT_KIND_METRE, 2.0, 0), 
                    (UNIT_KIND_SECOND, -1.0, 0)]
units['m3_per_s'] = [(UNIT_KIND_METRE, 3.0, 0), 
                    (UNIT_KIND_SECOND, -1.0, 0)]
units['mM']       = [(UNIT_KIND_MOLE, 1.0, 0), 
                    (UNIT_KIND_METRE, -3.0, 0)]
units['mM_s']       = [(UNIT_KIND_MOLE, 1.0, 0), (UNIT_KIND_SECOND, 1.0, 0),
                    (UNIT_KIND_METRE, -3.0, 0)]

units['per_mM']   = [(UNIT_KIND_METRE, 3.0, 0), 
                    (UNIT_KIND_MOLE, -1.0, 0)]
units['per_m2']   = [(UNIT_KIND_METRE, -2.0, 0)]
units['per_m3']   = [(UNIT_KIND_METRE, -3.0, 0)]
units['kg_per_m3']   = [(UNIT_KIND_KILOGRAM, 1.0, 0), 
                    (UNIT_KIND_METRE, -3.0, 0)]
units['m3_per_skg']   = [(UNIT_KIND_METRE, 3.0, 0), 
                    (UNIT_KIND_KILOGRAM, -1.0, 0), (UNIT_KIND_SECOND, -1.0, 0)]

##########################################################################
# Parameters
########################################################################## 
pars = [# id, value, unit, constant            
            ('scale_f',   0.31,   'per_m3',    True),
            ('REF_P',     1.0,      'mM',   True),
            ('deficiency',  0,      '-',    True),
]
##########################################################################
# Assignments
########################################################################## 
assignments = [# id, assignment, unit       
               ]
##########################################################################
# Rules
########################################################################## 
rules = [# id, name, rule, unit
            ('c__scale', 'scale_f * Vol_cell', '-'),   
            ('e__gal_tot', 'e__gal + e__galM', 'mM'),
            ('c__gal_tot', 'c__gal + c__galM', 'mM'),
]
##########################################################################
# Compartments
########################################################################## 
comp_df = DataFrame(columns=['id', 'name', 'size', 'spatialDimensions', 'constant'],
                       data=[
                             ['c', 'cytosol', 1.0, 3, False],
                             ['m', 'membrane', 1.0, 2, False],
                             ['e', 'extracellular', 1.0, 3, False],
                             ['n', 'none', 1.0, 3, False],
                            ])
comp_df.set_index(comp_df.id, inplace=True)



def create_units(model):
    pass
    
    
def create_parameters(model):
    pass


def create_species(model):
    # Create metabolites and proteins
    pass
    

def create_compartments(model, comp_df):
    ''' Create compartments based on compartment information. '''
    for index, row in comp_df.iterrows():
        c = model.createCompartment()
        c.setId(row['id'])
        c.setName(row['name'])
        c.setSize(row['size'])
        c.setSpatialDimensions(row['spatialDimensions'])
        c.setConstant(row['constant'])

def create_species(model):
    ''' Creates species for metabolite and protein counts '''
    pass


if __name__ == "__main__":
    # Load matrix information
    matrix_dir = os.path.join(RESULTS_DIR, 'fba_matrices')
    
    s_fba_df = pd.read_csv(os.path.join(matrix_dir, 's_fba.csv'), sep="\t", keep_default_na=False)
    s_fba_df.set_index('sid', inplace=True)

    r_fba_df = pd.read_csv(os.path.join(matrix_dir, 'r_fba.csv'), sep="\t")
    r_fba_df.set_index('rid', inplace=True)

    e_df = pd.read_csv(os.path.join(matrix_dir, 'e_fba.csv'), sep="\t")
    e_df.set_index('eid', inplace=True)
    e_df    
    
    mat_stoichiometry = pd.read_csv(os.path.join(matrix_dir, 'fbaReactionStoichiometryMatrix.csv'), sep="\t")
    mat_stoichiometry.set_index(s_fba_df.index, inplace=True)        
        
    
    mat_catalysis = pd.read_csv(os.path.join(matrix_dir, 'fbaReactionCatalysisMatrix.csv'), sep="\t")
    mat_catalysis.set_index(r_fba_df.index, inplace=True)  
    
    
    # ----------------------------------------
    # New SBML model with FBC support
    sbmlns = SBMLNamespaces(3,1,"fbc",1)
    
    doc = SBMLDocument(sbmlns)
    model = doc.createModel()

    # history & creators
    from sbml.model_history import set_history_information
    set_history_information(model)
    
    # compartments
    create_compartments(model, comp_df)    
    
    # species
    # <metabolites>
    for index, row in s_fba_df.iterrows():
        # print index
        s = model.createSpecies()
        s.setId(index)
        s.setName(row['name'])
        s.setConstant(False)
        s.setBoundaryCondition(False)
        s.setCompartment(row['compartment'])
        s.setHasOnlySubstanceUnits(False)
        
    # <proteins>
    for index, row in e_df.iterrows():        
        # p = model.createParameter()
        p = model.createSpecies()
        p.setId(index)
        name = row['name']
        if not pd.isnull(name):
            p.setName(name)
        p.setConstant(False)
        p.setBoundaryCondition(False)
        # TODO: proper way to find location of reactions & proteins
        p.setCompartment('c') # this is just fix
        p.setHasOnlySubstanceUnits(False)

    # <reactions>
    tol = 1E-8
    for index, row in r_fba_df.iterrows():
        r = model.createReaction()
        r.setId(index)
        name = row['name']
        if not pd.isnull(name):
            p.setName(name)
        r.setFast(False)
        
        # set proteins as modifiers from catalysis matrix       
        row = mat_catalysis.ix[index]
        row = row[row>tol]
        for eid, value in row.iteritems():
            mod = r.createModifier()
            mod.setSpecies(eid)        
                
        # stoichiometry from stoichiometric matrix # [376x504]
        # find the non-zero elements in the reaction column 
        col = mat_stoichiometry[index]
        col = col[abs(col)>tol]
        # create species references depending on stoichiometry
        for sid, stoichiometry in col.iteritems():
            if stoichiometry < 0:
                rt = r.createReactant()
                rt.setSpecies(sid)
                rt.setStoichiometry(abs(stoichiometry))
                rt.setConstant(True)
            if stoichiometry > 0:
                pt = r.createProduct()
                pt.setSpecies(sid)
                pt.setStoichiometry(stoichiometry)
                pt.setConstant(True)
    
        # reversibility from the initial reaction bounds
        # TODO: calculate from flux bounds
        r.setReversible(False)
        
        # set local parameters for calculating flux bounds
        klaw = r.createKineticLaw()        
        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds', 'lb_fbaEnzymeBounds', 'ub_fbaEnzymeBounds'):
            lp = klaw.createLocalParameter()
            lp.setId(p_name)
            lp.setValue(r_fba_df[p_name][index])
    
    
    # FBC package
    # set reaction flux bounds
    mplugin = model.getPlugin("fbc");
 
    bound = mplugin.createFluxBound();
    bound.setId("bound1");
    bound.setReaction("J0");
    bound.setOperation("equal");
    bound.setValue(10);
 
    # <gene associations>
    mplugin.createGeneAssociation()
        
        
    # <objective function>
    objective = mplugin.createObjective();
    objective.setId("obj1");
    objective.setType("maximize");
    mplugin.setActiveObjectiveId("obj1");
    
    fluxObjective = objective.createFluxObjective();
    fluxObjective.setReaction("J8");
    fluxObjective.setCoefficient(1);

    # <chemical formulas => for balance>


    # write sbml    
    sbml_out = os.path.join(RESULTS_DIR, "Metabolism_matrices_{}_L3V1.xml".format(VERSION))
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_out)
    print sbml_out   
    from sbml_tools.checks import check_sbml
    check_sbml(sbml_out) 
    
    
    # TODO: gene associations
    
    
    
    
    # validator = SBMLValidator()
    # validator.validate(doc)
    
    
    