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

    # handle the sodium NA id (not parsing as NaN)    
    s_fba_df = pd.read_csv(os.path.join(matrix_dir, 's_fba.csv'), sep="\t", 
                           keep_default_na=False, na_values=('nan'))
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
    
    
    ###########################################################################
    # New SBML model with FBC support
    sbmlns = SBMLNamespaces(3,1,"fbc",1)
    
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("fbc", False)
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
        
        # chemical formula and charge => for balance
        splugin = s.getPlugin("fbc");
        formula = row['formula']
        if not pd.isnull(formula):
            splugin.setChemicalFormula(formula)
        charge = row['charge']

        # string to int desaster due to NA        
        if not pd.isnull(charge) and len(charge)!=0:            
            splugin.setCharge(int(float(charge)))
    
    # <proteins>
    for index, row in e_df.iterrows():        
        # p = model.createParameter()
        # check if the protein is already a species (due to involvment in reaction)
        s = model.getSpecies(index)
        if (s is not None):
            print index, 'is already species.'
            continue
    
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

    # FBC support
    mplugin = model.getPlugin("fbc");
    
    # <reactions>
    tol = 1E-12
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

            # gene associations
            gene_str = e_df['genes'][eid]
            genes = [g.strip() for g in gene_str.split(',')]
            genes_formula = '*'.join(genes)
            
            # <gene associations>
            ga = mplugin.createGeneAssociation()
            ga.setId('ga__{}__{}'.format(index, eid))
            ga.setReaction(index)
            
            # ast_node = parseL3Formula('*'.join(genes))
            ass = Association_parseInfixAssociation(genes_formula)            
            # ass = Association(ast_node)
            ga.setAssociation(ass)
            # "GENE_ASSOCIATION: MG_271 and MG_272 and MG_273 and MG_274"    
            infix = ass.toInfix()
            print infix     
            ass = ga.getAssociation()
            infix = ass.toInfix()
            print infix
                

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
        # Manage if reversible in backward direction
        r.setReversible(False)
        
        # not possible to set local parameters (math is required)
        # prefix with reaction to get the actual parameters
        # klaw = r.createKineticLaw()        
        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds', 'lb_fbaEnzymeBounds', 'ub_fbaEnzymeBounds'):
            # lp = klaw.createLocalParameter()
            # lp.setId(p_name)
            # lp.setValue(r_fba_df[p_name][index])
            par = model.createParameter()
            par.setId('{}__{}'.format(index, p_name))
            par.setValue(r_fba_df[p_name][index])
            par.setConstant(True)
    
         
    # bound = mplugin.createFluxBound();
    # bound.setId("bound1");
    # bound.setReaction("J0");
    # bound.setOperation("equal");
    # bound.setValue(10);
 
        
    # <objective function>
    objective = mplugin.createObjective();
    objective.setId("growth")
    objective.setType("maximize");
    mplugin.setActiveObjectiveId("growth");
    
    for index, row in r_fba_df.iterrows():
        coeff = row['fbaObjective']
        if (not pd.isnull(coeff) and abs(coeff)>tol):
            fluxObjective = objective.createFluxObjective();
            fluxObjective.setReaction(index);
            fluxObjective.setCoefficient(coeff);

    # write sbml    
    sbml_out = os.path.join(RESULTS_DIR, "Metabolism_matrices_{}_L3V1.xml".format(VERSION))
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_out)
    print sbml_out   
    from sbml_tools.checks import check_sbml
    check_sbml(sbml_out) 
    
    # ---------
    
    validator = SBMLValidator(False)
    print validator.validate(sbml_out)
    
    from libsbml import ConversionProperties, LIBSBML_OPERATION_SUCCESS
    conversion_properties = ConversionProperties()
    conversion_properties.addOption("convert fbc to cobra", True, "Convert FBC model to Cobra model")
    result = doc.convert(conversion_properties)
    if result != LIBSBML_OPERATION_SUCCESS:
        raise(Exception("Conversion of SBML+fbc to COBRA failed"))    
    sbml_cobra = os.path.join(RESULTS_DIR, "Metabolism_matrices_cobra_{}_L3V1.xml".format(VERSION))
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_cobra)
    
    