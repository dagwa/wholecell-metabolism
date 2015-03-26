'''
Write SBML for metabolism process directly from database information.

@author: Matthias Koeng
@date: 2015-03-25
'''
from libsbml import *
from sbml_tools.validator import SBMLValidator
from metabolism_settings import RESULTS_DIR, VERSION        

import pandas as pd
from pandas import DataFrame

##########################################################################
# Compartments
########################################################################## 
# compartment information is hard coded.
# The additional compartment 'none' is added to account for pseudo-metabolites
# used in the FBA.
comp_df = DataFrame(columns=['id', 'name', 'size', 'spatialDimensions', 'constant'],
                       data=[
                             ['c', 'cytosol', 1.0, 3, False],
                             ['m', 'membrane', 1.0, 2, False],
                             ['e', 'extracellular', 1.0, 3, False],
                             ['n', 'none', 1.0, 3, False],
                            ])
comp_df.set_index(comp_df.id, inplace=True)

    
def create_compartments(model, comp_df):
    ''' Create compartments based on compartment information. '''
    for index, row in comp_df.iterrows():
        c = model.createCompartment()
        c.setId(row['id'])
        c.setName(row['name'])
        c.setSize(row['size'])
        c.setSpatialDimensions(row['spatialDimensions'])
        c.setConstant(row['constant'])


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
    # SBML Creation
    ###########################################################################
    tol = 1E-12         # within this tolerance matrix elements are considered != 0    
                        # for instance in stoichiometric matrix
    
    # SBML model with FBC support    
    sbmlns = SBMLNamespaces(3,1,"fbc",1)
    doc = SBMLDocument(sbmlns)
    doc.setPackageRequired("fbc", False)
    model = doc.createModel()
    mplugin = model.getPlugin("fbc");           

    # history & creators
    from sbml.model_history import set_history_information
    set_history_information(model)
    
    # compartments
    create_compartments(model, comp_df)    
    
    # <metabolites>
    # Metabolites in the FBA problem (rows) are encoded as species
    for index, row in s_fba_df.iterrows():
        s = model.createSpecies()
        s.setId(index)
        s.setName(row['name'])
        s.setConstant(False)
        s.setBoundaryCondition(False)
        s.setCompartment(row['compartment'])
        s.setHasOnlySubstanceUnits(False)
        s.setInitialAmount(0)  # fix to get rid of warnings  
        
        # chemical formula and charge => for balance
        splugin = s.getPlugin("fbc");
        formula = row['formula']
        if not pd.isnull(formula):
            splugin.setChemicalFormula(formula)
        charge = row['charge']

        # string to int desaster due to NA name for sodium
        # this is ugly but works        
        if not pd.isnull(charge) and len(charge)!=0:            
            splugin.setCharge(int(float(charge)))
    
    # <proteins> [104]
    # Proteins (ProteinMonomer & ProteinComplex) are encoded as species.
    # The main reasons are:
    #   1. The protein count is changing during the simulation and input of the dynamical flux bound
    #      calculation.
    #   2. The proteins can be encoded as modifiers of the respective reactions.
    #      This provides clarity for the reaction <- protein <- gene information
    #      And provides important information for possible visualization.
    
    def create_protein_species(sid, name):
        s = model.createSpecies()
        s.setId(sid)
        # check name
        if not pd.isnull(name):
            s.setName(name)
        s.setConstant(False)
        s.setBoundaryCondition(False)
        # TODO: proper way to find location of reactions & proteins
        # Not important for simulation, only for visualization
        s.setCompartment('c')             # this is just fix
        s.setHasOnlySubstanceUnits(False) # ? 
        s.setInitialAmount(0)  # fix to get rid of warnings  
    
    # Handle all Enzymes
    for index, row in e_df.iterrows():
        # check if the protein is already a species (due to involvment in reaction)
        s = model.getSpecies(index)
        if (s is not None):
            print index, 'is already species.'
        else:
            # create species for protein
            create_protein_species(sid=index, name=row['name'])
        
    # Handle the special case of ACP (MG_287_MONOMER - acyl carrier protein)
    # Neither substrate nor enzyme (no part of FBA, but part of fluxbound calculation)
    create_protein_species(sid='MG_287_MONOMER', name="acyl carrier protein")

    # <reactions>
    # Reactions are all columns in the FBA stoichiometric matrix.
    # This includes some pseudo-reactions (internal & external exchange) which
    # are not represented in the knowledgbase.
    for index, row in r_fba_df.iterrows():
        r = model.createReaction()
        r.setId(index)
        name = row['name']
        if not pd.isnull(name):
            r.setName(name)
        r.setFast(False)
        
        # set proteins as modifiers from catalysis matrix       
        row = mat_catalysis.ix[index]
        row = row[row>tol]
        for eid, value in row.iteritems():
            # set protein as modifier
            mod = r.createModifier()
            mod.setSpecies(eid) 

            # gene associations
            gene_str = e_df['genes'][eid]
            genes = [g.strip() for g in gene_str.split(',')]
            genes_formula = '*'.join(genes)
            ga = mplugin.createGeneAssociation()
            ga.setId('ga__{}__{}'.format(index, eid))
            ga.setReaction(index)
            ass = Association_parseInfixAssociation(genes_formula)            
            ga.setAssociation(ass)

        # stoichiometry from stoichiometric matrix  # [376x504]
        # find non-zero elements in the reaction column 
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
    
        # <reversibility>
        # The reversibility can be calculated from the reaction bounds. In some
        # cases the reversibility is in backward direction. This would require the 
        # change of reactants and products for encoding. The SBML is strictly reproducing
        # the FBA problem, so that no reversibilities are defined in the SBML. 
        # Reversibility is functional in the FBA via the actual flux bounds.
        r.setReversible(True) # some are irreversible via Flux bounds in forward or backward direction
        
        # <fluxbounds>
        # parameters for dynamical calculation of flux bounds
        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds', 'lb_fbaEnzymeBounds', 'ub_fbaEnzymeBounds'):
            par = model.createParameter()
            par.setId('{}__{}'.format(index, p_name))
            par.setValue(r_fba_df[p_name][index])
            par.setConstant(True) 
        # The reaction flux bounds are set as hard upper and lower flux bounds
        # These are NOT the dynamical flux bounds.
        
        for p_name in ('lb_fbaReactionBounds', 'ub_fbaReactionBounds'):
            # "lessEqual", "greaterEqual", "equal"
            bound = mplugin.createFluxBound();            
            bound.setReaction(index);
            bound.setValue(r_fba_df[p_name][index])
            if p_name.startswith('lb'):
                bound.setId('lb__{}'.format(index))
                bound.setOperation("greaterEqual")
            if p_name.startswith('ub'):
                bound.setId('ub__{}'.format(index))
                bound.setOperation("lessEqual")                
        
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
    