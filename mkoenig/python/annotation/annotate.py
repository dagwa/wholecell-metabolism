# -*- coding: utf-8 -*-
"""
Annotation of Metabolism.sbml with information from the 
publication supplement.

Annotation are divided into
- CV terms (annotated via CV terms)
- annotation strings (not implemented yet)

The CV terms are defined for metabolites and reactions in m_cvdf & r_cvdf.

@author: Matthias Koenig
@date: 2015-03-08
"""
from libsbml import *
from annotate_tools import check_sbml, check
from pandas import DataFrame
import pandas as pd

#######################################################################
DATA_DIR = "../../data"
RESULTS_DIR = "../../results"
VERSION = 1
sbml_raw = "{}/Metabolism.sbml".format(DATA_DIR)
sbml_out = "{}/Metabolism_annotated_{}.sbml".format(RESULTS_DIR, VERSION)
csv_metabolites = "{}/Table_S3G_metabolites.csv".format(DATA_DIR)
csv_reactions = "{}/Table_S3O_reactions.csv".format(DATA_DIR)
#######################################################################

def cvdf_from_resource_data(r_data):
    cv_df = DataFrame(columns=('ID', 'BQB', 'Qualifier', 'URI'))
    # TODO: fix in place elongation of df
    for k, data in enumerate(r_data):
        cv_df.loc[k] = data
    cv_df = cv_df.set_index(cv_df.ID)
    return cv_df

# Metabolite Resources
m_cv_data = [
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000194
    ['BioCyc', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/biocyc/'],
    # TODO: no resource found, official link for BiGG resources ??
    ['BiGG', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://bigg.ucsd.edu/bigg/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000013
    ['KEGG', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/kegg.compound/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000237
    ['CAS', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/cas/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000034
    # TODO: used pubchem.compound, there is also pubchem.substance ?? 
    ['PubChem', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/pubchem.compound/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000066
    ['3DMET', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/3dmet/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000002
    ['ChEBI', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/chebi/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000113
    ['PDB-CCD', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/pdb-ccd/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000271
    ['KNApSAcK', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/knapsack/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000052
    ['LIPID MAPS', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/lipidmaps/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000115
    ['LipidBank', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/lipidbank/'],
]
m_cvdf = cvdf_from_resource_data(m_cv_data)
print m_cvdf
#######################################################################
# Reaction Resources
r_cv_data = [
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000194
    ['BioCyc', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/biocyc/'],
    # TODO: no resource found, official link for BiGG resources ??
    ['BiGG', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://bigg.ucsd.edu/bigg/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000086
    # TODO: multiple Sabio entries (Forward & Backward)
    ['SABIO-RK ID (Forward)', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/sabiork.kineticrecord/'],
    ['SABIO-RK ID (Backward)', BIOLOGICAL_QUALIFIER, BQB_IS_VERSION_OF, 'http://identifiers.org/sabiork.kineticrecord/'],
]
r_cvdf = cvdf_from_resource_data(r_cv_data)
print r_cvdf
#######################################################################

def cid_from_sid(sid):
    '''
    Create the compound id in the supplement tables from ids in SBML.
    The id lookup for species is without M_ and _compartment.
    For instance M_A23CMP_c -> A23CMP in the annotation table.
    '''
    tokens = sid.split('_')
    return '_'.join(tokens[1:(len(tokens)-1)])

def cid_from_rid(sid):
    '''
    Get compound id for lookup from SBML reaction id.
    '''
    tokens = sid.split('_')
    return '_'.join(tokens[1:len(tokens)])

def annotate_model(filename):
    '''
    Annotate the model with given annotation data frames
    for metabolites and reactions.
    Uses the global cvdf datasets and the parsed table information.
    '''
    # read original model
    doc = readSBML(filename)
    m = doc.getModel()
    
    mid = 'Metabolism_annotated_{}'.format(VERSION)
    m.setId(mid)
    m.setName(mid)
    
    annotate_objects(m.getListOfSpecies(), m_df, m_cvdf, otype='SPECIES')
    annotate_objects(m.getListOfReactions(), r_df, r_cvdf, otype='REACTION')

    # save
    writeSBMLToFile(doc, sbml_out)

def annotate_objects(objects, o_df, o_cvdf, otype):
        # Metabolite annotation
    # no fancy things here: just take the ids of the species
    # and look up existing annotations. Than write the CV terms
    # from the defined list above.
    
    for s in objects:        
        sid = s.getId()
        
        if otype == 'SPECIES':
            cid = cid_from_sid(sid)
        elif otype == 'REACTION':
            cid = cid_from_rid(sid)
        print '***', sid, '->', cid, '***'
        
        # WTF - not working without meta id, but no proper warning
        # TODO: how to properly generate meta ids
        s.setMetaId('meta_{}'.format(sid))
        
        # get the annotation info & create all CV terms
        # check if in index
        if not cid in o_df.index:
            pass
            # print cid, 'not in annotation data'
        else:
            
            for cv_type in o_cvdf.ID:   
                cv_id = o_df.loc[cid, cv_type]

                # check the NaN in table
                if pd.isnull(cv_id):                    
                    # print cid, 'no id available'
                    continue
                
                # create cv term                                
                qt = int(o_cvdf.loc[cv_type, 'BQB'])
                bqt = int(o_cvdf.loc[cv_type, 'Qualifier'])
                uri = '{}{}'.format(o_cvdf.loc[cv_type, 'URI'], cv_id)
                print qt, bqt, uri
                
                cv = CVTerm();
                check(cv.setQualifierType(qt), 'setQualifier')
                check(cv.setBiologicalQualifierType(bqt), 'setQualifierType')
                check(cv.addResource(uri), 'setURI')
                check(s.addCVTerm(cv), 'addCVTerm');
                
                

if __name__ == "__main__":
    check_sbml(sbml_raw)
    
    # Load annotation data & index with ID for O(1) lookup    
    m_df = pd.io.parsers.read_csv(csv_metabolites, sep="\t")
    m_df = m_df.set_index(m_df.ID)
    # print m_df.head()
    # m_df.ix['A23CMP']
    
    r_df = pd.io.parsers.read_csv(csv_reactions, sep="\t")
    r_df = r_df.set_index(r_df.ID)
    # print r_df.head()
    
    annotate_model(sbml_raw)    
    check_sbml(sbml_raw)
    