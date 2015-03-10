# -*- coding: utf-8 -*-
"""
Annotation of Metabolism.sbml with information from the 
publication supplement.

Annotation are divided into
- CV terms (annotated via CV terms)
- annotation strings (not implemented yet)

The CV terms are defined for metabolites and reactions in m_cvdf & r_cvdf.

-----------------------------------
ChangeLog:
-----------------------------------
v4
- xml file ending
- BQB_QUALIFIER better specified ( BQB_IS, BQB_IS_PROPERTY_OF,
        BQB_IS_VERSION_OF)
-----------------------------------

@author: Matthias Koenig
@date: 2015-03-08
"""
from libsbml import *
from sbml_tools.checks import check_sbml, check
from pandas import DataFrame
import pandas as pd
from metabolism_settings import RESULTS_DIR, DATA_DIR, VERSION
from atk import Document

#######################################################################
sbml_raw = os.path.join(DATA_DIR, 'Metabolism.sbml')
sbml_out = os.path.join(RESULTS_DIR, "Metabolism_annotated_{}.xml".format(VERSION))
csv_metabolites = os.path.join(DATA_DIR, "Table_S3G_metabolites.csv")
csv_reactions = os.path.join(DATA_DIR, "Table_S3O_reactions.csv")
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
    ['BioCyc', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/biocyc/'],
    # TODO: no resource found, official link for BiGG resources ??
    ['BiGG', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://bigg.ucsd.edu/bigg/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000013
    ['KEGG', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/kegg.compound/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000237
    ['CAS', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/cas/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000034
    # TODO: used pubchem.compound, there is also pubchem.substance ?? 
    ['PubChem', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/pubchem.compound/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000066
    ['3DMET', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/3dmet/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000002
    ['ChEBI', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/chebi/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000113
    ['PDB-CCD', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/pdb-ccd/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000271
    ['KNApSAcK', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/knapsack/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000052
    ['LIPID MAPS', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/lipidmaps/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000115
    ['LipidBank', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/lipidbank/'],
]
m_cvdf = cvdf_from_resource_data(m_cv_data)
# print m_cvdf
#######################################################################
# Reaction Resources
r_cv_data = [
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000194
    ['BioCyc', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://identifiers.org/biocyc/'],
    # TODO: no resource found, official link for BiGG resources ??
    ['BiGG', BIOLOGICAL_QUALIFIER, BQB_IS, 'http://bigg.ucsd.edu/bigg/'],
    # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000086
    # TODO: multiple Sabio entries (Forward & Backward)
    ['SABIO-RK ID (Forward)', BIOLOGICAL_QUALIFIER, BQB_IS_PROPERTY_OF, 'http://identifiers.org/sabiork.kineticrecord/'],
    ['SABIO-RK ID (Backward)', BIOLOGICAL_QUALIFIER, BQB_IS_PROPERTY_OF, 'http://identifiers.org/sabiork.kineticrecord/'],
]
r_cvdf = cvdf_from_resource_data(r_cv_data)
# print r_cvdf
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

def annotate_model_cv(m):
    '''
    Annotate the model with given annotation data frames
    for metabolites and reactions.
    Uses the global cvdf datasets and the parsed table information.
    '''
    # read original model    
    annotate_objects(m.getListOfSpecies(), m_df, m_cvdf, otype='SPECIES')
    annotate_objects(m.getListOfReactions(), r_df, r_cvdf, otype='REACTION')


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
                

def annotate_model_sbo(m):
    from public.models import Entry
    from public.models import Reaction as DBReaction
    from django.core.exceptions import ObjectDoesNotExist
    print '* Annotate SBO *'
    
    # compartments
    for c in m.getListOfCompartments():
        sbo_id = "SBO:0000290"   # SBO:0000290 - physical compartment
        check(c.setSBOTerm(sbo_id), 'Set SBO')
    
    # species
    sbo_dict = {
            'Metabolite' : 'SBO:0000247',     # simple chemical
            'Gene' : 'SBO:0000243',           # gene
            'Protein' : 'SBO:0000245',        # macromolecule
            'ProteinMonomer' : 'SBO:0000245', # macromolecule
            'ProteinComplex' : 'SBO:0000297', # protein complex
            'Stimulus' : 'SBO:0000170'        # stimulation.
    }
    for s in m.getListOfSpecies():
        # lookup the Entry in the database
        sid = s.getId()
        wid = cid_from_sid(sid)  # get wid from sid
        try:
            e = Entry.objects.get(wid=wid)
            m_type = e.model_type
            sbo_id = sbo_dict[m_type]
            check(s.setSBOTerm(sbo_id), 'Set SBO')
        except ObjectDoesNotExist:
            print 'Warning - Entry not existing in DB, no SBO', sid, wid
    
    # reactions
    sbo_dict = {
            'TransportReaction' : 'SBO:0000185',  # transport reaction (find via compartments)
            'Reaction' : 'SBO:0000176',           # biochemical reaction
    }
    for r in m.getListOfReactions():
        # lookup the Entry in the database
        sid = r.getId()
        wid = cid_from_rid(sid)  # get wid from sid
        try:
            e = Entry.objects.get(wid=wid)
            m_type = e.model_type
            if m_type == 'Reaction':
                # check if multiple compartments, than transporter
                reaction = DBReaction.objects.get(wid=wid)
                comps = set([c.compartment for c in reaction.stoichiometry.all()])
                if len(comps)>1:
                    m_type = 'TransportReaction'
            sbo_id = sbo_dict[m_type]
            check(r.setSBOTerm(sbo_id), 'Set SBO')
        except ObjectDoesNotExist:
            print 'Warning - Entry not existing in DB, no SBO', sid, wid
        
        
        


def create_sbo_dict():
    """
    Necessary to define a dictionary of SBO.
    """
    pass
    '''
    sbo_dict = {
                # Compartment
                'Compartment' : 'SBO:0000290',   # physical compartment
                
                # Species
                'Metabolite' : 'SBO:0000247',    # simple chemical
                'Gene': 'SBO:0000243',
                'Protein':'SBO:0000297',         # protein (SBO:0000420 -> macromolecule)
                'ProteinComplex': 'SBO:0000418', # complex
                
                # Reactions
                'Reaction': '',
                # SBO:0000375 ! process
                # SBO:0000397 - omitted process
                # SBO:0000396 ! uncertain process.
                # SBO:0000177 ! non-covalent binding.
                # SBO:0000180 ! dissociation.
                # SBO:0000394 ! consumption.
                # SBO:0000393 ! production.
                
                # is spontaneuos
                SBO:0000185 - transport reaction (find via compartments)
                SBO:0000176 - biochemical reaction
                # SBO:0000396 - uncertain process
                # SBO:0000180 - dissociation
                
                # ? Stimulus
                # SBO:0000170 ! stimulation.
                # SBO:0000172 ! catalysis.
                
                # SpeciesReference
                SBO:0000019 - modifier
                # SBO:0000594 - neutral participant
                SBO:0000011 - product
                # SBO:0000598 - promoter
                SBO:0000010 - reactant
    ''' 
    


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
    
    # read model
    doc = readSBML(sbml_raw)
    doc.setLevelAndVersion(2,4,False, True)
    # write id and name
    m = doc.getModel()
    mid = 'Metabolism_annotated_{}'.format(VERSION)
    m.setId(mid)
    m.setName(mid)
    # annotate
    
    # annotate_model_cv(m)
    annotate_model_sbo(m)

    # convert to 3.1
    convert = False
    if convert:
        props = ConversionProperties()
        props.addOption("convert cobra", True, "Convert Cobra model")
        check(doc.convert(props), 'Convert COBRA')

    # save
    writeSBMLToFile(doc, sbml_out)
    check_sbml(sbml_out)
    print sbml_out
    
    
    
    