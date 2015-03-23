# -*- coding: utf-8 -*-
"""
Adding history information to SBML file based on a creator data frame.

@author: 'mkoenig'
@date: '2015-03-23'
"""
from libsbml import *

from pandas import DataFrame
from metabolism_settings import RESULTS_DIR, VERSION

from sbml_tools.checks import check 
from sbml.annotation import create_meta_id

###############################################################################
sbml_L3V1_history = os.path.join(RESULTS_DIR, "Metabolism_history_{}_L3V1.xml".format(VERSION))
creators = [
        ['Bergmann', 'Frank', 'fbergman@caltech.edu', 'Caltech'],
        ['Koenig', 'Matthias', 'konigmatt@googlemail.com', 'Charite Berlin'],
        ['Smallbone', 'Kieran', 'kieran.smallbone@manchester.ac.uk', 'University of Manchester'],
        ['Tokic', 'Milenko', 'milenko.tokic@epfl.ch', 'EPFL'],
        ['Costa', 'Rafael', '‎rcosta@kdbio.inesc-id.pt', 'University of Lisbon'],
        ['Baghalian', 'Kambiz', 'kambiz.baghalian@plants.ox.ac.uk', 'University of Oxford'],
]
creators_df = DataFrame(data=creators, columns=['FamilyName', 'GivenName', 'Email', 'Organization'])
    
cvterms = [
        [MODEL_QUALIFIER, BQM_IS_DERIVED_FROM, "http://identifiers.org/doi/10.1016/j.cell.2012.05.044"],
        [MODEL_QUALIFIER, BQM_IS_DESCRIBED_BY, "http://identifiers.org/pubmed/22817898"], 
        [MODEL_QUALIFIER, BQM_IS,  'http://identifiers.org/mamo/MAMO_0000040'], # metabolic network
        [MODEL_QUALIFIER, BQM_IS,  'http://identifiers.org/mamo/MAMO_0000009'], # constraint-based model
        [BIOLOGICAL_QUALIFIER, BQB_HAS_TAXON, "http://identifiers.org/taxonomy/243273"], # Mycoplasma genitalium G37
]
cvterms_df = DataFrame(data=cvterms, columns=['Qualifier', 'QualifierType', 'Resource'])
###############################################################################


def date_now():
    import datetime 
    time = datetime.datetime.now()
    timestr = time.strftime('%Y-%m-%dT%H:%M:%S')
    return Date(timestr);

def create_history(creators_df):
    h = ModelHistory();
 
    # add all creators
    for k in range(0, len(creators_df)):
        c = ModelCreator()
        c.setFamilyName(creators_df.FamilyName[k])
        c.setGivenName(creators_df.GivenName[k])
        c.setEmail(creators_df.Email[k])
        c.setOrganization(creators_df.Organization[k])
        check(h.addCreator(c), 'add creator');

    # create time is now            
    date = date_now()
    check(h.setCreatedDate(date), 'set creation date')
    check(h.setModifiedDate(date), 'set creation date')
    return h

def set_model_id_name(model, model_id='WCM_3_10', 
                      model_name = 'Whole Cell 2015 - Metabolism'):
    # model name and id for process     
    model.setId(model_id)
    model.setName(model_name)

def set_model_history(model, creators_df):
    if not model.isSetMetaId():
        print 'setting meta id for model'
        model.setMetaId(create_meta_id(model.getId()))      
    
    # set history
    h = create_history(creators_df)
    check(model.setModelHistory(h), 'set model history')
    
def set_cv_terms(model, cvterms_df):
    ''' Set model cv terms from DataFrame. '''
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

def set_history_information(model):
    set_model_id_name(model)
    set_model_history(model, creators_df)
    set_cv_terms(model, cvterms_df)


if __name__ == "__main__":

    doc = SBMLDocument(3,1)
    model = doc.createModel()
    
    set_history_information(model)
    
    # write sbml    
    sbml_out = os.path.join(RESULTS_DIR, 'history_test.xml')
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_out)
    print sbml_out    
    