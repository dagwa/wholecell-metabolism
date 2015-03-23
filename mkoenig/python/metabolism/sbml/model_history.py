# -*- coding: utf-8 -*-
"""
Adding history information to SBML file based on a creator data frame.

@author: 'mkoenig'
@date: '3/16/15'
"""
from libsbml import *

from pandas import DataFrame
from metabolism_settings import RESULTS_DIR, VERSION

from sbml_tools.checks import check 
from sbml.annotation import create_meta_id

#######################################################################
sbml_L3V1_history = os.path.join(RESULTS_DIR, "Metabolism_history_{}_L3V1.xml".format(VERSION))


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
        
        c = h.getCreator(k)
        print c
 
    # TODO: proper date now
    # date = datetime
    date = Date("1999-11-13T06:54:32");
    check(h.setCreatedDate(date), 'set creation date')
    check(h.setModifiedDate(date), 'set creation date')
 
    return h
    

def set_CV_terms(model, cvterms_df):
    #annotate model with reference to original model
    cv = CVTerm()
    cv.setQualifierType(MODEL_QUALIFIER)
    cv.setModelQualifierType(BQM_IS_DERIVED_FROM)
    cv.addResource("http://identifiers.org/doi/10.1016/j.cell.2012.05.044")
    model.addCVTerm(cv)

    #annotate model with taxonomic information
    cv=CVTerm()
    cv.setQualifierType(BIOLOGICAL_QUALIFIER)
    cv.setBiologicalQualifierType(BQB_HAS_TAXON)
    cv.addResource("http://identifiers.org/taxonomy/243273")
    model.addCVTerm(cv)


if __name__ == "__main__":

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

    # model name and id
    model_id = 'WCM_3_10' # Metabolism
    model_name = 'Whole Cell 2015 - Metabolism'
    
    # create test model
    doc = SBMLDocument(3,1)
    model = doc.createModel()
    model.setId(model_id)
    model.setName(model_name)
    model.setMetaId(create_meta_id(model.getId()))
    
    # set history
    h = create_history(creators_df)
    print h
   
    
    print h.getListCreators()
    for c in h.getListCreators():
        print c
    
    for k in range(h.getListCreators().getSize()):
        c = h.getCreator(k)
        print c
    
    check(model.setModelHistory(h), 'set model history');    

    # set cv terms


    # write sbml    
    sbml_out = os.path.join(RESULTS_DIR, 'history_test.xml')
    writer = SBMLWriter()
    writer.writeSBML(doc, sbml_out)
    print sbml_out    
    