# -*- coding: utf-8 -*-
"""
Adding history information to SBML file based on a creator data frame.

@author: 'mkoenig'
@date: '3/16/15'
"""
from libsbml import *

from pandas import DataFrame
from metabolism_settings import RESULTS_DIR, VERSION

#######################################################################
sbml_L3V1_history = os.path.join(RESULTS_DIR, "Metabolism_history_{}_L3V1.xml".format(VERSION))


def set_history_in_model(model, creators_df):
    #get model


def set_CV_term():
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
    print creators_df


    cvterms = [
        [MODEL_QUALIFIER, BQM_IS_DERIVED_FROM, "http://identifiers.org/doi/10.1016/j.cell.2012.05.044"],
        [BIOLOGICAL_QUALIFIER, BQB_HAS_TAXON, "http://identifiers.org/taxonomy/243273"],

        [MODEL_QUALIFIER, BQM_IS_DERIVED_FROM, "http://identifiers.org/doi/10.1016/j.cell.2012.05.044"],

    ]
    cvterms_df = DataFrame(data=cvterms, columns=['Qualifier', 'QualifierType', 'Resource')
    
    
    date = 
    model_id = 'WCM_3_
    model_name = 'Whole Cell 2015 - Metabolism'
    
    
    
    doc = SBMLDocument()
    model = 

