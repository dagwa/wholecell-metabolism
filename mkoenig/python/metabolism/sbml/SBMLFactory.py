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

from metabolism_settings import RESULTS_DIR, VERSION


def createSBML(model_id):
    '''
    Create the SBML from scratch using the database information.
    '''
    sbml_out = os.path.join(RESULTS_DIR, model_id)
    print '***', model_id, '***'
    print model_id
    
    try:
        document = SBMLDocument(3, 1)
    except ValueError:
        raise SystemExit('Could not create SBMLDocumention object')
    
    model = document.createModel()
    check(model, 'create model')
    model.setId(model_id)
    model.setName(model_id)
    
    # create compartments
    
    # create species
    # necessary to add proteins for the reactions, i.e. all the proteins used in the model
    
    return sbml_out

if __name__ == "__main__":
    mid = "Metabolism_annotated_{}.xml".format(VERSION)
    sbml_out = createSBML(mid)
    validator = SBMLValidator()
    validator.validate(sbml_out)
