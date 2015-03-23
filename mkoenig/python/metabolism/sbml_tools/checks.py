# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 08:33:17 2015

@author: mkoenig
"""
import time
from libsbml import *

def check_sbml(filename): 
    current = time.clock()
    doc = readSBML(filename)
    errors = doc.getNumErrors()
    
    print
    print(" filename: " + filename)
    print(" file size: " + str(os.stat(filename).st_size))
    print(" read time (ms): " + str(time.clock() - current))
    print(" validation error(s): " + str(errors))
    print
    doc.printErrors()
    return errors
    
def check(value, message):
    if value == None:
        print('LibSBML returned a null value trying to ' + message + '.')
        sys.exit(1)
    elif type(value) is int:
        if value == LIBSBML_OPERATION_SUCCESS:
            return
        else:
            print('Error encountered trying to ' + message + '.')
            print('LibSBML returned error code ' + str(value) + ': "' \
                + OperationReturnValue_toString(value).strip() + '"')
            print('Exiting.')
            sys.exit(1)
    else:
        return
    