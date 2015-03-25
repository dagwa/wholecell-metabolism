'''
Functions to work with matlab state dumps.
For instance get the variables from the state dump

@author: mkoenig
@date: 2015-03-25
'''

import scipy.io
import numpy as np

def read_state(state_file):
    ''' Reading matlab state dump and providing the names as fields. '''
    state = scipy.io.loadmat(state_file)
    return state

    
def print_state(state):
    ''' Print the content and dimensions of all state variables. '''
    for key, value in sorted(state.iteritems()):
        if isinstance(value, np.ndarray):
            print key, value.shape 
        else:
            print key