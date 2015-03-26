'''
Functions to work with matlab state dumps.
For instance get the variables from the state dump

@author: mkoenig
@date: 2015-03-25
'''

import scipy.io
import numpy as np

class AttrDict(dict):
    '''
        Access keys as attributes.
        http://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute-in-python
    '''
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def read_state(state_file):
    ''' Reading matlab state dump into AttrDict. 
        All fields are accessible as attributes.
    '''
    state_dict = scipy.io.loadmat(state_file)
    state = AttrDict(state_dict)
    return state

    
def print_state(state):
    ''' Print the content and dimensions of all state variables. '''
    for key, value in sorted(state.iteritems()):
        if isinstance(value, np.ndarray):
            print key, value.shape 
        else:
            print key