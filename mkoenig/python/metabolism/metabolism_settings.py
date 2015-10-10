'''
Created on Mar 10, 2015

@author: Matthias Koenig
@date: global settings for the python scripts to write to correct folders
'''
import os
VERSION = '7'
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
if (False):
    DATA_DIR = os.path.join(BASE_DIR, "../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../results")
else:
    DATA_DIR = os.path.join(BASE_DIR, "../../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../../results")
