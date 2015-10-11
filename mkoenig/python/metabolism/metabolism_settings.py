"""
Global settings for the python scripts to write and read data.

@author: Matthias Koenig
@date: 2015-10-11
"""
import os
VERSION = '8'
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
if True:
    DATA_DIR = os.path.join(BASE_DIR, "../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../results")
else:
    DATA_DIR = os.path.join(BASE_DIR, "../../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../../results")
