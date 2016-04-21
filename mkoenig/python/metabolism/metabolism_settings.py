"""
Global settings for the python scripts to write and read data.
"""
import os
VERSION = '09'
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
if True:
    DATA_DIR = os.path.join(BASE_DIR, "../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../results")
else:
    DATA_DIR = os.path.join(BASE_DIR, "../../data") 
    RESULTS_DIR = os.path.join(BASE_DIR, "../../results")

sbml_L3V1_history = os.path.join(RESULTS_DIR, "Metabolism_history_{}_L3V1.xml".format(VERSION))