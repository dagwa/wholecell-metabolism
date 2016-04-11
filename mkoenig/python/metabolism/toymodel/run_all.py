"""
Create all files and run the simulations.
"""
from __future__ import print_function, division
from sbmlutils import comp
from simulator import simulate_toy_model
from simsettings import *
import model_factory
import comp_factory

import os

os.chdir(out_dir)  # set working dir to models

###########################################
# Create single models
###########################################
model_factory.create_fba(fba_file)
model_factory.create_ode_bounds(ode_bounds_file)
model_factory.create_ode_update(ode_update_file)
model_factory.create_ode_model(ode_model_file)

###########################################
# Create top level model
###########################################
comp_factory.create_top_level_model(top_level_file)
# flatten the combined model

comp.flattenSBMLFile(top_level_file, output_file=flattened_file)

###########################################
# Simulate the comp models
###########################################
simulate_toy_model(tend=50.0, step_size=0.1)
