"""
Settings for the toy model,
like file names and directories.
"""
from __future__ import print_function, division
import os
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')

# FBA submodel
fba_file = os.path.join(out_dir, 'toy_fba.xml')
# ODE submodels
ode_bounds_file = os.path.join(out_dir, 'toy_ode_bounds.xml')
ode_update_file = os.path.join(out_dir, 'toy_ode_update.xml')
ode_model_file = os.path.join(out_dir, 'toy_ode_model.xml')
# comp ODE submodels
comp_ode_file = os.path.join(out_dir, 'toy_comp_ode.xml')
# comp ODE & FBA model
comp_full_file = os.path.join(out_dir, 'toy_comp_full.xml')

