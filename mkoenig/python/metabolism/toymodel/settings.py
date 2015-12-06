"""
Settings for the toy model,
like file names and directories.
"""

import os
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')

fba_file = os.path.join(out_dir, 'toy_fba.xml')
ode_bounds_file = os.path.join(out_dir, 'toy_ode_bounds.xml')
ode_update_file = os.path.join(out_dir, 'toy_ode_update.xml')
ode_model_file = os.path.join(out_dir, 'toy_ode_model.xml')
comp_file = os.path.join(out_dir, 'toy_comp.xml')
