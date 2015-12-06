"""
Create all files and run the simulations.
"""

###########################################
# Create single models
###########################################
from settings import *
import model_factory

model_factory.create_fba(fba_file)
model_factory.create_ode_bounds(ode_bounds_file)
model_factory.create_ode_update(ode_update_file, fba_file)
model_factory.create_ode_model(ode_model_file)

# store in local database
import multiscale.multiscalesite.simapp.db.api as db_api
db_api.create_model(fba_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_bounds_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_update_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_model_file, model_format=db_api.CompModelFormat.SBML)


###########################################
# Create comp models
###########################################
import comp_factory
comp_factory.create_comp_ode_model(comp_ode_file)
comp_factory.create_comp_full_model(comp_full_file)

###########################################
# Simulate the comp models
###########################################
# Simulate via manual connection approach
# comp_ode & fba
import simulator

df1 = simulator.simulate_manual(fba_sbml=fba_file, comp_ode_sbml=comp_ode_file,
                          tend=50.0, step_size=0.1, debug=False)
# df2 = simulate(tend=10.0, step_size=None, debug=False)

df1.plot(x='time', y=['submodel_update__R1',
                      'submodel_update__R2',
                      'submodel_update__R3',
                      'submodel_model__R4'])
df1.plot(x='time', y=['[submodel_update__A]',
                      '[submodel_update__B1]',
                      '[submodel_update__B2]',
                      '[C]',
                      '[submodel_model__D]'])