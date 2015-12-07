"""
Create all files and run the simulations.
"""
from settings import *
import model_factory
import comp_factory
import simulator
import multiscale.multiscalesite.simapp.db.api as db_api

###########################################
# Create single models
###########################################
model_factory.create_fba(fba_file)
model_factory.create_ode_bounds(ode_bounds_file)
model_factory.create_ode_update(ode_update_file, fba_file)
model_factory.create_ode_model(ode_model_file)

# store in local database & report
db_api.create_model(fba_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_bounds_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_update_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(ode_model_file, model_format=db_api.CompModelFormat.SBML)


###########################################
# Create comp models
###########################################
comp_factory.create_comp_ode_model(comp_ode_file)
comp_factory.create_comp_full_model(comp_full_file)

db_api.create_model(comp_ode_file, model_format=db_api.CompModelFormat.SBML)
db_api.create_model(comp_full_file, model_format=db_api.CompModelFormat.SBML)


###########################################
# Simulate the comp models
###########################################
# comp_ode & fba (Simulate via manual connection approach)
df_test = simulator.simulate_manual(fba_sbml=fba_file, comp_ode_sbml=comp_ode_file,
                                    tend=50.0, step_size=0.1, debug=False)

# simulate the complete model
df = simulator.simulate(mixed_sbml=comp_full_file, tend=50.0, step_size=0.1)
df.plot(x='time', y=['fba__R1', 'fba__R2', 'fba__R3', 'model__R4'])
df.plot(x='time', y=['[update__A]',
                     '[update__B1]',
                     '[update__B2]',
                      '[C]',
                      '[model__D]'])
