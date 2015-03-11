'''
Created on Mar 11, 2015

@author: mkoenig
'''
import os
from metabolism_settings import VERSION, RESULTS_DIR 
sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_{}_L3V1.xml".format(VERSION))

import cobra

model = cobra.io.read_sbml_model(sbml)

# Perform FBA
model.optimize()

# Solution status
model.solution.status

# objective function
{reaction: reaction.objective_coefficient for reaction in model.reactions
if reaction.objective_coefficient != 0}


# Print reactions & flux bounds
print '*** Reactions ***'
for r in model.reactions:
    print r, r.upper_bound, r.lower_bound
print '*'*80
