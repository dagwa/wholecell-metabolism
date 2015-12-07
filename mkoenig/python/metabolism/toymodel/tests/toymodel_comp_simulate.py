"""
Try to run a single step simulation.
@author: Matthias Koenig
"""
import roadrunner
print roadrunner.__version__
# 1.4.0; Compiler: gcc 4.8.4, C++ version: 199711; JIT Compiler: LLVM-3.4; Date: Oct 14 2015, 00:00:59

# Load model results in warnings about fbc v2
from toymodel.settings import comp_ode_file, comp_full_file
# rr = roadrunner.RoadRunner(comp_ode_file)
rr = roadrunner.RoadRunner(comp_full_file)

# boundary and floating species in selection

sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr.model.getFloatingSpeciesIds()]
rr.timeCourseSelections = sel

# simulate multiple steps
result = rr.simulate(start=0, end=10, steps=100)
print result
