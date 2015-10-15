"""
Trying to set selection and and run simulation.
@author: Matthias Koenig
"""
import roadrunner
print roadrunner.__version__
# 1.4.0; Compiler: gcc 4.8.4, C++ version: 199711; JIT Compiler: LLVM-3.4; Date: Oct 14 2015, 00:00:59

# Load model results in warnings about fbc v2
rr = roadrunner.RoadRunner("toy_fba.xml")

# create a selection array
sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr.model.getFloatingSpeciesIds()]
print sel

# new way is handled via timeCourseSelections
rr.timeCourseSelections = sel
result = rr.simulate(start=0, end=10, steps=5)
print result
