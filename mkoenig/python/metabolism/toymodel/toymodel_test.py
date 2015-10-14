"""
This changed the selection before, how is this possible in latest version of libRoadRunner ??

"""

import roadrunner
print roadrunner.__version__
# 1.4.0; Compiler: gcc 4.8.4, C++ version: 199711; JIT Compiler: LLVM-3.4; Date: Oct 14 2015, 00:00:59

# Load model results in warnings about fbc v2
rr = roadrunner.RoadRunner("toy_ode_update.xml")

# create a selection array
sel = ['time'] \
        + ["".join(["[", item, "]"]) for item in rr.model.getBoundarySpeciesIds()] \
        + ["".join(["[", item, "]"]) for item in rr.model.getFloatingSpeciesIds()]
print sel

# set the selection !!! NOT WORKING ANY MORE. Why is it not possible ????
rr.selections = sel
result = rr.simulate(start=0, end=10, steps=5)
print result
# wrong number of points ! this should be 6 points, and stop at 10
# but is 51 points and stops at 5 !!!! (SOMETHING IS REALLY WRONG HERE)
"""
    time, [B1], [B2]
 [[    0,    0,    0],
  [  0.1,    0,    0],
  [  0.2,    0,    0],
  [  0.3,    0,    0],
  [  0.4,    0,    0],
  [  0.5,    0,    0],
  [  0.6,    0,    0],
  [  0.7,    0,    0],
  [  0.8,    0,    0],
  [  0.9,    0,    0],
  [    1,    0,    0],
  [  1.1,    0,    0],
  [  1.2,    0,    0],
  [  1.3,    0,    0],
  [  1.4,    0,    0],
  [  1.5,    0,    0],
  [  1.6,    0,    0],
  [  1.7,    0,    0],
  [  1.8,    0,    0],
  [  1.9,    0,    0],
  [    2,    0,    0],
  [  2.1,    0,    0],
  [  2.2,    0,    0],
  [  2.3,    0,    0],
  [  2.4,    0,    0],
  [  2.5,    0,    0],
  [  2.6,    0,    0],
  [  2.7,    0,    0],
  [  2.8,    0,    0],
  [  2.9,    0,    0],
  [    3,    0,    0],
  [  3.1,    0,    0],
  [  3.2,    0,    0],
  [  3.3,    0,    0],
  [  3.4,    0,    0],
  [  3.5,    0,    0],
  [  3.6,    0,    0],
  [  3.7,    0,    0],
  [  3.8,    0,    0],
  [  3.9,    0,    0],
  [    4,    0,    0],
  [  4.1,    0,    0],
  [  4.2,    0,    0],
  [  4.3,    0,    0],
  [  4.4,    0,    0],
  [  4.5,    0,    0],
  [  4.6,    0,    0],
  [  4.7,    0,    0],
  [  4.8,    0,    0],
  [  4.9,    0,    0],
  [    5,    0,    0]]
 """

# trying more crazy things also not working to set selection
selObjects = [rr.createSelection(item) for item in sel]
rr.selections = selObjects
result = rr.simulate(start=0, end=10, steps=5)
print result
