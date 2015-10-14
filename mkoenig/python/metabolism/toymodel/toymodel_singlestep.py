"""
Try to run a single step simulation.
@author: Matthias Koenig
"""
import roadrunner
print roadrunner.__version__
# 1.4.0; Compiler: gcc 4.8.4, C++ version: 199711; JIT Compiler: LLVM-3.4; Date: Oct 14 2015, 00:00:59

# Load model results in warnings about fbc v2
rr = roadrunner.RoadRunner("toy_ode_bounds.xml")

# simulate multiple steps
result = rr.simulate(start=0, end=10, steps=5)
print result

# simulate multiple steps with variable steps
result = rr.simulate(start=0, end=1, variableStep=True)
print result

# simulate single setp with variable steps
# ! Not simulating a single step, but reusing the settings for end from previous simulation
# How to get a single setp variableStep simulation ????
result = rr.simulate(start=0, steps=1, variableStep=True)
print result
