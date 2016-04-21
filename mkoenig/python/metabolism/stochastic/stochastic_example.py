"""

Run stochastic model simulation with roadrunner
"""

from __future__ import print_function
import os
import roadrunner
import pandas

# load test file
test_sbml = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         './models/00001/00001-sbml-l3v1.xml')
rr = roadrunner.RoadRunner(test_sbml)

# deterministic
result = rr.simulate(start=0, end=50, steps=50, integrator="cvode")
rr.plot()

# stochastic
rr.reset()
result = rr.simulate(start=0, end=50, steps=50, integrator="gillespie")
rr.plot()

# repeated simulation
all_results = []
for k in xrange(0, 100):
    print(k)
    rr.reset()
    result = rr.simulate(start=0, end=50, steps=50, integrator="gillespie")
    # rr.plot()
    all_results.append(result['[X]'])

df = pandas.DataFrame(data=all_results)   
df = df.transpose()
df.plot()
