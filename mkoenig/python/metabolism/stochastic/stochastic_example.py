"""
Run a stochastic model simulation.

"""

import roadrunner
# rr = roadrunner.RoadRunner('./data/00001/00001-sbml-l3v1.xml')
# rr = roadrunner.RoadRunner("00001-sbml-l3v1.xml")
rr = roadrunner.RoadRunner("test.xml")

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
    print k
    rr.reset()
    result = rr.simulate(start=0, end=50, steps=50, integrator="gillespie")
    # rr.plot()
    all_results.append( result['[X]'] )

import pandas
df = pandas.DataFrame(data=all_results)   
df = df.transpose()

df.plot()
