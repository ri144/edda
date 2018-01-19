import numpy as np
import pyedda as edda


dm = edda.DistributionModeler(2)

dummy_data = np.random.rand(100)

dm.computeGMM(dummy_data,2,0)
dm.computeHistogram(dummy_data,32,1)

dm.printDistr()



r1,r2,r3,r4 = dm.getDistr(1)

print r1
print r2
print r3
print r4
