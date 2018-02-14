#make print in python 2, 3 compatible
from __future__ import print_function 
import numpy as np
import pyedda as edda


#Univariate Gaussian
print("//////////Univariate Gaussian///////")
dummy_data = np.random.rand(100)
gaussian = edda.Gaussian(100, 20)
print("gaussian.getMean():", gaussian.getMean())
print("gaussian.getVar():", gaussian.getVar())
print("gaussian.getPdf(105):", gaussian.getPdf(105))
print("gaussian.getSample():", gaussian.getSample())
print("gaussian.getCdf(105):", gaussian.getCdf(105))
print("gaussian.getCdfPrecise():", gaussian.getCdfPrecise(105))
print("Output gaussian:")
gaussian.output()
print()
