#make print in python 2, 3 compatible
from __future__ import print_function 
import numpy as np
import pyedda as edda

#Univariate histogram
print("//////////Univariate Histogram///////")
dummy_data = np.random.rand(100)
hist = edda.Histogram(dummy_data, 10)
print ("hist.getMean():", hist.getMean())
histCopy = edda.Histogram(hist)
print("histCopy.getMean():", histCopy.getMean())
print("hist.getVar():", hist.getVar())
print("hist.getPdf(0.5):",hist.getPdf(0.5))
print("hist.getCdf(1.0):",hist.getCdf(1.0))
print("hist.getSample():", hist.getSample())
print("Output histogram:")
hist.output()
print("hist.getBins():", hist.getBins())
print("hist.getMaxValue():", hist.getMaxValue())
print("hist.getMinValue():", hist.getMinValue())
print("hist.getBinValue(3):", hist.getBinValue(3))
print()