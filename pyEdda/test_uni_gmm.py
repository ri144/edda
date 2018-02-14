#make print in python 2, 3 compatible
from __future__ import print_function 
import numpy as np
import pyedda as edda

#Univariate GMM
print("//////////Univariate GMM///////")
gmm = edda.GMM(4)
models = np.array([[0.2, 10, 2],[0.5, 15, 4],[0.1, 8, 7],[0.2,6, 9]])
gmm = edda.GMM(models)
print("gmm.getNumComponents():", gmm.getNumComponents() )
print("gmm.getMean():", gmm.getMean())
print("gmm.getVar():", gmm.getVar())
print("gmm.getPdf(10):", gmm.getPdf(10))
print("gmm.getSample():", gmm.getSample())
print("gaussian.getCdf(10):", gmm.getCdf(10))
print("Output GMM:")
gmm.output()
print("Model and output a gmm by trainig samples:")
dummy_data = np.random.rand(100)
modelGmm = edda.GMM(dummy_data,3)
modelGmm.output()
print()
