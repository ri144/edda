#make print in python 2, 3 compatible
from __future__ import print_function 
import numpy as np
import pyedda as edda
from scipy import misc


#JointHistogram
print("Testing code for JointHistogram with an image example...")
print("This testing code requires the scipy package!")

im = misc.imread("../sample_data/test_img.bmp")
width, height, channel = (im.shape)
print ((im.dtype))

bz = 20	#block size
dsUs = width/bz
dsVs = height/bz

jointMeanImg = np.ndarray(shape=(width, height, channel), dtype=np.int)
jointSmpImg = np.ndarray(shape=(width, height, channel), dtype=np.int)
marginMeanImg = np.ndarray(shape=(width, height, channel), dtype=np.int)
marginSmpImg = np.ndarray(shape=(width, height, channel), dtype=np.int)

for dsU in range(dsUs):
	for dsV in range(dsVs):
		var0 = []
		var1 = []
		var2 = []

		var0 = im[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 0]
		var1 = im[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 1]
		var2 = im[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 2]

		npdata = []
		npdata.append(var0.flatten())
		npdata.append(var1.flatten())
		npdata.append(var2.flatten())
		npdata = np.array(npdata)

		# print(npdata.shape) # should be 3, blckSize*blockSize

		bins = np.array([30, 30, 30]) # each var has 30 bins
		mins = np.amin(npdata, axis=1) # min of each var
		maxs = np.amax(npdata, axis=1) # max of each var
		
		# construct joint histogram for the three variables (for the current block)
		jointHist = edda.JointHistogram(npdata, npdata.shape[1], mins, maxs, bins)

		##########################################################################
		#                  finish construction, start testing                    #
		##########################################################################
		# 1. test joint mean
		mean = jointHist.getJointMean() # return joint mean as a numpy array
		#mean = edda.getJointMean(jointHist); # another way to call it
		jointMeanImg[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz] = mean

		# 2. test joint sample
		for u in range(dsU*bz, (dsU+1)*bz):
			for v in range(dsV*bz, (dsV+1)*bz):
				jointSmpImg[u, v] = jointHist.getJointSample()

		# 3. test marginalization
		marvars = np.array([1, 2]) # marginalize to dimension 1 and 2
		marginHist = edda.marginalization(jointHist, marvars)
		marginMean = marginHist.getJointMean()
		# print ("The dimension fo the marginalized mean is {}".format(marginMean.shape))
		marginMeanImg[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 0] = 0
		marginMeanImg[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 1] = marginMean[0]
		marginMeanImg[dsU*bz:(dsU+1)*bz, dsV*bz:(dsV+1)*bz, 2] = marginMean[1]

		# for u in range(dsU*bz, (dsU+1)*bz):
		# 	for v in range(dsV*bz, (dsV+1)*bz):
		# 		marginSmp = marginHist.getJointSample()
		# 		marginSmpImg[u, v, 0] = 0
		# 		marginSmpImg[u, v, 1] = marginMean[0]
		# 		marginSmpImg[u, v, 2] = marginMean[1]

		# 4. test conditionalHistogram
		# orivars = 
		# condvars = 

misc.imsave('./joint_mean.bmp', jointMeanImg)
misc.imsave('./joint_sample.bmp', jointSmpImg)
# result for marginalized distribution
misc.imsave('./margin_mean.bmp', marginMeanImg)
misc.imsave('./margin_sample.bmp', marginSmpImg)