// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTRIBUTION_MODELER_H
#define DISTRIBUTION_MODELER_H

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include <distributions/gaussian_mixture.h>
#include <distributions/gmm.h>
#include "distributions/histogram.h"
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"

#include "distributions/joint_GMM.h"
#include "distributions/joint_gaussian.h"
#include "distributions/joint_histogram.h"

#include "dataset/dataset.h"
#include "dataset/distr_array.h"
#include "core/interpolator.h"
#include "io/edda_vtk_writer.h"

using namespace std;

namespace edda{

// Distribution modeler

class DistributionModeler {
protected:
	size_t len;
	dist::Variant* distArray;
public:
    DistributionModeler()
   	{
        len = 0;
	}
	DistributionModeler(int l)
   	{
        len = l;
        distArray = new dist::Variant[len];
	}
	void computeGMM(float *data, size_t size, size_t nGmm, size_t index)
	{
		//TODO: check for out of bound errors.
		double *dataD;		
  		dataD = new double[size];

  		for(size_t i=0; i<size; i++)
	 		dataD[i] = (double) data[i];

	 	dist::GMM new_distr;
		new_distr = eddaComputeGMM(dataD, size, nGmm);

		distArray[index] = new_distr;
	}
	void computeHistogram(float *data, size_t size, size_t index, int binCount)
	{
		dist::Histogram new_distr;
		new_distr = eddaComputeHistogram(data, size, binCount);

		distArray[index] = new_distr;
	}
	void computeJointGMM(std::vector<Real*>& data, size_t size, size_t nGmm, size_t index)
	{
		dist::JointGMM new_distr;
		new_distr = eddaComputeJointGMM(data, size, nGmm);
		
		distArray[index] = new_distr;
	}

	//temporary test function. will be deleted afterwards
	void setJointGMM(dist::JointGMM  new_d, size_t index)
	{
		//distArray[index] = new_distr;
		dist::JointGMM new_distr;
		new_distr = new_d;
		distArray[index] = new_distr;

		cout << "in modeler func, joint GMM No. " << index << ":" << endl;
		cout << distArray[index] << endl;
	}

	//temporary test function. will be deleted afterwards
	void setJointGMM(int nVar, int nComp, ublas_vector &w, ublas_matrix &m, ublas_matrix &covs, size_t index)
	{
		dist::JointGMM new_distr;
		new_distr.setGMM(nVar, nComp, w, m, covs);
		distArray[index] = new_distr;
	}

	void computeJointHistogram(std::vector<Real*>& data, size_t size, const std::vector<Real>& mins, const std::vector<Real>& maxs, const std::vector<int>& nBins, size_t index)
	{
		dist::JointHistogram new_distr;
		new_distr = dist::eddaComputeJointHistogram(data, size, mins, maxs, nBins);
		
		distArray[index] = new_distr;
	}
	DistrArray * getDistrArray()
	{
		shared_ary<dist::Variant> pArray (new dist::Variant[len], len);
		for(int i=0; i<len; i++)
		{
			pArray[i] = distArray[i];
		}
		DistrArray * abs_array = new ScalarDistrArray<dist::Variant>(pArray);
		return abs_array;
	}
	
};

}
#endif // DISTRIBUTION_MODELER_H
