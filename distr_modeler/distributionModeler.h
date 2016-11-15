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
#include "distributions/histogram.h"
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"
#include "dataset/dataset.h"
#include "dataset/distr_array.h"
#include "core/interpolator.h"
#include "io/edda_vtk_writer.h"

using namespace std;

namespace edda{

// Distribution modeler
template <typename T>
class DistributionModeler {
protected:
	shared_ptr<Dataset<Real> > dataset;
	size_t xSize, ySize, zSize;
	T* distArray;
	size_t binCount;
public:
    DistributionModeler()
   	{
        xSize = 1;
        ySize = 1;
        zSize = 1;
		binCount = 32;
	}
	void assignGrid(int X, int Y, int Z)
	{
		xSize = X;
		ySize = Y;
		zSize = Z;
		distArray = new T[xSize*ySize*zSize];
	}
	void initHistogram(size_t bin)
	{
		binCount = bin;
	}
	void computeGMM(float *data, size_t size, size_t index)
	{
		double *dataD;		
  		dataD = new double[size];

  		for(size_t i=0; i<size; i++)
	 		dataD[i] = (double) data[i];

	 	T new_distr;
		eddaComputeEM(dataD, size, &new_distr);

		distArray[index] = new_distr;
	}
	void computeHistogram(float *data, size_t size, size_t index)
	{
		T new_distr;
		new_distr = eddaComputeHistogram(data, size, binCount);

		distArray[index] = new_distr;
	}
	void model()
	{
		shared_ary<T> pArray (new T[xSize*ySize*zSize], xSize*ySize*zSize);
		for(int i=0; i<xSize*ySize*zSize; i++)
		{
			pArray[i] = distArray[i];
		}
		DistrArray * abs_array = new ScalarDistrArray<T>(pArray);
	 	dataset = make_Dataset<Real>(new RegularCartesianGrid(xSize, ySize, zSize), abs_array);
	}
	void writeToVTK(const string &edda_fname, const string &array_name_prefix){
		writeEddaVtkDataset(dataset, edda_fname, array_name_prefix);
	}
	shared_ptr<Dataset<Real>> getDataset(){
		return dataset;
	}


};

}
#endif // DISTRIBUTION_MODELER_H
