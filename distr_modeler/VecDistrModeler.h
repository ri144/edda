// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef VEC_DISTRIBUTION_MODELER_H
#define VEC_DISTRIBUTION_MODELER_H

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

#include "core/vector_matrix.h"
#include "core/shared_ary.h"


using namespace edda;
using namespace std;
using namespace edda::dist;


using namespace std;

namespace edda{

// Distribution modeler
template <typename T>
class VecDistrModeler {
protected:
	shared_ptr<Dataset<VECTOR3> > dataset;
	size_t xSize, ySize, zSize;
	T* distArray;
	size_t binCount;
	size_t numComp; //number of components.
public:
    VecDistrModeler()
   	{
        xSize = 1;
        ySize = 1;
        zSize = 1;
		binCount = 32;
		numComp = 3;
	}
	void initVectorComponent(size_t n)
	{
		numComp = n;
	}
	void assignGrid(int X, int Y, int Z)
	{
		xSize = X;
		ySize = Y;
		zSize = Z;
		distArray = new T[xSize*ySize*zSize*numComp];
	}
	void initHistogram(size_t bin)
	{
		binCount = bin;
	}
	void computeGMM(float *data, size_t size, size_t index)
	{
		for(size_t i=0; i<numComp; i++)
	 	{
	 		double *dataD;
	 		dataD = new double[size];
	 		for(size_t j=0; j<size; j++)
	 		{
	 			dataD[j] = (double) data[j*numComp + i];
	 		}
	 		T new_distr;
	 		eddaComputeEM(dataD, size, &new_distr);
	 		distArray[index*numComp + i] = new_distr;
	 	}
		

		
	}
	void computeHistogram(float *data, size_t size, size_t index)
	{

		for(size_t i=0; i<numComp; i++)
		{
			float *dataF;
	 		dataF = new float[size];
	 		for(size_t j=0; j<size; j++)
	 		{
	 			dataF[j] = data[j*numComp + i];
	 		}
	 		T new_distr;
	 		new_distr = eddaComputeHistogram(data, size, binCount);
	 		distArray[index*numComp + i] = new_distr;
		}


	}
	void model()
	{
		shared_ary<Vector<T,3> > pArray(new Vector<T, 3>[xSize*ySize*zSize], xSize*ySize*zSize);
		for(int i=0; i<xSize*ySize*zSize*3; i += 3)
		{
			Vector<T, 3> curVectorDistr;
			for(int j=0; j<3; j++)
			{
				curVectorDistr[j] = distArray[i+j];
			}
			pArray[i/3] = curVectorDistr;
		}
		DistrArray * abs_array = new VectorDistrArray<T,3>(pArray);
	 	dataset = make_Dataset<VECTOR3>(new RegularCartesianGrid(xSize, ySize, zSize), abs_array);
	}
	/*void writeToVTK(const string &edda_fname, const string &array_name_prefix){
		writeEddaVtkDataset(dataset, edda_fname, array_name_prefix);
	}*/
	shared_ptr<Dataset<VECTOR3>> getDataset(){
		return dataset;
	}


};

}
#endif // VEC_DISTRIBUTION_MODELER_H
