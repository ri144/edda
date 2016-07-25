// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DISTR_MODELER_H
#define DISTR_MODELER_H

#include <iostream>
#include <string>
#include <vector>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include <distributions/gaussian_mixture.h>
#include "distributions/histogram.h"
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"
#include "dataset/dataset.h"
#include "dataset/abstract_distr_array.h"
#include "dataset/distr_array.h"
#include "core/interpolator.h"
#include "io/edda_vtk_writer.h"

#include <netcdfcpp.h>

using namespace std;


namespace edda{

// Distribution modeler
class DistrModeler {
protected:
	DistrType dType;
	const int NC_ERR = 2;
	shared_ptr<Dataset<Real> > dataset;
	size_t binCount;
public:
        DistrModeler(DistrType type){
		//currently updating the type of distribution (GMM/HIST)
		//can have more variables in future to hold different settings while modelling,  like partitioning algortihm etc.
		dType = type;		
	}

	//loader() will be overloaded to handle different types of input data and partitioning styles
	void loader(string filename, string xDimName, string yDimName, string zDimName, string ensDimName, string varName)
	{
		switch(dType){
			case(GMM):
				loader<dist::DefaultGaussianMixture>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
			case(GMM2):
				loader<dist::GaussianMixture<2>>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
			case(GMM3):
				loader<dist::GaussianMixture<3>>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
			case(GMM4):
				loader<dist::GaussianMixture<4>>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
			case(GMM5):
				loader<dist::GaussianMixture<5>>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
			case(HIST):
				loader<dist::Histogram>(filename, xDimName, yDimName, zDimName, ensDimName, varName);
				break;
		}
	}

	template <typename T>
	void loader(string filename, string xDimName, string yDimName, string zDimName, string ensDimName, string varName)
	{
		//input data parameters
		size_t xDim(0), yDim(0), zDim(0), ensDim(0);

		//try loading the netcdf file.
		NcFile dataFile(filename.c_str(), NcFile::ReadOnly);
		if(!dataFile.is_valid())
  		{
    		cerr << "ERROR[" << NC_ERR << "]" << "reading file" << endl;
    		exit(NC_ERR);
  		}

  		//update the data dimensions by checking the nc file.
  		xDim = dataFile.get_dim(xDimName.c_str())->size();
  		yDim = dataFile.get_dim(yDimName.c_str())->size();
  		zDim = dataFile.get_dim(zDimName.c_str())->size();
  		ensDim = dataFile.get_dim(ensDimName.c_str())->size();
                ensDim = 10;

  		NcVar *var = dataFile.get_var(varName.c_str());

  		cout << "Dimensions : (" << xDim << "," << yDim << "," << zDim << ")" << endl;
  		cout << "Ensemble Count : " << ensDim << endl;
  	
  		float *data;
  		

  		data = new float[ensDim];

  		//based on the dType provided by the user create appropiate array type.  		
  		shared_ary<T> pArray (new T[xDim*yDim*zDim], xDim*yDim*zDim);
  		

	 	for(size_t z=0; z<zDim; z++)
	 	{
	 		for(size_t y=0; y<yDim; y++)
	 		{
	 			for(size_t x=0; x<xDim; x++)
	 			{
	 				if(!var->set_cur(0, x, y, z))
	 				{
	 					cerr << "ERROR[" << NC_ERR << "]" << "setting current pointer" << endl;
	 					exit(NC_ERR);	
	 				}
	 				if(!var->get(data, ensDim, 1, 1, 1))
	 				{
	 					cerr << "ERROR[" << NC_ERR << "]" << "getting data" << endl;
	 					exit(NC_ERR);
	 				}
	        		//TODO::Will have to create an API to handle Histogram
	 				
	 				

	 				T new_distr;
	 				computeDistribution(data, ensDim, new_distr);
	 				
	        		pArray[z*xDim*yDim + y*xDim + x] = new_distr;


	 			}
	 		}
	 	}
	 	AbstractDistrArray * abs_array = new DistrArray<T>(pArray);
	 	dataset = make_Dataset<Real>(new RegularCartesianGrid(xDim,yDim,zDim), abs_array);

	}

	template <int GMs> //GMs: Gaussian models
	void computeDistribution(float *data, size_t ensDim, dist::GaussianMixture<GMs> &new_distr)
	{
		double *dataD;		
  		dataD = new double[ensDim];

  		for(size_t i=0; i<ensDim; i++)
	 		dataD[i] = (double) data[i];

		eddaComputeEM(dataD,ensDim, &new_distr);
	}

	void computeDistribution(float *data, size_t ensDim, dist::Histogram &new_distr)
	{
		new_distr = eddaComputeHistogram(data, ensDim, -10, 10, 20);
	}

	DistrType getType(){
		return dType;
	}
	
	void writeToVTK(const string &edda_fname, const string &array_name_prefix){
		writeEddaVtkDataset(dataset, edda_fname, array_name_prefix);
	}

	shared_ptr<Dataset<Real>> getDataset(){
		return dataset;
	}


};

}
#endif // DISTR_MODELER_H
