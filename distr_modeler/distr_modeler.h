#ifndef DISTR_MODELER_H
#define DISTR_MODELER_H

#include <iostream>
#include <string>
#include <vector>
#include "edda.h"
#include "distributions/gaussian.h"
#include "distributions/distribution.h"
#include <distributions/gaussian_mixture.h>
#include <distributions/estimate_gmm.h>
#include "distributions/variant.h"
#include "dataset/dataset.h"
#include "dataset/abstract_distr_array.h"
#include "dataset/distr_array.h"
#include "core/interpolator.h"

#include <netcdfcpp.h>

using namespace std;


namespace edda{

template <typename T>
class DistrModeler {
protected:
	DistrType dType;
	const int NC_ERR = 2;
	shared_ptr<Dataset<T> > dataset;
public:
	DistrModeler(DistrType type){
		//currently updating the type of distribution (GMM/HIST)
		//can have more variables in future to hold different settings while modelling like partitioning algortihm etc.
		dType = type;
		
	}
	//loader() will be overloaded to handle different types of input data and partitioning styles
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
  		double *dataD;

  		data = new float[ensDim];
  		dataD = new double[ensDim];
  		
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
	 				T new_gmm;
	 				for(size_t i=0; i<ensDim; i++)
	 					dataD[i] = (double) data[i];


    				eddaComputeEM(dataD,ensDim, &new_gmm);
    				//cout << "[" <<z*xDim*yDim + y*xDim + x << "]" <<new_gmm << endl;
	 				
	        		pArray[z*xDim*yDim + y*xDim + x] = new_gmm;


	 			}
	 		}
	 	}

	 	//cout << "some test results:\n" ;
	 	//cout << "length = " << pArray.getLength() << endl;
	 	//cout << "value = " << pArray[11255]<< endl;

	 	AbstractDistrArray * abs_array = new DistrArray<T>(pArray);
	 	dataset = make_Dataset<T>(new RegularCartesianGrid(xDim,yDim,zDim), abs_array);

	}

	DistrType getType(){
		return dType;
	}
	/*T getDistributionAt(int i, int j, int k){
		return dataset->at_comp_distr(i,j,k);
	}*/

	void writer(){
		//TODO: write the dataset to .vti files
	}
};

}
#endif // DISTR_MODELER_H
