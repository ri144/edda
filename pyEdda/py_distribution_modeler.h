#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "distributions/distribution_modeler.h"


using namespace std;
using namespace edda;


namespace py = pybind11;


class PyDistributionModeler {
protected:
	DistributionModeler* dm;
	size_t len;
public:
	/// \brief Constructor with specified length
	PyDistributionModeler(int l){
		len = l;
        dm = new DistributionModeler(l);
	}
	void computeGMM(py::array_t<Real> data, size_t nGmm, size_t index){
		//first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = data.request();
		auto ptr = static_cast<Real *>(info.ptr);

		int ndim = info.ndim;
		if(ndim != 1)
			throw std::runtime_error("Number of dimensions must be one");

		int size = info.shape[0];

		
		dm->computeGMM(ptr, size, nGmm, index);
		
	}
	void computeHistogram(py::array_t<Real> data, int binCount, size_t index ){

		//first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = data.request();
		auto ptr = static_cast<Real *>(info.ptr);

		int ndim = info.ndim;
		if(ndim != 1)
			throw std::runtime_error("Number of dimensions must be one");

		int size = info.shape[0];

		

		dm->computeHistogram(ptr, size, index, binCount);

	}
	void printDistr(){
		DistrArray * darray = dm->getDistrArray(); 
		for(int i=0; i<len; i++)
		{
			cout << "DistrArray[" << i << "]" << darray->getDistr(i) << endl << endl;
		}
	}
	auto getDistr(int index){
		DistrArray * darray = dm->getDistrArray(); 
		dist::Variant curDist = darray->getDistr(index);
		string s = getName(curDist);

		if (s.compare(0, 15, "GaussianMixture") == 0) {
			// Only compare the first 15 chars because this string ends with the number of Gaussian models
			// distr type. 1: GaussianMixture. 2: Histogram
			// if GaussianMixture:
			int distrTypeNumber = 1;
			
			// number of gaussian components
			int GMs = stoi(s.substr(15));

			py::array_t<float> rv1(GMs);
			auto r1 = rv1.mutable_unchecked<1>();

			py::array_t<float> rv2(GMs);
			auto r2 = rv2.mutable_unchecked<1>();

			py::array_t<float> rv3(GMs);
			auto r3 = rv3.mutable_unchecked<1>();
			
			//	GMM parameters
			std::vector<dist::GMMTuple> gmmtv = boost::get<dist::GMM >(curDist).models;
			int c1=0,c2=0,c3=0;
			for (int i = 0; i < GMs * 3; i++){
				if (i % 3 == 0)
				{
					r1[c1] = gmmtv[i / 3].m;
					c1++;
				}
				else if (i % 3 == 1)
				{
					r2[c2] = gmmtv[i / 3].v;
					c2++;
				}
				else
				{
					r3[c3] = gmmtv[i / 3].w;
					c3++;
				}
			}
			
			//return multiple values in python call by returning a tuple in c++ side
			return std::make_tuple(rv1,rv2,rv3,s);
			
		}
		else if (s.compare(0, 15, "Histogram") == 0) {
			//	if Histogram:
			int distrTypeNumber = 2;
			
			dist::Histogram curHist = boost::get<dist::Histogram>(curDist);

			int nbins = curHist.getBins();
			float minv = curHist.getMinValue();
			float maxv = curHist.getMaxValue();

			py::array_t<float> rv1(1);
			auto r1 = rv1.mutable_unchecked<1>();

			py::array_t<float> rv2(1);
			auto r2 = rv2.mutable_unchecked<1>();

			py::array_t<float> rv3(nbins);
			auto r3 = rv3.mutable_unchecked<1>();

			r1[0] = minv;
			r2[0] = maxv;
			
			//	bin values
			for (int b = 0; b < nbins; b++)
			{
				r3[b] = curHist.getBinValue(b);
			}

			//return multiple values in python call by returning a tuple in c++ side
			return std::make_tuple(rv1,rv2,rv3,s);
			
		}



	}
	
	
};


