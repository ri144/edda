#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "distributions/gmm.h"

using namespace std;
using namespace edda::dist;

namespace py = pybind11;

class PyGMM {
protected:
	GMM* gmm;

public:
    PyGMM()
    {
        gmm = new GMM();
    }

    PyGMM(int gms)
    {
        gmm = new GMM(gms);
    }

    PyGMM(py::array_t<Real> p_model)
    {
        //first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = p_model.request();

        //Check the input numpy array's dimensionality
        int ndim = info.ndim;
		if(ndim != 2)
			throw std::runtime_error("Number of dimensions must be 2");
        
        if( info.shape[1] != 3 )
            throw std::runtime_error("Length of the 2nd dimensions must be 3");

        std::vector<GMMTuple> c_Models;
        auto parameter = p_model.unchecked<2>();
        for( int i=0; i<info.shape[0]; i++ ){
            GMMTuple tuple;
            tuple.w = parameter(i,0);
            tuple.m = parameter(i,1);
            tuple.v = parameter(i,2);
            c_Models.push_back(tuple);
        }
        
        gmm = new GMM(c_Models);
    }

    PyGMM(py::array_t<double> data, int GMs)
    {
        //first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = data.request();
		auto ptr = static_cast<double *>(info.ptr);

        //Check the input numpy array's dimensionality
        int ndim = info.ndim;
		if(ndim != 1)
			throw std::runtime_error("Number of dimensions must be one");

        int nElements = info.shape[0];

        GMM tmpGmm = eddaComputeGMM(ptr, nElements, GMs);
        gmm = new GMM(tmpGmm.models);
        
    }

    void assign(py::array_t<Real> p_model)
    {
        //first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = p_model.request();

        //Check the input numpy array's dimensionality
        int ndim = info.ndim;
		if(ndim != 2)
			throw std::runtime_error("Number of dimensions must be 2");
        
        if( info.shape[1] != 3 )
            throw std::runtime_error("Length of the 2nd dimensions must be 3");

        std::vector<GMMTuple> c_Models;
        auto parameter = p_model.unchecked<2>();
        for( int i=0; i<info.shape[0]; i++ ){
            GMMTuple tuple;
            tuple.w = parameter(i,0);
            tuple.m = parameter(i,1);
            tuple.v = parameter(i,2);
            c_Models.push_back(tuple);
        }
        
        gmm->assign(c_Models);
    }

    int getNumComponents()
    {
        return gmm->getNumComponents();
    }

    void normalizeWeights()
    {
        gmm->normalizeWeights();
    }

    double getGMMMean()
    {
        return getMean(*gmm);
    }

    double getGMMVar()
    {
        return getVar(*gmm);
    }

    double getGMMPdf(const double x)
    {
        return getPdf(*gmm, x);
    }

    double getGMMSample()
    {
        return getSample(*gmm);
    }

    double getGMMCdf(const double x)
    {
        return getCdf(*gmm, x);
    }

    void output()
    {
        std::cout << *gmm << std::endl;
    }
};