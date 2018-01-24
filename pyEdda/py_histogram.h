#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "distributions/histogram.h"

using namespace std;
using namespace edda::dist;

namespace py = pybind11;

class PyHistogram {
protected:
	Histogram* hist;

public:
    PyHistogram(py::array_t<Real> data, const int _nBins)
    {
        //first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = data.request();
		auto ptr = static_cast<Real *>(info.ptr);

        //Check the input numpy array's dimensionality
        int ndim = info.ndim;
		if(ndim != 1)
			throw std::runtime_error("Number of dimensions must be one");

        int nElements = info.shape[0];

        hist = new Histogram(ptr, nElements, _nBins);
    }

    PyHistogram(py::array_t<Real> data, const int _nBins, 
                const Real _minValue, const Real _maxValue)
    {
        //first convert the numpy array (data) to a pointer of array to be passed to the c++ functions
		py::buffer_info info = data.request();
		auto ptr = static_cast<Real *>(info.ptr);

        //Check the input numpy array's dimensionality
        int ndim = info.ndim;
		if(ndim != 1)
			throw std::runtime_error("Number of dimensions must be one");

        int nElements = info.shape[0];

        hist = new Histogram(ptr, nElements, _nBins, _minValue, _maxValue);
    }

    PyHistogram(PyHistogram* pyHist)
    {
        hist = pyHist->hist;
    }

    Real getMean() const
    {
        return hist->getMean();
    }

    Real getVar() const
    {
        return hist->getVar();
    }

    Real getPdf(const double x) const
    {
        return hist->getPdf(x);
    }

    Real getCdf(const double x) const
    {
        return hist->getCdf(x);
    }

    Real getSample() const
    {
        return hist->getSample();
    }

    void output() const
    {
        std::cout << *hist << std::endl;
    }

    int getBins()
    {
        return hist->getBins();
    }

    float getMaxValue(){
	    return hist->getMaxValue();
    }

    float getMinValue(){
        return hist->getMinValue();
    }

    float getBinValue(int b){
	    return hist->getBinValue(b);
    }
};