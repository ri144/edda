#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "distributions/joint_histogram.h"

using namespace std;
using namespace edda::dist;

namespace py = pybind11;

class PyJointHistogram {
public:
	JointHistogram* hist;

public:
    /// constructors
    PyJointHistogram(JointHistogram* pyHist) {
        hist = pyHist;
    }

    PyJointHistogram(PyJointHistogram* pyHist) {
        hist = pyHist->hist;
    }

    PyJointHistogram(py::array_t<Real> data, int nElements, py::array_t<Real> mins, py::array_t<Real> maxs, py::array_t<int> nBins) {
        // convert data from 2D numpy array to vector of Real pointers, i.e. vector<Real*>
        py::buffer_info info = data.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 2)
            throw std::runtime_error("Number of dimensions for the data array must be two, one for variables, one for elemetns of the vars");
        int nVars = info.shape[0];
        nElements = info.shape[1]; // these two values should be the same, if not, use numpy array's second dimension

        std::vector<Real*> c_data;
        for(int i = 0; i < nVars; i++) {
            auto ptr = static_cast<Real *>((Real*)info.ptr + i * nElements);
            c_data.push_back(ptr);
        }
        // convert mins array from numpy array of floats to vector of floats
        info = mins.request();
        ptr = static_cast<Real *>(info.ptr);
        ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the mins array must be one");
        int nMins = info.shape[0];
        std::vector<Real> c_mins {ptr, ptr + nMins};


        // convert maxs array from numpy array of floats to vector of floats
        info = maxs.request();
        ptr = static_cast<Real *>(info.ptr);
        ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the maxs array must be one");
        int nMaxs = info.shape[0];
        std::vector<Real> c_maxs {ptr, ptr + nMaxs};

        // convert nBins array from numpy array of ints to vector of ints
        info = maxs.request();
        auto ptrint = static_cast<int *>(info.ptr);
        ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the nBins array must be one");
        int nnBins = info.shape[0];
        std::vector<int> c_nBins {ptrint, ptrint + nnBins};

        hist = new JointHistogram(c_data, nElements, c_mins, c_maxs, c_nBins);
    }

    /// seter and geter
    void setMinVals(py::array_t<Real> mins) {
        py::buffer_info info = mins.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the mins array must be one");
        int nMins = info.shape[0];
        std::vector<Real> c_mins {ptr, ptr + nMins};
        hist->setMinVals(c_mins);
    }

    py::array_t<Real> getMinVals() {
        std::vector<Real> v = hist->getMinVals();
        return py::array_t<Real>(v.size(), v.data());
    }

    void setMaxVals(py::array_t<Real> maxs) {
        py::buffer_info info = maxs.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the mins array must be one");
        int nMaxs = info.shape[0];
        std::vector<Real> c_maxs {ptr, ptr + nMaxs};
        hist->setMinVals(c_maxs);
    }

    py::array_t<Real> getMaxVals() {
        std::vector<Real> v = hist->getMaxVals();
        return py::array_t<Real>(v.size(), v.data());
    }

    void setBinWidths(py::array_t<Real> widths) {
        py::buffer_info info = widths.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the mins array must be one");
        int nWidths = info.shape[0];
        std::vector<Real> c_widths {ptr, ptr + nWidths};
        hist->setBinWidths(c_widths);
    } 

    py::array_t<Real> getBinWidths() {
        std::vector<Real> v = hist->getBinWidths();
        return py::array_t<Real>(v.size(), v.data());
    }

    void setNumBins(py::array_t<int> nBins) {
        py::buffer_info info = nBins.request();
        auto ptr = static_cast<int *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions for the mins array must be one");
        int nnBins = info.shape[0];
        std::vector<int> c_nBins {ptr, ptr + nnBins};
        hist->setNumBins(c_nBins);
    } 

    py::array_t<int> getNumBins() {
        std::vector<int> v = hist->getNumBins();
        return py::array_t<int>(v.size(), v.data());
    }

    void setNumVars(int nVars) {
        hist->setNumVars(nVars);
    } 

    int getNumVars() {
        return hist->getNumVars();
    }

    py::array_t<Real> getJointMean() {
        std::vector<Real> v = hist->getJointMean();
        return py::array_t<Real>(v.size(), v.data());
    }

    py::array_t<Real> getJointSample() {
        std::vector<Real> v = hist->getJointSample();
        return py::array_t<Real>(v.size(), v.data());
    }
};


py::array_t<Real> getJointMean_py(PyJointHistogram& hist) {
    return hist.getJointMean();
}

PyJointHistogram marginalization(PyJointHistogram& pyHist, py::array_t<int> vars) {
    std::unordered_set<int> marvars;

    py::buffer_info info = vars.request();
    auto ptr = static_cast<int *>(info.ptr);
    int ndim = info.ndim;
    if(ndim != 1)
        throw std::runtime_error("Variables to be marginalized must be in 1D");
    int num_vars = info.shape[0];
    for(int i=0; i<num_vars; i++) {
        marvars.insert(*(ptr+i));
    }

    PyJointHistogram mhist(pyHist.hist->marginalization(marvars));
    return mhist;
}

PyJointHistogram conditionalHist(PyJointHistogram& pyHist, py::array_t<int> orivars, py::array_t<int> condvars, py::array_t<int> condbins) {
    std::unordered_set<int> c_orivars;
    py::buffer_info info = orivars.request();
    auto ptr = static_cast<int *>(info.ptr);
    int ndim = info.ndim;
    if(ndim != 1)
        throw std::runtime_error("Dimensions of Variables must be in 1D");
    int num_vars = info.shape[0];
    for(int i=0; i<num_vars; i++) {
        c_orivars.insert(*(ptr+i));
    }

    std::vector<int> c_condvars;
    info = condvars.request();
    ptr = static_cast<int *>(info.ptr);
    ndim = info.ndim;
    if(ndim != 1)
        throw std::runtime_error("Dimensions of conditional variables must be in 1D");
    num_vars = info.shape[0];
    for(int i=0; i<num_vars; i++) {
        c_condvars.push_back(*(ptr+i));
    }

    std::vector<std::pair<int,int>> c_binrange;
    info = condbins.request();
    ptr = static_cast<int *>(info.ptr);
    ndim = info.ndim;
    if(ndim != 2)
        throw std::runtime_error("Dimension of conditional range must be in 2D");
    num_vars = info.shape[0];
    int bin_range = info.shape[1]; // must be 2, representing start and end index of the bin
    if(bin_range != 2)
        throw std::runtime_error("Dimension of a bin range must be in 2D");
    for(int i=0; i<num_vars; i++) {
        int rg_start = *(ptr+i*bin_range+0);
        int rg_end = *(ptr+i*bin_range+1);
        std::pair<int, int> rg(rg_start, rg_end);
        c_binrange.push_back(rg);
    }

    PyJointHistogram chist(pyHist.hist->conditionalHist(c_orivars, c_condvars, c_binrange));
    return chist;
}