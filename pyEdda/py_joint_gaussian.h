#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "distributions/joint_gaussian.h"

using namespace std;
using namespace edda::dist;

namespace py = pybind11;

class PyJointGaussian {
public:
	JointGaussian* hist;

public:
    /// constructors
    PyJointGaussian(JointGaussian* pyHist) {
        hist = pyHist;
    }

    PyJointGaussian(PyJointGaussian* pyHist) {
        hist = pyHist->hist;
    }

    PyJointGaussian(py::array_t<Real> mean, py::array_t<Real> covs) {
        py::buffer_info info = mean.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("The mean of a PyJointGaussian should be a 1D numpy array");
        int dim = info.shape[0];
        ublas_vector c_mean(dim);
        for(int i=0; i<dim; i++) {
            c_mean(i) = *(ptr+i);
        }

        // convert covariance matrix
        info = covs.request();
        ptr = static_cast<Real *>(info.ptr);
        ndim = info.ndim;
        if(ndim != 2)
            throw std::runtime_error("The covariance matrix of a PyJointGaussian should be a 2D numpy array");
        int dim0 = info.shape[0];
        int dim1 = info.shape[1];
        if(dim0 != dim1)
            throw std::runtime_error("The covariance matrix should be a square matrix");

        ublas_matrix c_covs(dim0, dim1);
        for(int i=0; i<dim0; i++) {
            for(int j=0; j<dim1; j++) {
                c_covs(i, j) = *(ptr + i*dim1+j);
            }
        }

        hist = new JointGaussian(c_mean, c_covs);
    }


    /// seter and geter
    void setMatrices(py::array_t<Real> mat) {
        py::buffer_info info = mat.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 2)
            throw std::runtime_error("The covariance matrix of a PyJointGaussian should be a 2D numpy array");
        int dim0 = info.shape[0];
        int dim1 = info.shape[1];
        ublas_matrix c_covs(dim0, dim1);
        for(int i=0; i<dim0; i++) {
            for(int j=0; j<dim1; j++) {
                c_covs(i, j) = *(ptr + i*dim1+j);
            }
        }
        hist->setMatrices(c_covs);
    }

    py::array_t<Real> getJointSample() {
        std::vector<Real> v = hist->getJointSample();
        return py::array_t<Real>(v.size(), v.data());
    }

    Real getJointLogPdf(py::array_t<Real> x) {
        py::buffer_info info = x.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions should be 1");
        int dim = info.shape[0];
        std::vector<Real> c_x {ptr, ptr + dim};
        
        return hist->getJointLogPdf(c_x);
    }

    Real getJointPdf(py::array_t<Real> x) {
        py::buffer_info info = x.request();
        auto ptr = static_cast<Real *>(info.ptr);
        int ndim = info.ndim;
        if(ndim != 1)
            throw std::runtime_error("Number of dimensions should be 1");
        int dim = info.shape[0];
        std::vector<Real> c_x {ptr, ptr + dim};
        
        return hist->getJointPdf(c_x);
    }

    py::array_t<Real> getMean() {
        ublas_vector m = hist->getMean();
        std::vector<Real> v{m.begin(), m.end()};
        return py::array_t<Real>(v.size(), v.data());
    }
};


py::array_t<Real> getJointMean_Gaussian(PyJointGaussian& pyhist) {
    return pyhist.getMean();
}

py::array_t<Real> getJointSample_Gaussian(PyJointGaussian& pyhist) {
    return pyhist.getJointSample();
}

Real getJointPdf_Gaussian(PyJointGaussian& pyhist, py::array_t<Real> x) {
    // py::buffer_info info = x.request();
    // auto ptr = static_cast<Real *>(info.ptr);
    // int ndim = info.ndim;
    // if(ndim != 1)
    //     throw std::runtime_error("Number of dimensions for the sample x must be one");
    // int dim = info.shape[0];
    // std::vector<Real> c_x {ptr, ptr + dim};

    // return pyhist.hist->getJointPdf(c_x);
    return pyhist.getJointPdf(x);   // a simpler way
}
