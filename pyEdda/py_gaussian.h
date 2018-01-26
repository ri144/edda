#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "distributions/gaussian.h"

using namespace std;
using namespace edda::dist;

namespace py = pybind11;

class PyGaussian {
protected:
	Gaussian* gaussian;

public:
    PyGaussian()
    {
        gaussian = new Gaussian();
    }

    PyGaussian(Real m, Real s)
    {
        gaussian = new Gaussian(m, s);
    }

    double getGaussianMean()
    {
        return getMean(*gaussian);
    }

    double getGaussianVar()
    {
        return getVar(*gaussian);
    }

    double getGaussianPdf(const double x)
    {
        return getPdf(*gaussian, x);
    }

    double getGaussianSample()
    {
        return getSample(*gaussian);
    }

    double getGaussianCdf(const double x)
    {
        return getCdf(*gaussian, x);
    }

    double getGaussianCdfPrecise(const double x)
    {
        return getCdfPrecise(*gaussian, x);
    }

    void output()
    {
        std::cout << *gaussian << std::endl;
    }
};