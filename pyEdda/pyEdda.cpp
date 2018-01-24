#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "py_distribution_modeler.h"
#include "py_histogram.h"
#include "py_gaussian.h"
#include "py_gmm.h"

using namespace std;


namespace py = pybind11;


PYBIND11_MODULE(pyedda, m) {
    py::class_<PyDistributionModeler>(m, "DistributionModeler")
       .def(py::init<int >())
       .def("computeGMM", &PyDistributionModeler::computeGMM)
       .def("computeHistogram", &PyDistributionModeler::computeHistogram)
       .def("printDistr", &PyDistributionModeler::printDistr)
       .def("getDistr", &PyDistributionModeler::getDistr)
       .def("unit_test1", []() {});
    
    //Univariate Histogram
    py::class_<PyHistogram>(m, "Histogram")
    .def(py::init<PyHistogram*>())
    .def(py::init<py::array_t<Real>, const int>())
    .def(py::init<py::array_t<Real>, const int, const Real, const Real >())
    .def("getMean", &PyHistogram::getMean)
    .def("getVar", &PyHistogram::getVar)
    .def("getPdf", &PyHistogram::getPdf)
    .def("getCdf", &PyHistogram::getCdf)
    .def("getSample", &PyHistogram::getSample)
    .def("output", &PyHistogram::output)
    .def("getBins", &PyHistogram::getBins)
    .def("getMaxValue", &PyHistogram::getMaxValue)
    .def("getMinValue", &PyHistogram::getMinValue)
    .def("getBinValue", &PyHistogram::getBinValue);

    //Univariate Gaussian
    py::class_<PyGaussian>(m, "Gaussian")
    .def(py::init<>())
    .def(py::init<Real, Real>())
    .def("getMean", &PyGaussian::getGaussianMean)
    .def("getVar", &PyGaussian::getGaussianVar)
    .def("getPdf", &PyGaussian::getGaussianPdf)
    .def("getSample", &PyGaussian::getGaussianSample)
    .def("getCdf", &PyGaussian::getGaussianCdf)
    .def("getCdfPrecise", &PyGaussian::getGaussianCdfPrecise)
    .def("output", &PyGaussian::output);

    //Univariate GMM
    py::class_<PyGMM>(m, "GMM")
    .def(py::init<>())
    .def(py::init<int>())
    .def(py::init<py::array_t<Real>>())
    .def(py::init<py::array_t<Real>, int>())
    .def("assign", &PyGMM::assign)
    .def("getNumComponents", &PyGMM::getNumComponents)
    .def("normalizeWeights", &PyGMM::normalizeWeights)
    .def("getMean", &PyGMM::getGMMMean)
    .def("getVar", &PyGMM::getGMMVar)
    .def("getPdf", &PyGMM::getGMMPdf)
    .def("getSample", &PyGMM::getGMMSample)
    .def("getCdf", &PyGMM::getGMMCdf)
    .def("output", &PyGMM::output);
}
