#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "py_distribution_modeler.h"
#include "py_histogram.h"
#include "py_gaussian.h"
#include "py_gmm.h"
#include "py_joint_histogram.h"
#include "py_joint_gaussian.h"

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


    //JointHistogram
    py::class_<PyJointHistogram>(m, "JointHistogram")
        .def(py::init<PyJointHistogram*>())
        .def(py::init<py::array_t<Real>, int, py::array_t<Real>, py::array_t<Real>, py::array_t<int>>())
        .def("setMinVals", &PyJointHistogram::setMinVals)
        .def("getMinVals", &PyJointHistogram::getMinVals)
        .def("setMaxVals", &PyJointHistogram::setMaxVals)
        .def("getMaxVals", &PyJointHistogram::getMaxVals)
        .def("setBinWidths", &PyJointHistogram::setBinWidths)
        .def("getBinWidths", &PyJointHistogram::getBinWidths)
        .def("setNumBins", &PyJointHistogram::setNumBins)
        .def("setNumVars", &PyJointHistogram::setNumVars)
        .def("getNumBins", &PyJointHistogram::getNumBins)
        .def("getNumVars", &PyJointHistogram::getNumVars)
        .def("getJointMean", &PyJointHistogram::getJointMean)
        .def("getJointSample", &PyJointHistogram::getJointSample);

    // C-like functions for JointHistogram
    m.def("getJointMean", &getJointMean_py);
    m.def("marginalization", &marginalization, py::return_value_policy::reference);
    m.def("conditionalHist", &conditionalHist, py::return_value_policy::reference);


    //JointGaussian
    py::class_<PyJointGaussian>(m, "JointGaussian")
        .def(py::init<PyJointGaussian*>())
        .def(py::init<py::array_t<Real>, py::array_t<Real>>())
        .def("getMean", &PyJointGaussian::getMean)
        .def("getJointSample", &PyJointGaussian::getJointSample)
        .def("getJointPdf", &PyJointGaussian::getJointPdf);

    // C-like functions for JointGaussian
    m.def("getJointMean_Gaussian", &getJointMean_Gaussian); // TODO: need to find a way for function overloading in pybind11
    m.def("getJointSample_Gaussian", &getJointSample_Gaussian); // TODO: need to find a way for function overloading in pybind11
    m.def("getJointPdf_Gaussian", &getJointPdf_Gaussian); // TODO: need to find a way for function overloading in pybind11
}
