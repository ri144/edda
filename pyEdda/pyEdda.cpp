#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "py_distribution_modeler.h"

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
    
    //m.def("tf", &test_f);
    //m.def("createDM", &createDM);
    //m.def("test_array", &test_array);
    //m.def("multi_ret", &multi_ret);


}
