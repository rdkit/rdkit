#include <boost/python.hpp>
#include "classA.h"
#include "classC.h"

namespace python = boost::python;

void testCrossA(classA *A) {
  A->printA();
};

void testCrossC(classC *C) {
  C->printC();
};


BOOST_PYTHON_MODULE(moduleB) 
{
  python::def("testCrossA", testCrossA);
  python::def("testCrossC", testCrossC);
  
}
