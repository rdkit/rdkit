#include <boost/python.hpp>
#include "classA.h"
#include "classC.h"
#include <DataStructs/ExplicitBitVect.h>
namespace python = boost::python;

void testCrossA(classA *A) {
  A->printA();
};

void testCrossC(classC *C) {
  C->printC();
};

classC *getClassC() {
  classC *nc = new classC();
  return nc;
}

ExplicitBitVect *getEBV() {
  ExplicitBitVect *ebv = new  ExplicitBitVect(20);
  return ebv;
}

BOOST_PYTHON_MODULE(moduleB) 
{
  python::def("testCrossA", testCrossA);
  python::def("testCrossC", testCrossC);
  python::def("GetClassC", getClassC,
              python::return_value_policy<python::manage_new_object>());
  python::def("GetEBV", getEBV,
              python::return_value_policy<python::manage_new_object>());
}
