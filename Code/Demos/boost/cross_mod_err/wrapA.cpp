#include "classA.h"
#include <boost/python.hpp>

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4788) // warning: long & complicated stl warning
#pragma warning (disable: 4660)
#pragma warning (disable: 4275) // warning: non dll-interface class used as...
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float' 
#endif

namespace python = boost::python;

struct A_wrapper {
  static void wrap() {
    python::class_<classA>("classA", python::init<>())
      .def("printA", &classA::printA)
      ;
  };
};

void wrap_classA() {
  A_wrapper::wrap();
}

BOOST_PYTHON_MODULE(moduleA)
{
  wrap_classA();
}

