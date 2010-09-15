#include "classC.h"
#include <boost/python.hpp>

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4788) // warning: long & complicated stl warning
#pragma warning (disable: 4660)
#pragma warning (disable: 4275) // warning: non dll-interface class used as...
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float' 
#endif

namespace python = boost::python;

struct C_wrapper {
  static void wrap() {
    python::class_<classC>("classC", python::init<>())
      .def("printC", &classC::printC)
      ;
  };
};

void wrap_classC() {
  C_wrapper::wrap();
}

BOOST_PYTHON_MODULE(moduleC)
{
  wrap_classC();
}

