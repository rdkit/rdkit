//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include "mods.h"

namespace python = boost::python;

BOOST_PYTHON_MODULE(moda)
{
    python::class_<ClassA>("ClassA")
      .def("Get4",&ClassA::get4)
      ;

}
