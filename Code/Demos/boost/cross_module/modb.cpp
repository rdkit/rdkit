//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include "mods.h"

namespace python = boost::python;

BOOST_PYTHON_MODULE(modb)
{
    python::class_<ClassB>("ClassB")
      .def("Get3",&ClassB::get3)
      .def("ReturnOther",&ClassB::returnOther,
	   python::return_value_policy<python::manage_new_object>())
      .def("AcceptOther",&ClassB::acceptOther)
      ;

}
