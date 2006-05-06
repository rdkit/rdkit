//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include <boost/python.hpp>
namespace python = boost::python;


// ------
//
// This one is relatively easy:
//   expose a function with a default argument
//
// ------
int func2(int v1,int plus=3);
int func2(int v1,int plus) {return v1+plus;};
BOOST_PYTHON_FUNCTION_OVERLOADS(f2_overloads, func2, 1, 2)


// ------
//
// More complex:
//   expose a templated function with a default argument
//
// ------
template <typename T>
T func(T v1,int plus=3);

template <typename T>
T func(T v1,int plus) {
  return v1+plus;
}

int (*f1_int)(int,int=3)=func; // gotta love that syntax!
BOOST_PYTHON_FUNCTION_OVERLOADS(f1_int_overloads, f1_int, 1, 2)
float (*f1_float)(float,int=3)=func;
BOOST_PYTHON_FUNCTION_OVERLOADS(f1_float_overloads, f1_float, 1, 2)



BOOST_PYTHON_MODULE(overloads)
{
  python::def("f2",func2,f2_overloads(python::args("v1","plus")));

  python::def("f1",f1_int,f1_int_overloads(python::args("v1","plus")));
  python::def("f1",f1_float,f1_float_overloads(python::args("v1","plus")));

}
