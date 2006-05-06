//
//  Copyright (C) 2003 Rational Discovery LLC
//

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
namespace python = boost::python;



// ----------
//
//  In both cases here, the restriction on the object passed in is
//  solely that it support the functions used.
//
// ----------
int seq_len(python::object seq){
  return python::len(seq);
}

int sum_first2(python::object seq){
  int sum;
  sum = python::extract<int>(seq[0]) + python::extract<int>(seq[1]);
  
  return sum;
}



BOOST_PYTHON_MODULE(python_objs)
{
  python::def("seq_len",seq_len);
  python::def("sum_first2",sum_first2);
}
