// $Id: DataStructs.cpp 4946 2006-02-17 01:44:04Z glandrum $
//
// Copyright (c) 2001-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//

#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include "DataStructs.h"

namespace python = boost::python;

void wrap_SBV();
void wrap_EBV();
void wrap_BitOps();
void wrap_Utils();
void wrap_discreteValVect();


BOOST_PYTHON_MODULE(cDataStructs)
{
  python::scope().attr("__doc__") =
    "Module containing an assortment of functionality for basic data structures.\n"
    "\n"
    "At the moment the data structures defined are:\n"
    "  Bit Vector classes (for storing signatures, fingerprints and the like:\n"
    "    - ExplicitBitVect:  class for relatively small (10s of thousands of bits) or\n"
    "                        dense bit vectors.\n"
    "    - SparseBitVect:  class for large, sparse bit vectors\n"
    ;
  
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  wrap_Utils();
  wrap_SBV();
  wrap_EBV();
  wrap_BitOps();
  wrap_discreteValVect();
}
