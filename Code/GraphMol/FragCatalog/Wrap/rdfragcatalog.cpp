// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "rdfragcatalog.h"
#include <boost/python.hpp>

namespace python = boost::python;

void wrap_fragcat();
void wrap_fragparams();
void wrap_fragcatgen() ;
void wrap_fragFPgen() ;

BOOST_PYTHON_MODULE(rdfragcatalog)
{
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  wrap_fragcat();
  wrap_fragparams();
  wrap_fragcatgen();
  wrap_fragFPgen() ;
}
