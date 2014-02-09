// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <RDGeneral/types.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/SanitException.h>

namespace python = boost::python;
using namespace RDKit;

void wrap_queries();


BOOST_PYTHON_MODULE(rdqueries)
{
  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for querying molecules."
    ;
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  wrap_queries();
}


