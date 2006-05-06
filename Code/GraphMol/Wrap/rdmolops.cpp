// $Id: rdmolops.cpp 4978 2006-02-18 00:59:33Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "rdmolops.h"
#include <boost/python.hpp>

#include "Numeric/arrayobject.h"

#include <RDGeneral/types.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/SanitException.h>

namespace python = boost::python;
using namespace RDKit;

void rdSanitExceptionTranslator(RDKit::MolSanitizeException const& x){
  std::ostringstream ss;
  ss << "Sanitization error: " << x.message();
  PyErr_SetString(PyExc_ValueError,ss.str().c_str());
}


void wrap_molops();


BOOST_PYTHON_MODULE(rdmolops)
{

  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for manipulating and querying molecules."
    ;
  import_array();
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  python::register_exception_translator<RDKit::MolSanitizeException>(&rdSanitExceptionTranslator);

  // ******************************
  // Functions from MolOps
  //****************************
  wrap_molops();

}


