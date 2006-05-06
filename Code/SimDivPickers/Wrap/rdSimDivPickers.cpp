// $Id: rdSimDivPickers.cpp 4959 2006-02-17 23:53:31Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define PY_ARRAY_UNIQUE_SYMBOL rdpicker_array_API
#include <RDBoost/Wrap.h>
#include "Numeric/arrayobject.h"
#include <RDGeneral/types.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <vector>
namespace python = boost::python;
using namespace RDPickers;

void wrap_maxminpick();
void wrap_HierarchCP();

BOOST_PYTHON_MODULE(rdSimDivPickers)
{
  python::scope().attr("__doc__") =
    "Module containing the diversity and similarity pickers"
    ;

  import_array();
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);


  wrap_maxminpick();
  wrap_HierarchCP();
}

