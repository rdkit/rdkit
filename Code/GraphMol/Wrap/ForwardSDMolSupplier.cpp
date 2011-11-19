// $Id: SDMolSupplier.cpp 1625 2011-01-13 04:22:56Z glandrum $
//
//  Copyright (C) 2003-2010  Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

//ours
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/python_streambuf.h>

#include "MolSupplier.h"

namespace python = boost::python;

  using boost_adaptbx::python::streambuf;
namespace {

  class LocalForwardSDMolSupplier : public RDKit::ForwardSDMolSupplier
  {
  public:
    LocalForwardSDMolSupplier(python::object &input,
                              bool sanitize,bool removeHs){
      // FIX: minor leak here
      streambuf *sb=new streambuf(input);
      dp_inStream=new streambuf::istream(*sb);
      df_owner=true;
      df_sanitize=sanitize;
      df_removeHs=removeHs;
      POSTCONDITION(dp_inStream,"bad instream");
    }
  };

  LocalForwardSDMolSupplier *FwdMolSupplIter(LocalForwardSDMolSupplier *self){
    return self;
  }
}

namespace RDKit {

  std::string fsdMolSupplierClassDoc="A class which supplies molecules from file-like object containing SD data.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n\
       >>> suppl = ForwardSDMolSupplier(file('in.sdf'))\n\
       >>> for mol in suppl:\n\
       ...    if mol is not None: mol.GetNumAtoms()\n\
\n\
    2) we can also read from compressed files: \n\
       >>> import gzip\n\
       >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))\n\
       >>> for mol in suppl:\n \
       ...   if mol is not None: print mol.GetNumAtoms()\n\
\n\
  Properties in the SD file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
  struct forwardsdmolsup_wrap {
    static void wrap() {
      python::class_<LocalForwardSDMolSupplier,
        boost::noncopyable>("ForwardSDMolSupplier",
                            fsdMolSupplierClassDoc.c_str(),
                            python::init<python::object &,bool,bool>
                            ((python::arg("fileobj"),
                              python::arg("sanitize")=true,
                              python::arg("removeHs")=true))
                            [python::with_custodian_and_ward_postcall<0,2>()])
        .def("next", (ROMol *(*)(LocalForwardSDMolSupplier *))&MolSupplNext,
	     "Returns the next molecule in the file.  Raises _StopIteration_ on EOF.\n",
	     python::return_value_policy<python::manage_new_object>())
	.def("atEnd", &ForwardSDMolSupplier::atEnd,
	     "Returns whether or not we have hit EOF.\n")
	.def("__iter__", &FwdMolSupplIter,
	     python::return_internal_reference<1>() )
	;
    };
  };
}

void wrap_forwardsdsupplier() {
  RDKit::forwardsdmolsup_wrap::wrap();
}
