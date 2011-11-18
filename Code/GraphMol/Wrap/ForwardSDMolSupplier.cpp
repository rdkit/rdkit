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

namespace {
  using boost_adaptbx::python::streambuf;
  RDKit::ForwardSDMolSupplier *createForwardMolSupplier(streambuf& input,
                                                 bool sanitize, bool removeHs){
    streambuf::istream *is = new streambuf::istream(input);
    return new RDKit::ForwardSDMolSupplier(is,true,sanitize,removeHs);
  }

  struct python_streambuf_wrapper
  {
    typedef boost_adaptbx::python::streambuf wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<wt, boost::noncopyable>("streambuf", no_init)
        .def(init<object&, std::size_t>((
          arg("python_file_obj"),
          arg("buffer_size")=0),"documentation")[with_custodian_and_ward<1,2>()])
        .def_readwrite(
          "default_buffer_size", wt::default_buffer_size,
          "The default size of the buffer sitting "
          "between a Python file object and a C++ stream.")
      ;
    }
  };

  struct python_ostream_wrapper
  {
    typedef boost_adaptbx::python::ostream wt;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<std::ostream, boost::noncopyable>("std_ostream", no_init);
      class_<wt, boost::noncopyable, bases<std::ostream> >("ostream", no_init)
        .def(init<object&, std::size_t>((
          arg("python_file_obj"),
          arg("buffer_size")=0)))
      ;
    }
  };

  RDKit::ForwardSDMolSupplier *FwdMolSupplIter(RDKit::ForwardSDMolSupplier *self){
    return self;
  }
}

namespace RDKit {

  std::string fsdMolSupplierClassDoc="A class which supplies molecules from an SD file.\n \
\n \
  Usage examples:\n \
\n \
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n \
       >>> suppl = SDMolSupplier('in.sdf')\n \
       >>> for mol in suppl:\n \
       ...    mol.GetNumAtoms()\n \
\n \
    2) Lazy evaluation 2:\n \
       >>> suppl = SDMolSupplier('in.sdf')\n \
       >>> mol1 = suppl.next()\n \
       >>> mol2 = suppl.next()\n \
       >>> suppl.reset()\n \
       >>> mol3 = suppl.next()\n \
       # mol3 and mol1 are the same: \n \
       >>> MolToSmiles(mol3)==MolToSmiles(mol1)\n \
\n \
    3) Random Access:\n \
       >>> suppl = SDMolSupplier('in.sdf')\n \
       >>> mol1 = suppl[0] \n \
       >>> mol2 = suppl[1] \n \
       NOTE: this will generate an IndexError if the supplier doesn't have that many\n \
       molecules.\n \
\n \
    4) Random Access 2:  looping over all molecules \n \
       >>> suppl = SDMolSupplier('in.sdf')\n \
       >>> nMols = len(suppl)\n \
       >>> for i in range(nMols):\n \
       ...   suppl[i].GetNumAtoms()\n \
\n \
  Properties in the SD file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
  struct forwardsdmolsup_wrap {
    static void wrap() {
      python::class_<ForwardSDMolSupplier,boost::noncopyable>("ForwardSDMolSupplier",
						       fsdMolSupplierClassDoc.c_str(),
						       python::no_init)
	.def("next", (ROMol *(*)(ForwardSDMolSupplier *))&MolSupplNext,
	     "Returns the next molecule in the file.  Raises _StopIteration_ on EOF.\n",
	     python::return_value_policy<python::manage_new_object>())
	.def("__iter__", &FwdMolSupplIter,
	     python::return_internal_reference<1>() )
	;
      python::def("CreateForwardSDMolSupplier",createForwardMolSupplier,
                  (python::arg("stream"),
                   python::arg("sanitize")=true,
                   python::arg("removeHs")=true),
                  "",
                  python::return_value_policy<python::manage_new_object>());

      python_streambuf_wrapper::wrap();
      python_ostream_wrapper::wrap();
    };
  };
}

void wrap_forwardsdsupplier() {
  RDKit::forwardsdmolsup_wrap::wrap();
}
