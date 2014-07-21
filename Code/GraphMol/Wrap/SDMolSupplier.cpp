// $Id$
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
#include <RDBoost/iterator_next.h>

#include "MolSupplier.h"

namespace python = boost::python;

namespace RDKit {
  void setStreamIndices(SDMolSupplier &self,python::object arg){
    std::vector<std::streampos> loc;
    PySequenceHolder<int> seq(arg);
    loc.reserve(seq.size());
    for(unsigned int i=0;i<seq.size();++i){
      loc.push_back(static_cast<std::streampos>(seq[i]));
    }
    self.setStreamIndices(loc);
  }


  std::string sdMolSupplierClassDoc="A class which supplies molecules from an SD file.\n \
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
  struct sdmolsup_wrap {
    static void wrap() {
      python::class_<SDMolSupplier,boost::noncopyable>("SDMolSupplier",
						       sdMolSupplierClassDoc.c_str(),
						       python::init<>())
	.def(python::init<std::string,bool,bool,bool>((python::arg("fileName"),
                                                       python::arg("sanitize")=true,
                                                       python::arg("removeHs")=true,
                                                       python::arg("strictParsing")=true)))
	.def("__iter__", (SDMolSupplier *(*)(SDMolSupplier *))&MolSupplIter,
	     python::return_internal_reference<1>() )
	.def(NEXT_METHOD, (ROMol *(*)(SDMolSupplier *))&MolSupplNextAcceptNullLastMolecule,
	     "Returns the next molecule in the file.  Raises _StopIteration_ on EOF.\n",
	     python::return_value_policy<python::manage_new_object>())
	.def("__getitem__", (ROMol *(*)(SDMolSupplier *,int))&MolSupplGetItem,
	     python::return_value_policy<python::manage_new_object>())
	.def("reset", &SDMolSupplier::reset,
	     "Resets our position in the file to the beginning.\n")
	.def("__len__", &SDMolSupplier::length)
	.def("SetData", &SDMolSupplier::setData,
	     "Sets the text to be parsed",
	     (python::arg("self"),python::arg("data"),python::arg("sanitize")=true,
              python::arg("removeHs")=true))
	.def("_SetStreamIndices", setStreamIndices,
	     "Sets the locations of mol beginnings in the input stream. Be *very* careful with this method.",
	     (python::arg("self"),python::arg("locs")))
	.def("GetItemText", &SDMolSupplier::getItemText,
	     "returns the text for an item",
	     (python::arg("self"),python::arg("index")))
	.def("atEnd", &SDMolSupplier::atEnd,
	     "Returns whether or not we have hit EOF.\n")
	;
    };
  };
}

void wrap_sdsupplier() {
  RDKit::sdmolsup_wrap::wrap();
}
