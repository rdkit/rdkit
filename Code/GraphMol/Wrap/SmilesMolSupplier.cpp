// $Id$
//
//  Copyright (C) 2003-2008  Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <RDBoost/iterator_next.h>
#include <string>

//ours
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>
#include "MolSupplier.h"

namespace python = boost::python;

namespace RDKit {

  SmilesMolSupplier *SmilesSupplierFromText(std::string text,
				      std::string delimiter=" ",
				      int smilesColumn=0,
				      int nameColumn=1,
				      bool titleLine=true,
				      bool sanitize=true){
    SmilesMolSupplier *res=new SmilesMolSupplier();
    res->setData(text,delimiter,smilesColumn,
		 nameColumn,titleLine,
		 sanitize);
    return res;
  }


  std::string smilesMolSupplierClassDoc="A class which supplies molecules from a text file.\n \
\n \
  Usage examples:\n \
\n \
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n \
       >>> suppl = SmilesMolSupplier('in.smi')\n \
       >>> for mol in suppl:\n \
       ...    mol.GetNumAtoms()\n \
\n \
    2) Lazy evaluation 2:\n \
       >>> suppl = SmilesMolSupplier('in.smi')\n \
       >>> mol1 = suppl.next()\n \
       >>> mol2 = suppl.next()\n \
       >>> suppl.reset()\n \
       >>> mol3 = suppl.next()\n \
       # mol3 and mol1 are the same: \
       >>> MolToSmiles(mol3)==MolToSmiles(mol1)\n \
\n \
    3) Random Access:  all molecules are constructed as soon as we ask for the\n \
       length:\n \
       >>> suppl = SmilesMolSupplier('in.smi')\n \
       >>> nMols = len(suppl)\n \
       >>> for i in range(nMols):\n \
       ...   suppl[i].GetNumAtoms()\n \
\n \
  If the input file has a title line and more than two columns (smiles and id), the\n\
  additional columns will be used to set properties on each molecule.  The properties\n\
  are accessible using the mol.GetProp(propName) method.\n\
\n";

  std::string smsDocStr="Constructor\n \n \
  ARGUMENTS: \n \
\n \
    - fileName: name of the file to be read\n \
\n \
    - delimiter: (optional) text delimiter (a string).  Defauts to ' '.\n \
\n \
    - smilesColumn: (optional) index of the column containing the SMILES\n \
      data.  Defaults to 0.\n \
\n \
    - nameColumn: (optional) index of the column containing molecule names.\n \
      Defaults to 1.\n \
\n \
    - titleLine: (optional) set this toggle if the file contains a title line.\n \
      Defaults to 1.\n \
\n \
    - sanitize: (optional) toggles sanitization of molecules as they are read.\n \
      Defaults to 1.\n \
\n";
  struct smimolsup_wrap {
    static void wrap() {
      python::class_<SmilesMolSupplier,boost::noncopyable>("SmilesMolSupplier",
			  smilesMolSupplierClassDoc.c_str(),
			  python::init<std::string, std::string, int, int, bool, bool>(
					   (python::arg("data"),
					    python::arg("delimiter")=" ",
					    python::arg("smilesColumn")=0,
					    python::arg("nameColumn")=1,
					    python::arg("titleLine")=true,
					    python::arg("sanitize")=true),smsDocStr.c_str()))
	.def(python::init<>())
	.def("__iter__", (SmilesMolSupplier *(*)(SmilesMolSupplier *))&MolSupplIter,
	     python::return_internal_reference<1>() )
	.def(NEXT_METHOD, (ROMol *(*)(SmilesMolSupplier *))&MolSupplNext,
	     "Returns the next molecule in the file.  Raises _StopIteration_ on EOF.\n",
	     python::return_value_policy<python::manage_new_object>())
	.def("__getitem__", (ROMol *(*)(SmilesMolSupplier *,int))&MolSupplGetItem,
	     python::return_value_policy<python::manage_new_object>())
	.def("reset", &SmilesMolSupplier::reset,
	     "Resets our position in the file to the beginning.\n")
	.def("__len__", &SmilesMolSupplier::length)
	.def("SetData", &SmilesMolSupplier::setData,
	     "Sets the text to be parsed",
	     (python::arg("self"),python::arg("data"),python::arg("delimiter")=" ",
	      python::arg("smilesColumn")=0,python::arg("nameColumn")=1,python::arg("titleLine")=true,
	      python::arg("sanitize")=true))
	.def("GetItemText", &SmilesMolSupplier::getItemText,
	     "returns the text for an item",
	     (python::arg("self"),python::arg("index")))
	;

      python::def("SmilesMolSupplierFromText",SmilesSupplierFromText,
		  (python::arg("text"),python::arg("delimiter")=" ",
		   python::arg("smilesColumn")=0,python::arg("nameColumn")=1,python::arg("titleLine")=true,
		   python::arg("sanitize")=true),
		  python::return_value_policy<python::manage_new_object>());

    }
  };
}

void wrap_smisupplier() {
  RDKit::smimolsup_wrap::wrap();
}
