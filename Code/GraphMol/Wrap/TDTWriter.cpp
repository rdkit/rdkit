// $Id: TDTWriter.cpp 5061 2006-03-08 00:36:29Z glandrum $
//
//  Copyright (C) 2005-2006  Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

//ours
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include "rdchem.h"
#include <RDBoost/PySequenceHolder.h>

namespace python = boost::python;

namespace RDKit {
  void SetTDTWriterProps(TDTWriter &writer, python::object props) {
    // convert the python list to a STR_VECT
    STR_VECT propNames;
    PySequenceHolder<std::string> seq(props);
    for (unsigned int i = 0; i < seq.size(); i++) {
      propNames.push_back(seq[i]);
    }
    writer.setProps(propNames);
  }
  struct tdtwriter_wrap {
    static void wrap() {
      python::class_<TDTWriter>("TDTWriter",
				"A class for writing molecules to TDT files.\n",
				python::init<std::string>(python::args("fileName")))
	.def("SetProps", SetTDTWriterProps,
	     "Sets the properties to be written to the output file\n\n"
	     "  ARGUMENTS:\n\n"
	     "    - props: a list or tuple of property names\n\n")
	.def("write", &TDTWriter::write,
             (python::arg("self"), python::arg("mol"), python::arg("confId")=defaultConfId),
	     "Writes a molecule to the output file.\n\n"
	     "  ARGUMENTS:\n\n"
	     "    - mol: the Mol to be written\n"
             "    - confId: (optional) ID of the conformation to write\n\n")
	.def("flush", &TDTWriter::flush,
	     "Flushes the output file (forces the disk file to be updated).\n\n"
	     )
	.def("NumMols", &TDTWriter::numMols,
	     "Returns the number of molecules written so far.\n\n"
	     )
	.def("SetWrite2D",&TDTWriter::setWrite2D,(python::arg("self"),python::arg("state")=true),
	     "causes 2D conformations to be written (default is 3D conformations)")
	.def("GetWrite2D",&TDTWriter::getWrite2D)
	.def("SetWriteNames",&TDTWriter::setWriteNames,
	     (python::arg("self"),python::arg("state")=true),
	     "causes names to be written to the output file as NAME records")
	.def("GetWriteNames",&TDTWriter::getWriteNames)
	.def("SetNumDigits",&TDTWriter::setNumDigits,
	     "sets the number of digits to be written for coordinates")
	.def("GetNumDigits",&TDTWriter::getNumDigits)
	;
    };
  };
}

void wrap_tdtwriter() {
  RDKit::tdtwriter_wrap::wrap();
}
