// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include "rdchem.h"
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/python_streambuf.h>

namespace python = boost::python;

namespace RDKit {
  using boost_adaptbx::python::streambuf;
  SmilesWriter *getSmilesWriter(python::object &fileobj,
                                std::string delimiter=" ",
                                std::string nameHeader="Name",
                                bool includeHeader=true,
                                bool isomericSmiles=false,
                                bool kekuleSmiles=false){
    // FIX: minor leak here
    streambuf *sb=new streambuf(fileobj);
    streambuf::ostream *ost=new streambuf::ostream(*sb);
    return new SmilesWriter(ost,delimiter,nameHeader,includeHeader,true,isomericSmiles,kekuleSmiles);
  }

  void SetSmiWriterProps(SmilesWriter &writer, python::object props) {
    // convert the python list to a STR_VECT
    STR_VECT propNames;

    PySequenceHolder<std::string> seq(props);
    for (unsigned int i = 0; i < seq.size(); i++) {
      propNames.push_back(seq[i]);
    }
    writer.setProps(propNames);
  }
  std::string swDocStr="Constructor.\n\n"
  "   ARGUMENTS:\n\n"
  "     - fileName: name of the output file. ('-' to write to stdout)\n"
  "     - delimiter: (optional) delimiter to be used to separate entries on each line.\n"
  "     - nameHeader: (optional) text to use for the name column in the header line.\n"
  "                   If this is blank, names will not be included in the output.\n"
  "     - includeHeader: (optional) toggles inclusion of a header line in the output file.\n"
  "     - isomericSmiles: (optional) toggles output of isomeric smiles (includes stereochem information).\n"
  "     - kekuleSmiles: (optional) toggles output of kekule smiles (no aromatic bonds for molecules that have been kekulized).\n\n";
  struct smiwriter_wrap {
    static void wrap() {
      python::class_<SmilesWriter>("SmilesWriter",
				   "A class for writing molecules to text files.",
                                   python::no_init)
        .def("__init__",python::make_constructor(&getSmilesWriter,python::default_call_policies(),
                                                 (python::arg("fileObj"),
                                                  python::arg("delimiter")=" ",
                                                  python::arg("nameHeader")="Name",
                                                  python::arg("includeHeader")=true,
                                                  python::arg("isomericSmiles")=false,
                                                  python::arg("kekuleSmiles")=false)))
        .def(python::init<std::string,std::string,std::string,bool,bool,bool>((python::arg("fileName"),
                                                                               python::arg("delimiter")=" ",
                                                                               python::arg("nameHeader")="Name",
                                                                               python::arg("includeHeader")=true,
                                                                               python::arg("isomericSmiles")=false,
                                                                               python::arg("kekuleSmiles")=false),
                                                                              swDocStr.c_str()))
	.def("SetProps", SetSmiWriterProps,
	     "Sets the properties to be written to the output file\n\n"
	     "  ARGUMENTS:\n\n"
	     "    - props: a list or tuple of property names\n\n" )
	.def("write", &SmilesWriter::write,
             (python::arg("self"), python::arg("mol"), python::arg("confId")=-1),
	     "Writes a molecule to the output file.\n\n"
	     "  ARGUMENTS:\n\n"
	     "    - mol: the Mol to be written\n"
	     "    - confId: (optional) ignored \n\n")
	.def("flush", &SmilesWriter::flush,
	     "Flushes the output file (forces the disk file to be updated).\n\n"
	     )
	.def("close", &SmilesWriter::close,
	     "Flushes the output file and closes it. The Writer cannot be used after this.\n\n"
	     )
	.def("NumMols", &SmilesWriter::numMols,
	     "Returns the number of molecules written so far.\n\n"
	     )
	;
    };
  };
}

void wrap_smiwriter() {
  RDKit::smiwriter_wrap::wrap();
}
