// $Id: SDMolSupplier.cpp 585 2008-03-30 13:36:56Z glandrum $
//
//  Copyright (C) 2009 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <string>

#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <RDStreams/streams.h>

// ours
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Wrap/ContextManagers.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/Wrap.h>

#include "MolSupplier.h"

namespace python = boost::python;

namespace RDKit {

// ForwardSDMolSupplier cannot (yet?) be reset, so we have to override
// the template that was defined in MolSupplier.h.
// Note that this returns a pointer to the supplier itself, so be careful
// that it doesn't get deleted by python!
template <>
ForwardSDMolSupplier *MolSupplIter(ForwardSDMolSupplier *suppl) {
  return suppl;
}

ROMol *MolSupplNext(ForwardSDMolSupplier *suppl) {
  ROMol *res = 0;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (...) {
      res = 0;
    }
  }
  if (!res && suppl->atEnd()) {
    PyErr_SetString(PyExc_StopIteration, "End of supplier hit");
    throw boost::python::error_already_set();
  }
  return res;
}

ForwardSDMolSupplier *createForwardSupplier(std::string filename, bool sanitize,
                                            bool removeHs) {
  std::vector<std::string> splitName;
  boost::split(splitName, filename, boost::is_any_of("."));
  const std::string &ext = splitName.back();
  std::istream *strm = nullptr;
  if (ext == "sdf") {
    strm = new std::ifstream(filename.c_str(), std::ios_base::binary);
  } else if (ext == "gz") {
    strm = new gzstream(filename);
  } else if (ext == "bz2") {
#ifdef RDK_USE_BZIP2
    strm = new bz2stream(filename);
#else
    throw_value_error("bzip2 support not enabled");
#endif
  } else {
    throw_value_error("Unrecognized extension: " + ext);
  }
  if (!strm->good()) {
    delete strm;
    throw_value_error("could not open file: " + filename);
  }
  return new ForwardSDMolSupplier(strm, true, sanitize, removeHs);
}

std::string csdMolSupplierClassDoc =
    "A class which supplies molecules from an SD file.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n\n\
       >>> suppl = SDMolSupplier('in.smi')\n\
       >>> for mol in suppl:\n\
       ...    mol.GetNumAtoms()\n\
\n\
  Properties in the SD file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
struct compressedsdmolsup_wrap {
  static void wrap() {
    python::class_<ForwardSDMolSupplier, boost::noncopyable>(
        "_CompressedSDMolSupplier", csdMolSupplierClassDoc.c_str(),
        python::no_init)
        .def("__iter__", &MolSupplIter<ForwardSDMolSupplier>,
             python::return_internal_reference<1>())
        .def("__enter__", &MolIOEnter<ForwardSDMolSupplier>,
             python::return_internal_reference<>())
        .def("__exit__", &MolIOExit<ForwardSDMolSupplier>)
        .def("__next__", &MolSupplNext<ForwardSDMolSupplier>,
             "Returns the next molecule in the file.  Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>());
    python::def("CompressedSDMolSupplier", createForwardSupplier,
                (python::arg("fileName"), python::arg("sanitize") = true,
                 python::arg("removeHs") = true),
                python::return_value_policy<python::manage_new_object>());
  };
};
}  // namespace RDKit

void wrap_compressedsdsupplier() { RDKit::compressedsdmolsup_wrap::wrap(); }
