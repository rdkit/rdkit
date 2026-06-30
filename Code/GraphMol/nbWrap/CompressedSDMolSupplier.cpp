//
//  Copyright (C) 2009-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/algorithm/string.hpp>

namespace io = boost::iostreams;

// ours
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include "ContextManagers.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
// ForwardSDMolSupplier cannot be reset, so __iter__ just returns self
ForwardSDMolSupplier *ForwardMolSupplIter(ForwardSDMolSupplier *suppl) {
  return suppl;
}

ROMol *ForwardMolSupplNext(ForwardSDMolSupplier *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  }
  if (!res && suppl->atEnd()) {
    throw nb::stop_iteration();
  }
  return res;
}

}  // namespace

ForwardSDMolSupplier *createForwardSupplier(std::string filename, bool sanitize,
                                            bool removeHs) {
  std::vector<std::string> splitName;
  boost::split(splitName, filename, boost::is_any_of("."));
  std::unique_ptr<io::filtering_istream> strm(new io::filtering_istream());
  if (splitName.back() == "sdf") {
  } else if (splitName.back() == "gz") {
#ifndef RDK_NOGZIP
    strm->push(io::gzip_decompressor());
#else
    throw nb::value_error("gzip support not enabled");
#endif
  } else if (splitName.back() == "bz2") {
#ifndef RDK_NOBZIP2
    strm->push(io::bzip2_decompressor());
#else
    throw nb::value_error("bzip2 support not enabled");
#endif
  } else {
    std::string errorTxt = "Unrecognized extension: " + splitName.back();
    throw nb::value_error(errorTxt.c_str());
  }
  io::file_source fileSource(filename);
  if (!fileSource.is_open()) {
    std::string errorTxt = "could not open file: " + filename;
    throw nb::value_error(errorTxt.c_str());
  }
  strm->push(fileSource);

  return new ForwardSDMolSupplier(strm.release(), true, sanitize, removeHs);
}

std::string csdMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from a compressed SD file.

  Usage examples:

    1) Lazy evaluation: the molecules are not constructed until we ask for them:

       >>> suppl = CompressedSDMolSupplier('in.sdf.gz')
       >>> for mol in suppl:
       ...    mol.GetNumAtoms()

  Properties in the SD file are used to set properties on each molecule.
  The properties are accessible using the mol.GetProp(propName) method.
)DOC";

struct compressedsdmolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<ForwardSDMolSupplier>(m, "_CompressedSDMolSupplier",
                                     csdMolSupplierClassDoc.c_str())
        .def("__iter__", &ForwardMolSupplIter,
             nb::rv_policy::reference_internal)
        .def("__enter__", &MolIOEnter<ForwardSDMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<ForwardSDMolSupplier>,
             "excType"_a = nb::none(), "excValue"_a = nb::none(),
             "traceback"_a = nb::none())
        .def(
            "__next__", &ForwardMolSupplNext, nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC");
    m.def("CompressedSDMolSupplier", createForwardSupplier, "fileName"_a,
          "sanitize"_a = true, "removeHs"_a = true,
          nb::rv_policy::take_ownership);
  };
};
}  // namespace RDKit

void wrap_compressedsdsupplier(nb::module_ &m) {
  RDKit::compressedsdmolsup_wrap::wrap(m);
}
