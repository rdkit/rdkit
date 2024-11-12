//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/MolProcessing/MolProcessing.h>
#include <GraphMol/FileParsers/GeneralFileReader.h>

namespace python = boost::python;
using namespace RDKit;

namespace {
template <typename OutputType>
python::tuple getFingerprintsHelper(
    const std::string &fileName, python::object pyGenerator,
    const GeneralMolSupplier::SupplierOptions &options) {
  FingerprintGenerator<OutputType> *generator = nullptr;
  if (pyGenerator) {
    generator =
        python::extract<FingerprintGenerator<OutputType> *>(pyGenerator);
  }

  std::vector<std::unique_ptr<ExplicitBitVect>> fps;
  {
    NOGIL gil;
    fps = MolProcessing::getFingerprintsForMolsInFile(fileName, options,
                                                      generator);
  }
  python::list pyFingerprints;
  boost::python::manage_new_object::apply<ExplicitBitVect *>::type converter;
  for (auto &fp : fps) {
    // transfer ownership to python
    python::handle<> handle(converter(fp.release()));
    pyFingerprints.append(handle);
  }

  return python::tuple(pyFingerprints);
}
}  // namespace

BOOST_PYTHON_MODULE(rdMolProcessing) {
  python::scope().attr("__doc__") =
      "Module containing functions for working with groups of molecules";

  python::class_<GeneralMolSupplier::SupplierOptions>(
      "SupplierOptions", "Supplier Options", python::init<>())
      .def_readwrite("numThreads",
                     &GeneralMolSupplier::SupplierOptions::numWriterThreads,
                     "the number of threads to use while working")
      .def_readwrite("sanitize", &GeneralMolSupplier::SupplierOptions::sanitize)
      .def_readwrite("removeHs", &GeneralMolSupplier::SupplierOptions::removeHs)
      .def_readwrite("strictParsing",
                     &GeneralMolSupplier::SupplierOptions::strictParsing)
      .def_readwrite("delimiter",
                     &GeneralMolSupplier::SupplierOptions::delimiter,
                     "used for SMILES files")
      .def_readwrite("smilesColumn",
                     &GeneralMolSupplier::SupplierOptions::smilesColumn,
                     "used for SMILES files")
      .def_readwrite("nameColumn",
                     &GeneralMolSupplier::SupplierOptions::nameColumn,
                     "used for SMILES files")
      .def_readwrite("titleLine",
                     &GeneralMolSupplier::SupplierOptions::titleLine,
                     "used for SMILES files")
      .def_readwrite("nameRecord",
                     &GeneralMolSupplier::SupplierOptions::nameRecord,
                     "used for TDT files")
      .def_readwrite("confId2D", &GeneralMolSupplier::SupplierOptions::confId2D,
                     "used for TDT files")
      .def_readwrite("confId3D", &GeneralMolSupplier::SupplierOptions::confId3D,
                     "used for TDT files");

  python::def(
      "GetFingerprintsForMolsInFile",
      (python::tuple(*)(const std::string &, python::object,
                        const GeneralMolSupplier::SupplierOptions &))
          getFingerprintsHelper<std::uint32_t>,
      (python::arg("filename"), python::arg("generator") = python::object(),
       python::arg("options") = GeneralMolSupplier::SupplierOptions()),
      "returns the fingerprints for the molecules in a file (32 bit version)");
  python::def(
      "GetFingerprintsForMolsInFile",
      (python::tuple(*)(const std::string &, python::object,
                        const GeneralMolSupplier::SupplierOptions &))
          getFingerprintsHelper<std::uint64_t>,
      (python::arg("filename"), python::arg("generator") = python::object(),
       python::arg("options") = GeneralMolSupplier::SupplierOptions()),
      "returns the fingerprints for the molecules in a file (64 bit version)");
}
