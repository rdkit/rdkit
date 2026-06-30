//
//  Copyright (C) 2024-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/MolProcessing/MolProcessing.h>
#include <GraphMol/FileParsers/GeneralFileReader.h>
#include <RDBoost/Wrap.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {
template <typename OutputType>
nb::tuple getFingerprintsHelper(
    const std::string &fileName, nb::object pyGenerator,
    const GeneralMolSupplier::SupplierOptions &options) {
  FingerprintGenerator<OutputType> *generator = nullptr;
  if (!pyGenerator.is_none()) {
    try {
      generator = nb::cast<FingerprintGenerator<OutputType> *>(pyGenerator);
    } catch (const nb::cast_error &) {
      throw nb::next_overload();
    }
  }

  std::vector<std::unique_ptr<ExplicitBitVect>> fps;
  {
    NOGIL gil;
    fps = MolProcessing::getFingerprintsForMolsInFile(fileName, options,
                                                      generator);
  }
  nb::list pyFingerprints;
  for (auto &fp : fps) {
    pyFingerprints.append(
        nb::cast(fp.release(), nb::rv_policy::take_ownership));
  }
  return nb::tuple(pyFingerprints);
}
}  // namespace

NB_MODULE(rdMolProcessing, m) {
  m.doc() = "Module containing functions for working with groups of molecules";

  nb::class_<GeneralMolSupplier::SupplierOptions>(m, "SupplierOptions",
                                                  "Supplier Options")
      .def(nb::init<>())
      .def_rw("numThreads",
              &GeneralMolSupplier::SupplierOptions::numWriterThreads,
              "the number of threads to use while working")
      .def_rw("sanitize", &GeneralMolSupplier::SupplierOptions::sanitize)
      .def_rw("removeHs", &GeneralMolSupplier::SupplierOptions::removeHs)
      .def_rw("strictParsing",
              &GeneralMolSupplier::SupplierOptions::strictParsing)
      .def_rw("delimiter", &GeneralMolSupplier::SupplierOptions::delimiter,
              "used for SMILES files")
      .def_rw("smilesColumn",
              &GeneralMolSupplier::SupplierOptions::smilesColumn,
              "used for SMILES files")
      .def_rw("nameColumn", &GeneralMolSupplier::SupplierOptions::nameColumn,
              "used for SMILES files")
      .def_rw("titleLine", &GeneralMolSupplier::SupplierOptions::titleLine,
              "used for SMILES files")
      .def_rw("nameRecord", &GeneralMolSupplier::SupplierOptions::nameRecord,
              "used for TDT files")
      .def_rw("confId2D", &GeneralMolSupplier::SupplierOptions::confId2D,
              "used for TDT files")
      .def_rw("confId3D", &GeneralMolSupplier::SupplierOptions::confId3D,
              "used for TDT files");

  m.def("GetFingerprintsForMolsInFile", getFingerprintsHelper<std::uint32_t>,
        "filename"_a, "generator"_a = nb::none(),
        "options"_a = GeneralMolSupplier::SupplierOptions(),
        "returns the fingerprints for the molecules in a file (32 bit version)");
  m.def("GetFingerprintsForMolsInFile", getFingerprintsHelper<std::uint64_t>,
        "filename"_a, "generator"_a = nb::none(),
        "options"_a = GeneralMolSupplier::SupplierOptions(),
        "returns the fingerprints for the molecules in a file (64 bit version)");
}
