//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/MonomerInfo.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
struct monomerinfo_wrapper {
  static void wrap(nb::module_ &m) {
    std::string classDoc =
        R"DOC(The class to store monomer information attached to Atoms
)DOC";
    nb::class_<AtomMonomerInfo>(m, "AtomMonomerInfo", classDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<AtomMonomerInfo::AtomMonomerType, const std::string &,
                      const std::string &, int, const std::string &,
                      const std::string &>(),
             "type"_a, "name"_a = "", "residueName"_a = "", "resNum"_a = 0,
             "chainId"_a = "", "monomerClass"_a = "")
        .def("GetName", &AtomMonomerInfo::getName)
        .def("GetMonomerType", &AtomMonomerInfo::getMonomerType)
        .def("GetResidueName", &AtomMonomerInfo::getResidueName)
        .def("GetResidueNumber", &AtomMonomerInfo::getResidueNumber)
        .def("GetChainId", &AtomMonomerInfo::getChainId)
        .def("GetMonomerClass", &AtomMonomerInfo::getMonomerClass)

        .def("SetName", &AtomMonomerInfo::setName, "nm"_a)
        .def("SetMonomerType", &AtomMonomerInfo::setMonomerType, "typ"_a)
        .def("SetResidueName", &AtomMonomerInfo::setResidueName, "val"_a)
        .def("SetResidueNumber", &AtomMonomerInfo::setResidueNumber, "val"_a)
        .def("SetChainId", &AtomMonomerInfo::setChainId, "val"_a)
        .def("SetMonomerClass", &AtomMonomerInfo::setMonomerClass, "val"_a);

    nb::enum_<AtomMonomerInfo::AtomMonomerType>(m, "AtomMonomerType")
        .value("UNKNOWN", AtomMonomerInfo::UNKNOWN)
        .value("PDBRESIDUE", AtomMonomerInfo::PDBRESIDUE)
        .value("OTHER", AtomMonomerInfo::OTHER)
        .export_values();

    classDoc =
        R"DOC(The class to store PDB residue information attached to Atoms
)DOC";
    nb::class_<AtomPDBResidueInfo, AtomMonomerInfo>(m, "AtomPDBResidueInfo",
                                                    classDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<std::string, int, std::string, std::string, int,
                      std::string, std::string, double, double, bool,
                      unsigned int, unsigned int, std::string>(),
             "atomName"_a, "serialNumber"_a = 1, "altLoc"_a = "",
             "residueName"_a = "", "residueNumber"_a = 0, "chainId"_a = "",
             "insertionCode"_a = "", "occupancy"_a = 1.0, "tempFactor"_a = 0.0,
             "isHeteroAtom"_a = false, "secondaryStructure"_a = 0,
             "segmentNumber"_a = 0, "monomerClass"_a = "")
        .def("GetSerialNumber", &AtomPDBResidueInfo::getSerialNumber)
        .def("GetAltLoc", &AtomPDBResidueInfo::getAltLoc)
        .def("GetResidueName", &AtomPDBResidueInfo::getResidueName)
        .def("GetResidueNumber", &AtomPDBResidueInfo::getResidueNumber)
        .def("GetChainId", &AtomPDBResidueInfo::getChainId)
        .def("GetInsertionCode", &AtomPDBResidueInfo::getInsertionCode)
        .def("GetOccupancy", &AtomPDBResidueInfo::getOccupancy)
        .def("GetTempFactor", &AtomPDBResidueInfo::getTempFactor)
        .def("GetIsHeteroAtom", &AtomPDBResidueInfo::getIsHeteroAtom)
        .def("GetSecondaryStructure",
             &AtomPDBResidueInfo::getSecondaryStructure)
        .def("GetSegmentNumber", &AtomPDBResidueInfo::getSegmentNumber)
        .def("GetMonomerClass", &AtomPDBResidueInfo::getMonomerClass)

        .def("SetSerialNumber", &AtomPDBResidueInfo::setSerialNumber, "val"_a)
        .def("SetAltLoc", &AtomPDBResidueInfo::setAltLoc, "val"_a)
        .def("SetResidueName", &AtomPDBResidueInfo::setResidueName, "val"_a)
        .def("SetResidueNumber", &AtomPDBResidueInfo::setResidueNumber, "val"_a)
        .def("SetChainId", &AtomPDBResidueInfo::setChainId, "val"_a)
        .def("SetInsertionCode", &AtomPDBResidueInfo::setInsertionCode, "val"_a)
        .def("SetOccupancy", &AtomPDBResidueInfo::setOccupancy, "val"_a)
        .def("SetTempFactor", &AtomPDBResidueInfo::setTempFactor, "val"_a)
        .def("SetIsHeteroAtom", &AtomPDBResidueInfo::setIsHeteroAtom, "val"_a)
        .def("SetSecondaryStructure",
             &AtomPDBResidueInfo::setSecondaryStructure, "val"_a)
        .def("SetSegmentNumber", &AtomPDBResidueInfo::setSegmentNumber, "val"_a)
        .def("SetMonomerClass", &AtomPDBResidueInfo::setMonomerClass, "val"_a);
  };
};
}  // namespace RDKit

void wrap_monomerinfo(nb::module_ &m) { RDKit::monomerinfo_wrapper::wrap(m); }
