// $Id$
//
//  Copyright (C) 2013 Greg Landrum
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
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/MonomerInfo.h>

namespace python = boost::python;

namespace RDKit {
struct monomerinfo_wrapper {
  static void wrap() {
    std::string classDoc =
        "The class to store monomer information attached to Atoms\n";
    python::class_<AtomMonomerInfo>("AtomMonomerInfo", classDoc.c_str(),
                                    python::init<>(python::args("self")))
        .def(
            python::init<AtomMonomerInfo::AtomMonomerType, const std::string &,
                         const std::string &, int,
                         const std::string &, const std::string &>(
                (python::arg("self"), python::arg("type"),
                 python::arg("name") = "", python::arg("residueName") = "",
                 python::arg("resNum") = 0, python::arg("chainId") = "",
                 python::arg("monomerClass") = "")))
        .def("GetName", &AtomMonomerInfo::getName,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetMonomerType", &AtomMonomerInfo::getMonomerType,
             python::args("self"))
        .def("GetResidueName", &AtomMonomerInfo::getResidueName,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetResidueNumber", &AtomMonomerInfo::getResidueNumber,
                 python::args("self"))
        .def("GetChainId", &AtomMonomerInfo::getChainId,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetMonomerClass", &AtomMonomerInfo::getMonomerClass,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))

        .def("SetName", &AtomMonomerInfo::setName, python::args("self", "nm"))
        .def("SetMonomerType", &AtomMonomerInfo::setMonomerType,
             python::args("self", "typ"))
        .def("SetResidueName", &AtomMonomerInfo::setResidueName,
             python::args("self", "val"))
        .def("SetResidueNumber", &AtomMonomerInfo::setResidueNumber,
             python::args("self", "val"))
        .def("SetChainId", &AtomMonomerInfo::setChainId,
             python::args("self", "val"))
        .def("SetMonomerClass", &AtomMonomerInfo::setMonomerClass,
             python::args("self", "val"))
        ;

    python::enum_<AtomMonomerInfo::AtomMonomerType>("AtomMonomerType")
        .value("UNKNOWN", AtomMonomerInfo::UNKNOWN)
        .value("PDBRESIDUE", AtomMonomerInfo::PDBRESIDUE)
        .value("OTHER", AtomMonomerInfo::OTHER);

    classDoc = "The class to store PDB residue information attached to Atoms\n";
    python::class_<AtomPDBResidueInfo, python::bases<AtomMonomerInfo>>(
        "AtomPDBResidueInfo", classDoc.c_str(),
        python::init<>(python::args("self")))
        .def(python::init<std::string, int, std::string, std::string, int,
                          std::string, std::string, double, double, bool,
                          unsigned int, unsigned int, std::string>(
            (python::arg("self"), python::arg("atomName"),
             python::arg("serialNumber") = 1, python::arg("altLoc") = "",
             python::arg("residueName") = "", python::arg("residueNumber") = 0,
             python::arg("chainId") = "", python::arg("insertionCode") = "",
             python::arg("occupancy") = 1.0, python::arg("tempFactor") = 0.0,
             python::arg("isHeteroAtom") = false,
             python::arg("secondaryStructure") = 0,
             python::arg("segmentNumber") = 0,
             python::arg("monomerClass") = "")))
        .def("GetSerialNumber", &AtomPDBResidueInfo::getSerialNumber,
             python::args("self"))
        .def("GetAltLoc", &AtomPDBResidueInfo::getAltLoc,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetResidueName", &AtomPDBResidueInfo::getResidueName,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetResidueNumber", &AtomPDBResidueInfo::getResidueNumber,
             python::args("self"))
        .def("GetChainId", &AtomPDBResidueInfo::getChainId,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetInsertionCode", &AtomPDBResidueInfo::getInsertionCode,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))
        .def("GetOccupancy", &AtomPDBResidueInfo::getOccupancy,
             python::args("self"))
        .def("GetTempFactor", &AtomPDBResidueInfo::getTempFactor,
             python::args("self"))
        .def("GetIsHeteroAtom", &AtomPDBResidueInfo::getIsHeteroAtom,
             python::args("self"))
        .def("GetSecondaryStructure",
             &AtomPDBResidueInfo::getSecondaryStructure, python::args("self"))
        .def("GetSegmentNumber", &AtomPDBResidueInfo::getSegmentNumber,
             python::args("self"))
        .def("GetMonomerClass", &AtomPDBResidueInfo::getMonomerClass,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"))

        .def("SetSerialNumber", &AtomPDBResidueInfo::setSerialNumber,
             python::args("self", "val"))
        .def("SetAltLoc", &AtomPDBResidueInfo::setAltLoc,
             python::args("self", "val"))
        .def("SetResidueName", &AtomPDBResidueInfo::setResidueName,
             python::args("self", "val"))
        .def("SetResidueNumber", &AtomPDBResidueInfo::setResidueNumber,
             python::args("self", "val"))
        .def("SetChainId", &AtomPDBResidueInfo::setChainId,
             python::args("self", "val"))
        .def("SetInsertionCode", &AtomPDBResidueInfo::setInsertionCode,
             python::args("self", "val"))
        .def("SetOccupancy", &AtomPDBResidueInfo::setOccupancy,
             python::args("self", "val"))
        .def("SetTempFactor", &AtomPDBResidueInfo::setTempFactor,
             python::args("self", "val"))
        .def("SetIsHeteroAtom", &AtomPDBResidueInfo::setIsHeteroAtom,
             python::args("self", "val"))
        .def("SetSecondaryStructure",
             &AtomPDBResidueInfo::setSecondaryStructure,
             python::args("self", "val"))
        .def("SetSegmentNumber", &AtomPDBResidueInfo::setSegmentNumber,
             python::args("self", "val"))
        .def("SetMonomerClass", &AtomPDBResidueInfo::setMonomerClass,
             python::args("self", "val"))

        ;
  };
};
}  // namespace RDKit

void wrap_monomerinfo() { RDKit::monomerinfo_wrapper::wrap(); }
