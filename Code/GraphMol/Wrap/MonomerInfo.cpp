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
#include <boost/python.hpp>
#include <string>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/MonomerInfo.h>

namespace python = boost::python;

namespace RDKit{
  struct monomerinfo_wrapper {
    static void wrap() {
      std::string classDoc="The class to store monomer information attached to Atoms\n";
      python::class_<AtomMonomerInfo>("AtomMonomerInfo",classDoc.c_str(),python::init<>())
        .def(python::init<AtomMonomerInfo::AtomMonomerType,const std::string &>((python::arg("type"),
                                                                                 python::arg("name")="")))
        .def("GetName", &AtomMonomerInfo::getName,
             python::return_value_policy<python::copy_const_reference>())
        .def("GetMonomerType", &AtomMonomerInfo::getMonomerType)
        .def("SetName", &AtomMonomerInfo::setName)
        .def("SetMonomerType", &AtomMonomerInfo::setMonomerType)
        ;

      python::enum_<AtomMonomerInfo::AtomMonomerType>("AtomMonomerType")
        .value("UNKNOWN",AtomMonomerInfo::UNKNOWN)
        .value("PDBRESIDUE",AtomMonomerInfo::PDBRESIDUE)
        .value("OTHER",AtomMonomerInfo::OTHER)
        ;

      classDoc="The class to store PDB residue information attached to Atoms\n";
      python::class_<AtomPDBResidueInfo,python::bases<AtomMonomerInfo> >("AtomPDBResidueInfo",classDoc.c_str(),python::init<>())
        .def(python::init<std::string,int,std::string,std::string,int,std::string,std::string,
                          double,double,bool,unsigned int,unsigned int>((python::arg("atomName"),
                                                            python::arg("serialNumber")=1,
                                                            python::arg("altLoc")="",
                                                            python::arg("residueName")="",
                                                            python::arg("residueNumber")=0,
                                                            python::arg("chainId")="",
                                                            python::arg("insertionCode")="",
                                                            python::arg("occupancy")=1.0,
                                                            python::arg("tempFactor")=0.0,
                                                            python::arg("isHeteroAtom")=false,
                                                            python::arg("secondaryStructure")=0,
                                                            python::arg("segmentNumber")=0)))
        .def("GetSerialNumber",&AtomPDBResidueInfo::getSerialNumber)
        .def("GetAltLoc",&AtomPDBResidueInfo::getAltLoc,python::return_value_policy<python::copy_const_reference>())
        .def("GetResidueName",&AtomPDBResidueInfo::getResidueName,python::return_value_policy<python::copy_const_reference>())
        .def("GetResidueNumber",&AtomPDBResidueInfo::getResidueNumber)
        .def("GetChainId",&AtomPDBResidueInfo::getChainId,python::return_value_policy<python::copy_const_reference>())
        .def("GetInsertionCode",&AtomPDBResidueInfo::getInsertionCode,python::return_value_policy<python::copy_const_reference>())
        .def("GetOccupancy",&AtomPDBResidueInfo::getOccupancy)
        .def("GetTempFactor",&AtomPDBResidueInfo::getTempFactor)
        .def("GetIsHeteroAtom",&AtomPDBResidueInfo::getIsHeteroAtom)
        .def("GetSecondaryStructure",&AtomPDBResidueInfo::getSecondaryStructure)
        .def("GetSegmentNumber",&AtomPDBResidueInfo::getSegmentNumber)

        .def("SetSerialNumber",&AtomPDBResidueInfo::setSerialNumber)
        .def("SetAltLoc",&AtomPDBResidueInfo::setAltLoc)
        .def("SetResidueName",&AtomPDBResidueInfo::setResidueName)
        .def("SetResidueNumber",&AtomPDBResidueInfo::setResidueNumber)
        .def("SetChainId",&AtomPDBResidueInfo::setChainId)
        .def("SetInsertionCode",&AtomPDBResidueInfo::setInsertionCode)
        .def("SetOccupancy",&AtomPDBResidueInfo::setOccupancy)
        .def("SetTempFactor",&AtomPDBResidueInfo::setTempFactor)
        .def("SetIsHeteroAtom",&AtomPDBResidueInfo::setIsHeteroAtom)
        .def("SetSecondaryStructure",&AtomPDBResidueInfo::setSecondaryStructure)
        .def("SetSegmentNumber",&AtomPDBResidueInfo::setSegmentNumber)

        ;

    };
  };
}

void wrap_monomerinfo() {
  RDKit::monomerinfo_wrapper::wrap();
}
      

      
