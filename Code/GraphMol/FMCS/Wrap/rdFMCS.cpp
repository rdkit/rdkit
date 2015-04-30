//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/FMCS/FMCS.h>

namespace python = boost::python;

namespace RDKit {
   

    MCSResult *FindMCSWrapper(python::object mols,bool maximizeBonds,double threshold,
                              unsigned timeout,bool verbose,
                              bool matchValences,
                              bool ringMatchesRingOnly,bool completeRingsOnly,
                              bool matchChiralTag,
                              AtomComparator atomComp, BondComparator bondComp) {
        std::vector<ROMOL_SPTR> ms;
        unsigned int nElems=python::extract<unsigned int>(mols.attr("__len__")());
        ms.resize(nElems);
        for(unsigned int i=0; i<nElems; ++i) {
            if(!mols[i]) throw_value_error("molecule is None");
            ms[i] = python::extract<ROMOL_SPTR>(mols[i]);
        }

        MCSResult *res= new MCSResult(findMCS(ms,maximizeBonds,threshold,timeout,verbose,
                                               matchValences,ringMatchesRingOnly,completeRingsOnly,matchChiralTag,
                                               atomComp,bondComp));
        return res;
    }
}

namespace {
    struct mcsresult_wrapper {
        static void wrap() {
            python::class_<RDKit::MCSResult>("MCSResult","used to return MCS results",python::no_init)
              .def_readonly("numAtoms",&RDKit::MCSResult::NumAtoms,"number of atoms in MCS")
              .def_readonly("numBonds",&RDKit::MCSResult::NumBonds,"number of bonds in MCS")
              .def_readonly("smartsString",&RDKit::MCSResult::SmartsString,"SMARTS string for the MCS")
              .def_readonly("canceled",&RDKit::MCSResult::Canceled,"if True, the MCS calculation did not finish")
            ;
        }
    };
}

BOOST_PYTHON_MODULE(rdFMCS) {
    python::scope().attr("__doc__") =
        "Module containing a C++ implementation of the FMCS algorithm";
    mcsresult_wrapper::wrap();

    python::enum_<RDKit::AtomComparator>("AtomCompare")
    .value("CompareAny",RDKit::AtomCompareAny)
    .value("CompareElements",RDKit::AtomCompareElements)
    .value("CompareIsotopes",RDKit::AtomCompareIsotopes)
    ;
    python::enum_<RDKit::BondComparator>("BondCompare")
    .value("CompareAny",RDKit::BondCompareAny)
    .value("CompareOrder",RDKit::BondCompareOrder)
    .value("CompareOrderExact",RDKit::BondCompareOrderExact)
    ;

    std::string docString = "Find the MCS for a set of molecules";
    python::def("FindMCS", RDKit::FindMCSWrapper,
                (python::arg("mols"),
                 python::arg("maximizeBonds")=true,
                 python::arg("threshold")=1.0,
                 python::arg("timeout")=3600,
                 python::arg("verbose")=false,
                 python::arg("matchValences")=false,
                 python::arg("ringMatchesRingOnly")=false,
                 python::arg("completeRingsOnly")=false,
                 python::arg("matchChiralTag")=false,
                 python::arg("atomCompare")=RDKit::AtomCompareElements,
                 python::arg("bondCompare")=RDKit::BondCompareOrder
                ),
                python::return_value_policy<python::manage_new_object>(),
                docString.c_str());

}
