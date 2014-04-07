// $Id$
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
                            unsigned timeout,bool verbose){
    std::vector<ROMOL_SPTR> ms;
    unsigned int nElems=python::extract<unsigned int>(mols.attr("__len__")());
    ms.resize(nElems);
    for(unsigned int i=0;i<nElems;++i){
      ms[i] = python::extract<ROMOL_SPTR>(mols[i]);
    }
    MCSParameters *ps=new MCSParameters();
    ps->MaximizeBonds=maximizeBonds;
    ps->Threshold=threshold;
    ps->Timeout=timeout;
    ps->Verbose=verbose;
    MCSResult *res=new MCSResult(findMCS(ms,ps));
    delete ps;
    return res;
  }
}

namespace {
  struct mcsresult_wrapper {
    static void wrap(){
      python::class_<RDKit::MCSResult>("MCSResult","stores MCS results",python::no_init)
        .def_readonly("numAtoms",&RDKit::MCSResult::NumAtoms)
        .def_readonly("numBonds",&RDKit::MCSResult::NumBonds)
        .def_readonly("smartsString",&RDKit::MCSResult::SmartsString)
        .def_readonly("canceled",&RDKit::MCSResult::Canceled)
        ;
    }
  };
}

BOOST_PYTHON_MODULE(rdFMCS) {
  
  python::scope().attr("__doc__") =
    "Module containing a C++ implementation of the FMCS algorithm";
  mcsresult_wrapper::wrap();
   
  std::string docString = "Find the MCS for a set of molecules";
  python::def("FindMCS", RDKit::FindMCSWrapper,
              (python::arg("mols"),
               python::arg("maximizeBonds")=true,
               python::arg("threshold")=1.0,
               python::arg("timeout")=3600,
               python::arg("verbose")=false
               ),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());
}
