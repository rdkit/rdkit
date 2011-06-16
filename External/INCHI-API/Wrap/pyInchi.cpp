// $Id$
//
// Copyright (C) 2011-2011 Novartis Institutes for BioMedical Research, Inc
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include "../inchi.h"

namespace {
  boost::python::tuple MolToInchi(const RDKit::ROMol &mol, std::string options)
  {
    std::string message, log, auxInfo;
    RDKit::ExtraInchiReturnValues rv;
    rv.messagePtr = &message;
    rv.logPtr = &log;
    rv.auxInfoPtr = &auxInfo;
    const char* _options = NULL;
    if (options.size())
      _options = options.c_str();
    std::string inchi = RDKit::MolToInchi(mol, rv, _options);
    return boost::python::make_tuple(inchi, rv.returnCode, message, log,
                                     auxInfo);
  }

  boost::python::tuple InchiToMol(const std::string &inchi, bool sanitize,
                                  bool removeHs)
  {
    std::string message, log;
    RDKit::ExtraInchiReturnValues rv;
    rv.messagePtr = &message;
    rv.logPtr = &log;
    RDKit::ROMol* mol = RDKit::InchiToMol(inchi, rv, sanitize, removeHs);
    if (mol == NULL)
      return boost::python::make_tuple(boost::python::object(), rv.returnCode,
                                       message, log);
    else {
      return boost::python::make_tuple(RDKit::ROMOL_SPTR(mol), rv.returnCode,
                                       message, log);
    }
  }
}

BOOST_PYTHON_MODULE(rdinchi) {
  std::string docString = "return a ROMol for a InChI string";
  boost::python::def("InchiToMol", InchiToMol,
                     (boost::python::arg("inchi"),
                      boost::python::arg("sanitize")=true,
                      boost::python::arg("removeHs")=true),
                     docString.c_str()
                     );
  docString = "return the InChI for a ROMol molecule. If options is an empty"
    "string, standard InChI string is returned";
  boost::python::def("MolToInchi", MolToInchi,
                     (boost::python::arg("mol"),
                      boost::python::arg("options")=std::string()),
                     docString.c_str()
                     );
  docString = "return the InChI key for an InChI string";
  boost::python::def("InchiToInchiKey", RDKit::InchiToInchiKey,
                     (boost::python::arg("inchi")),
                     docString.c_str()
                    );
}

