//
//  Copyright (c) 2011-2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <boost/python.hpp>
#include "../inchi.h"

namespace {
  boost::python::tuple MolToInchi(const RDKit::ROMol &mol, std::string options)
  {
    RDKit::ExtraInchiReturnValues rv;
    const char* _options = nullptr;
    if (options.size()) {
      _options = options.c_str();
    }
    std::string inchi = RDKit::MolToInchi(mol, rv, _options);
    return boost::python::make_tuple(inchi, rv.returnCode, rv.messagePtr, rv.logPtr,
                                     rv.auxInfoPtr);
  }

  boost::python::tuple MolBlockToInchi(const std::string &molblock, std::string options)
  {
    RDKit::ExtraInchiReturnValues rv;
    const char* _options = nullptr;
    if (options.size()) {
      _options = options.c_str();
    }
    std::string inchi = RDKit::MolBlockToInchi(molblock, rv, _options);
    return boost::python::make_tuple(inchi, rv.returnCode, rv.messagePtr, rv.logPtr,
                                     rv.auxInfoPtr);
  }
  boost::python::tuple InchiToMol(const std::string &inchi, bool sanitize,
                                  bool removeHs)
  {
    RDKit::ExtraInchiReturnValues rv;
    RDKit::ROMol* mol = RDKit::InchiToMol(inchi, rv, sanitize, removeHs);
    if (mol == nullptr) {
      return boost::python::make_tuple(boost::python::object(), rv.returnCode,
                                       rv.messagePtr, rv.logPtr);
    } else {
      return boost::python::make_tuple(RDKit::ROMOL_SPTR(mol), rv.returnCode,
                                       rv.messagePtr, rv.logPtr);
    }
  }
}

BOOST_PYTHON_MODULE(rdinchi) {
  std::string docString = "return a ROMol for a InChI string\n\
  Returns:\n\
    a tuple with:\n\
      - the molecule\n\
      - the return code from the InChI conversion\n\
      - a string with any messages from the InChI conversion\n\
      - a string with any log messages from the InChI conversion\n";
  boost::python::def("InchiToMol", InchiToMol,
                     (boost::python::arg("inchi"),
                      boost::python::arg("sanitize")=true,
                      boost::python::arg("removeHs")=true),
                     docString.c_str()
                     );
  docString = "return the InChI for a ROMol molecule.\n\
\n\
  Arguments:\n\
    - mol: the molecule to use.\n\
    - options: the InChI generation options.\n\
      Options should be prefixed with either a - or a /\n\
      Available options are explained in the InChI technical FAQ:\n\
      http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14\n\
      and the User Guide:\n\
      http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf\n\
  Returns:\n\
    a tuple with:\n\
      - the InChI\n\
      - the return code from the InChI conversion\n\
      - a string with any messages from the InChI conversion\n\
      - a string with any log messages from the InChI conversion\n\
      - a string with the InChI AuxInfo\n";
  boost::python::def("MolToInchi", MolToInchi,
                     (boost::python::arg("mol"),
                      boost::python::arg("options")=std::string()),
                     docString.c_str()
                     );
  docString = "return the InChI for a ROMol molecule.\n\
\n\
  Arguments:\n\
    - molblock: the mol block to use.\n\
    - options: the InChI generation options.\n\
      Options should be prefixed with either a - or a /\n\
      Available options are explained in the InChI technical FAQ:\n\
      http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14\n\
      and the User Guide:\n\
      http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf\n\
  Returns:\n\
    a tuple with:\n\
      - the InChI\n\
      - the return code from the InChI conversion\n\
      - a string with any messages from the InChI conversion\n\
      - a string with any log messages from the InChI conversion\n\
      - a string with the InChI AuxInfo\n";
  boost::python::def("MolBlockToInchi", MolBlockToInchi,
                     (boost::python::arg("molblock"),
                      boost::python::arg("options")=std::string()),
                     docString.c_str()
                     );
  docString = "return the InChI key for an InChI string";
  boost::python::def("InchiToInchiKey", RDKit::InchiToInchiKey,
                     (boost::python::arg("inchi")),
                     docString.c_str()
                    );
  docString = "return the InChI key for a ROMol molecule.\n\
\n\
  Arguments:\n\
    - mol: the molecule to use.\n\
    - options: the InChI generation options.\n\
      Options should be prefixed with either a - or a /\n\
      Available options are explained in the InChI technical FAQ:\n\
      https://www.inchi-trust.org/technical-faq-2/#15.14\n\
      and the User Guide available from:\n\
      https://www.inchi-trust.org/downloads/\n\
  Returns: the InChI key\n";
  boost::python::def("MolToInchiKey", RDKit::MolToInchiKey,
                     (boost::python::arg("mol"),
                      boost::python::arg("options")=""),
                     docString.c_str()
                     );

}
