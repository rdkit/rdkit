// $Id$
//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior
//       written permission.
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
// Created by Greg Landrum, September 2006
//
#include <boost/python.hpp>
#include <GraphMol/SLNParse/SLNParse.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>

namespace python = boost::python;


void rdSLNParseExceptionTranslator(RDKit::SLNParseException const& x){
  std::ostringstream ss;
  ss << "SLNParseException: " << x.message();
  PyErr_SetString(PyExc_ValueError,ss.str().c_str());
}

namespace RDKit {
  ROMol *MolFromSLN(std::string sln,bool sanitize=1,bool debugParser=false){
    RWMol *newM = SLNToMol(sln,sanitize,debugParser);
    return static_cast<ROMol *>(newM);
  }
  ROMol *MolFromQuerySLN(std::string sln,bool mergeHs=1,bool debugParser=false){
    RWMol *newM = SLNQueryToMol(sln,mergeHs,debugParser);
    return static_cast<ROMol *>(newM);
  }
}

BOOST_PYTHON_MODULE(rdSLNParse) {
  python::scope().attr("__doc__") =
    "Module containing classes and functions for working with Sybyl line notation (SLN)."
    ;

  python::register_exception_translator<RDKit::SLNParseException>(&rdSLNParseExceptionTranslator);
    
  std::string docString;
   
  docString="Construct a molecule from an SLN string.\n\n\
    ARGUMENTS:\n\
\n\
    - SLN: the SLN string\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n\
  NOTE: the SLN should not contain query information or properties. To build a\n\
    query from SLN, use MolFromQuerySLN.\n\
\n";  
  python::def("MolFromSLN",RDKit::MolFromSLN,
        (python::arg("SLN"),
         python::arg("sanitize")=true,
         python::arg("debugParser")=false),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());
  
  docString="Construct a query molecule from an SLN string.\n\n\
  ARGUMENTS:\n\
\n\
    - SLN: the SLN string\n\
\n\
    - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached\n\
      heavy atoms. Defaults to False.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object suitable for using in substructure queries, None on failure.\n\
\n";  
  python::def("MolFromQuerySLN",RDKit::MolFromQuerySLN,
        (python::arg("SLN"),
         python::arg("mergeHs")=true,
         python::arg("debugParser")=false),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());
              
}
              
