/*
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following
*       disclaimer in the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc.
*       nor the names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%{
#include <GraphMol/MolOps.h>
%}

%template(MolSanitizeException_Vect) std::vector<boost::shared_ptr<RDKit::MolSanitizeException>>;

%newobject RDKit::MolOps::renumberAtoms;
%newobject RDKit::MolOps::removeHs;
%newobject RDKit::MolOps::addHs;
%newobject RDKit::MolOps::removeAllHs;
%newobject RDKit::MolOps::mergeQueryHs;
%newobject RDKit::MolOps::adjustQueryProperties;

%ignore RDKit::MolOps::detectChemistryProblems;
%include <GraphMol/MolOps.h>
%ignore RDKit::MolOps::sanitizeMol(RWMol &,unsigned int &,unsigned int &);

%inline %{
  int sanitizeMol(RDKit::RWMol &mol,int sanitizeOps){
    unsigned int opThatFailed;
    try{
      RDKit::MolOps::sanitizeMol(mol,opThatFailed,
                                 static_cast<unsigned int>(sanitizeOps));
    } catch(...) {

    }
    return static_cast<int>(opThatFailed);
  };
  std::vector<boost::shared_ptr<RDKit::MolSanitizeException>> detectChemistryProblems(RDKit::ROMol &mol,int sanitizeOps=RDKit::MolOps::SANITIZE_ALL){
    std::vector<boost::shared_ptr<RDKit::MolSanitizeException>> res;
    auto probs = RDKit::MolOps::detectChemistryProblems(mol,sanitizeOps);
    for(const auto &exc_ptr : probs) {
      res.push_back(boost::shared_ptr<RDKit::MolSanitizeException>(exc_ptr->copy()));
    }
    return res;
  };
%}
