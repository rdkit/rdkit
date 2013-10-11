/* 
* $Id$
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
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Bond.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
%}

%ignore RDKit::RWMol::addAtom(Atom *atom);
%ignore RDKit::RWMol::addAtom(Atom *atom,bool updateLabel);
%ignore RDKit::RWMol::addBond(Atom *beginAtom, Atom *endAtom, Bond::BondType order);
%ignore RDKit::RWMol::addBond(Atom *beginAtom, Atom *endAtom);
%ignore RDKit::RWMol::addBond(Bond *bond);

%shared_ptr(RDKit::RWMol)
%include "enums.swg"
#if swifjava
%javaconst(1);
#endif
%include <GraphMol/FileParsers/FileParsers.h>
%include <GraphMol/RWMol.h>

%extend RDKit::RWMol {  
  static RDKit::RWMOL_SPTR MolFromSmiles(std::string smi,int debugParse=0,bool sanitize=1,
                                         std::map<std::string,std::string> *replacements=0){
    return RDKit::RWMOL_SPTR(RDKit::SmilesToMol(smi, debugParse, sanitize,replacements));
  }
  static RDKit::RWMOL_SPTR MolFromSmarts(std::string sma,int debugParse=0,bool mergeHs=false,
                                         std::map<std::string,std::string> *replacements=0){
    return RDKit::RWMOL_SPTR(RDKit::SmartsToMol(sma, debugParse, mergeHs,replacements));
  }
static RDKit::RWMOL_SPTR MolFromMolBlock(std::string molB,
                                  bool sanitize=true,bool removeHs=true){
  RDKit::RWMol *mol=0;
    mol=RDKit::MolBlockToMol(molB,sanitize,removeHs);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMolFile(std::string filename,
                                 bool sanitize=true,bool removeHs=true){
  RDKit::RWMol *mol=0;
    mol=RDKit::MolFileToMol(filename,sanitize,removeHs);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromTPLFile(std::string fName,bool sanitize=true,
                      bool skipFirstConf=false) {
  RDKit::RWMol *mol=0;
    mol=RDKit::TPLFileToMol(fName, sanitize, skipFirstConf);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMol2File(std::string fName,bool sanitize=true,bool removeHs=true,
                       RDKit::Mol2Type variant=RDKit::CORINA) {
  RDKit::RWMol *mol=0;
  mol=RDKit::Mol2FileToMol(fName, sanitize, removeHs, variant);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMol2Block(const std::string &molBlock,bool sanitize=true,bool removeHs=true,
                        RDKit::Mol2Type variant=RDKit::CORINA) {
  RDKit::RWMol *mol=0;
    mol=RDKit::Mol2BlockToMol(molBlock, sanitize, removeHs, variant);
  return RDKit::RWMOL_SPTR(mol);
}

static RDKit::RWMOL_SPTR MolFromPDBBlock(std::string molB,
                                         bool sanitize=true,bool removeHs=true,
                                         unsigned int flavor=0){
  RDKit::RWMol *mol=0;
  mol=RDKit::PDBBlockToMol(molB,sanitize,removeHs,flavor);
  return RDKit::RWMOL_SPTR(mol);
}

static RDKit::RWMOL_SPTR MolFromPDBFile(std::string fName,
                                        bool sanitize=true,bool removeHs=true,
                                        unsigned int flavor=0){
  RDKit::RWMol *mol=0;
  mol=RDKit::PDBFileToMol(fName,sanitize,removeHs,flavor);
  return RDKit::RWMOL_SPTR(mol);
}

 
  /* Methods from MolFileStereoChem.h */
  void DetectAtomStereoChemistry(const RDKit::Conformer *conf) {
	RDKit::DetectAtomStereoChemistry(*($self), conf);
  };

   /* From Kekulize.cpp, MolOps.h */
    void Kekulize(bool markAtomsBonds=true, unsigned int maxBackTracks=100) {
	RDKit::MolOps::Kekulize(*($self), markAtomsBonds, maxBackTracks);
   }

  /* MolOps.h */
  void sanitizeMol() {
	RDKit::MolOps::sanitizeMol(*($self));
  }

}
