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
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/Bond.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
%}

%template(RWMol_Vect) std::vector< boost::shared_ptr<RDKit::RWMol> >;

// ignore the methods that allow the molecule to take ownership of atoms/Bonds
// (instead of copying them). This just leads to memory problems with Java
%ignore RDKit::RWMol::addAtom(Atom *atom,bool updateLabel,bool takeOwnership);
%ignore RDKit::RWMol::addBond(Bond *bond,bool takeOwnership);

%newobject RDKit::SmilesToMol;
%newobject RDKit::SmartsToMol;
%newobject RDKit::MolBlockToMol;
%newobject RDKit::MolFileToMol;
%newobject RDKit::MolFromMolFile;
%newobject RDKit::MolFromTPLFIle;
%newobject RDKit::MolFromMol2File;
%newobject RDKit::MolFromMol2Block;
%newobject RDKit::MolFromPDBBlock;
%newobject RDKit::MolFromPDBFile;
%newobject RDKit::MolFromSequence;
%newobject RDKit::MolFromFasta;


%shared_ptr(RDKit::RWMol)
%include "enums.swg"
#if swifjava
%javaconst(1);
#endif
%include <GraphMol/FileParsers/FileParsers.h>
%include <GraphMol/SmilesParse/SmilesParse.h>
%include <GraphMol/RWMol.h>

%extend RDKit::RWMol {
  static RDKit::RWMOL_SPTR MolFromSmiles(const std::string &smi,int debugParse=0,bool sanitize=1,
                                         std::map<std::string,std::string> *replacements=0){
    return RDKit::RWMOL_SPTR(RDKit::SmilesToMol(smi, debugParse, sanitize,replacements));
  }
  static RDKit::RWMOL_SPTR MolFromSmiles(const std::string &smi, const RDKit::SmilesParserParams &params){
    return RDKit::RWMOL_SPTR(RDKit::SmilesToMol(smi, params));
  }
  static RDKit::RWMOL_SPTR MolFromSmarts(const std::string &sma,int debugParse=0,bool mergeHs=false,
                                         std::map<std::string,std::string> *replacements=0){
    return RDKit::RWMOL_SPTR(RDKit::SmartsToMol(sma, debugParse, mergeHs,replacements));
  }
static RDKit::RWMOL_SPTR MolFromMolBlock(const std::string &molB,
                                  bool sanitize=true,bool removeHs=true,bool strictParsing=true){
  RDKit::RWMol *mol=RDKit::MolBlockToMol(molB,sanitize,removeHs,strictParsing);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMolFile(const std::string &filename,
                                 bool sanitize=true,bool removeHs=true,bool strictParsing=true){
  RDKit::RWMol *mol=RDKit::MolFileToMol(filename,sanitize,removeHs,strictParsing);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromTPLFile(const std::string &fName,bool sanitize=true,
                      bool skipFirstConf=false) {
  RDKit::RWMol *mol=0;
    mol=RDKit::TPLFileToMol(fName, sanitize, skipFirstConf);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMol2File(const std::string &fName,bool sanitize=true,bool removeHs=true,
                       RDKit::Mol2Type variant=RDKit::CORINA, bool cleanupSubstructures=true) {
  RDKit::RWMol *mol=0;
  mol=RDKit::Mol2FileToMol(fName, sanitize, removeHs, variant, cleanupSubstructures);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromMol2Block(const std::string &molBlock,bool sanitize=true,bool removeHs=true,
                        RDKit::Mol2Type variant=RDKit::CORINA, bool cleanupSubstructures=true) {
  RDKit::RWMol *mol=0;
    mol=RDKit::Mol2BlockToMol(molBlock, sanitize, removeHs, variant, cleanupSubstructures);
  return RDKit::RWMOL_SPTR(mol);
}

static RDKit::RWMOL_SPTR MolFromPDBBlock(const std::string &molB,
                                         bool sanitize=true,bool removeHs=true,
                                         unsigned int flavor=0,bool proximityBonding=true){
  RDKit::RWMol *mol=0;
  mol=RDKit::PDBBlockToMol(molB,sanitize,removeHs,flavor,proximityBonding);
  return RDKit::RWMOL_SPTR(mol);
}

static RDKit::RWMOL_SPTR MolFromPDBFile(const std::string &fName,
                                        bool sanitize=true,bool removeHs=true,
                                        unsigned int flavor=0,bool proximityBonding=true){
  RDKit::RWMol *mol=0;
  mol=RDKit::PDBFileToMol(fName,sanitize,removeHs,flavor,proximityBonding);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromSequence(const std::string &text,
                                  bool sanitize=true,int flavor=0){
  RDKit::RWMol *mol=0;
  mol=RDKit::SequenceToMol(text,sanitize,flavor);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromFASTA(const std::string &text,
                                  bool sanitize=true,int flavor=0){
  RDKit::RWMol *mol=0;
  mol=RDKit::FASTAToMol(text,sanitize,flavor);
  return RDKit::RWMOL_SPTR(mol);
}
static RDKit::RWMOL_SPTR MolFromHELM(const std::string &text,
                                  bool sanitize=true){
  RDKit::RWMol *mol=0;
  mol=RDKit::HELMToMol(text,sanitize);
  return RDKit::RWMOL_SPTR(mol);
}

static std::vector<RDKit::RWMOL_SPTR> MolsFromCDXML(const std::string &text,
						     bool sanitize=true){
  auto res = RDKit::CDXMLToMols(text, sanitize);
  std::vector<RDKit::RWMOL_SPTR> mols;
  for(auto &mol: res) {
    mols.emplace_back(mol.release());
  }
  return mols;

}

static std::vector<RDKit::RWMOL_SPTR> MolsFromCDXMLFile(const std::string &text,
							 bool sanitize=true){
  auto res = RDKit::CDXMLFileToMols(text, sanitize);
  std::vector<RDKit::RWMOL_SPTR> mols;
  for(auto &mol: res) {
    mols.emplace_back(mol.release());
  }
  return mols;
}


/* Methods from MolFileStereoChem.h */
void DetectAtomStereoChemistry(const RDKit::Conformer *conf) {
  RDKit::DetectAtomStereoChemistry(*($self), conf);
}
void DetectBondStereoChemistry(const RDKit::Conformer *conf) {
  RDKit::DetectBondStereoChemistry(*($self), conf);
}
void ClearSingleBondDirFlags() {
 RDKit::ClearSingleBondDirFlags(*($self));
};
void reapplyMolBlockWedging() {
  RDKit::reapplyMolBlockWedging(*($self));
}
void clearMolBlockWedgingInfo() {
  RDKit::clearMolBlockWedgingInfo(*($self));
}
void invertMolBlockWedgingInfo() {
  RDKit::invertMolBlockWedgingInfo(*($self));
}
void markUnspecifiedStereoAsUnknown(int confId) {
  RDKit::markUnspecifiedStereoAsUnknown(*($self), confId);
}

/* From Kekulize.cpp, MolOps.h */
void Kekulize(bool markAtomsBonds=true, unsigned int maxBackTracks=100) {
  RDKit::MolOps::Kekulize(*($self), markAtomsBonds, maxBackTracks);
}

/* MolOps.h */
void sanitizeMol() {
  RDKit::MolOps::sanitizeMol(*($self));
}

}
