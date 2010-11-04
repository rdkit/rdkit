// $Id$
//
// Copyright (C) 2008-2010 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <RDGeneral/versions.h>

RDKit::ROMOL_SPTR MolFromSmiles(std::string smi){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::SmilesToMol(smi));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromSmarts(std::string sma){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::SmartsToMol(sma));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromMolBlock(std::string molB,
                                  bool sanitize=true,bool removeHs=true){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::MolBlockToMol(molB,sanitize,removeHs));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ROMOL_SPTR MolFromMolFile(std::string filename,
                                 bool sanitize=true,bool removeHs=true){
  RDKit::ROMol *mol=0;
  try{
    mol=static_cast<RDKit::ROMol *>(RDKit::MolFileToMol(filename,sanitize,removeHs));
  } catch (...){
    mol=0;
  }
  return RDKit::ROMOL_SPTR(mol);
};
RDKit::ChemicalReaction *ReactionFromSmarts(std::string sma){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnSmartsToChemicalReaction(sma);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnBlock(std::string block){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnBlockToChemicalReaction(block);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnFile(std::string filename){
  RDKit::ChemicalReaction *res=0;
  try {
    res=RDKit::RxnFileToChemicalReaction(filename);
    if(res) res->initReactantMatchers();
  } catch (...){
    res=0;
  }
  return res;
};

std::string MolToSmiles(RDKit::ROMOL_SPTR mol,bool doIsomericSmiles=false,
                        bool doKekule=false, int rootedAtAtom=-1){
  return RDKit::MolToSmiles(*mol,doIsomericSmiles,doKekule,rootedAtAtom);
};
std::string MolToMolBlock(RDKit::ROMOL_SPTR mol, bool includeStereo=true, 
                          int confId=-1) {
  return RDKit::MolToMolBlock(*mol,includeStereo,confId);
}

std::vector<int> MolToBinary(RDKit::ROMOL_SPTR mol){
  std::string sres;
  RDKit::MolPickler::pickleMol(*mol,sres);
  std::vector<int> res(sres.length());
  std::copy(sres.begin(),sres.end(),res.begin());
  return res;
};
RDKit::ROMOL_SPTR MolFromBinary(std::vector<int> pkl){
  std::string sres;
  sres.resize(pkl.size());
  std::copy(pkl.begin(),pkl.end(),sres.begin());
  RDKit::ROMol *res=new RDKit::ROMol(sres);
  return RDKit::ROMOL_SPTR(res);
};


std::vector<char> RxnToBinary(RDKit::ChemicalReaction *rxn){
  std::string sres;
  RDKit::ReactionPickler::pickleReaction(rxn,sres);
  std::vector<char> res(sres.length());
  std::copy(sres.begin(),sres.end(),res.begin());
  return res;
};
RDKit::ChemicalReaction *RxnFromBinary(std::vector<char> pkl){
  std::string sres;
  sres.resize(pkl.size());
  std::copy(pkl.begin(),pkl.end(),sres.begin());
  RDKit::ChemicalReaction *res=new RDKit::ChemicalReaction(sres);
  return res;
};

std::string ReactionToSmarts(RDKit::ChemicalReaction *rxn){
  return RDKit::ChemicalReactionToRxnSmarts(*rxn);
};

std::string rdkitVersion(){
  return RDKit::rdkitVersion;
}

unsigned int compute2DCoords(RDKit::ROMol &mol,bool canonOrient=false,
                             bool clearConfs=true){
  return RDDepict::compute2DCoords(mol,0,canonOrient,clearConfs);
}
