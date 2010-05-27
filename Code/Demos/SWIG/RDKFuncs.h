// $Id$
//
// Copyright (C) 2008-2010 Greg Landrum
// All Rights Reserved
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
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
  RDKit::ChemicalReaction *res=RDKit::RxnSmartsToChemicalReaction(sma);
  if(res) res->initReactantMatchers();
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnBlock(std::string block){
  RDKit::ChemicalReaction *res=RDKit::RxnBlockToChemicalReaction(block);
  if(res) res->initReactantMatchers();
  return res;
};
RDKit::ChemicalReaction *ReactionFromRxnFile(std::string filename){
  RDKit::ChemicalReaction *res=RDKit::RxnFileToChemicalReaction(filename);
  if(res) res->initReactantMatchers();
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

std::string MolToBinary(RDKit::ROMOL_SPTR mol){
  std::string res;
  RDKit::MolPickler::pickleMol(*mol,res);
  return res;
};


std::string rdkitVersion(){
  return RDKit::rdkitVersion;
}

