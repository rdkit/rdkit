// $Id$
//
// Copyright (C) 2008-2009 Greg Landrum
// All Rights Reserved
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>

RDKit::ROMOL_SPTR MolFromSmiles(std::string smi){
  return RDKit::ROMOL_SPTR(RDKit::SmilesToMol(smi));
};
RDKit::ROMOL_SPTR MolFromSmarts(std::string sma){
  return RDKit::ROMOL_SPTR(RDKit::SmartsToMol(sma));
};
RDKit::ROMOL_SPTR MolFromMolBlock(std::string molB,
                                  bool sanitize=true,bool removeHs=true){
  return RDKit::ROMOL_SPTR(RDKit::MolBlockToMol(molB,sanitize,removeHs));
};
RDKit::ROMOL_SPTR MolFromMolFile(std::string filename,
                                 bool sanitize=true,bool removeHs=true){
  return RDKit::ROMOL_SPTR(RDKit::MolFileToMol(filename,sanitize,removeHs));
};
RDKit::ChemicalReaction *ReactionFromSmarts(std::string sma){
  return RDKit::RxnSmartsToChemicalReaction(sma);
};
RDKit::ChemicalReaction *ReactionFromRxnBlock(std::string block){
  return RDKit::RxnBlockToChemicalReaction(block);
};
RDKit::ChemicalReaction *ReactionFromRxnFile(std::string filename){
  return RDKit::RxnFileToChemicalReaction(filename);
};

std::string MolToSmiles(RDKit::ROMOL_SPTR mol,bool doIsomericSmiles=false,
                        bool doKekule=false, int rootedAtAtom=-1){
  return RDKit::MolToSmiles(*mol,doIsomericSmiles,doKekule,rootedAtAtom);
};
std::string MolToMolBlock(RDKit::ROMOL_SPTR mol, bool includeStereo=true, 
                          int confId=-1) {
  return RDKit::MolToMolBlock(*mol,includeStereo,confId);
}



