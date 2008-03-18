// $Id$
//
// Copyright (C) 2008 Greg Landrum
// All Rights Reserved
//

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include "RDKFuncs.h"

using namespace RDKit;
//ROMol *MolFromSmiles(std::string smi) {
ROMOL_SPTR MolFromSmiles(std::string smi) {
  return ROMOL_SPTR(SmilesToMol(smi));;
};
ROMOL_SPTR MolFromSmarts(std::string sma) {
  return ROMOL_SPTR(SmartsToMol(sma));
};
std::string MolToSmiles(ROMOL_SPTR mol,bool doIsomericSmiles,
                        bool doKekule, int rootedAtAtom) {
  return MolToSmiles(*mol,doIsomericSmiles,doKekule,rootedAtAtom);
};

ChemicalReaction *ReactionFromSmarts(std::string sma) {
  return RxnSmartsToChemicalReaction(sma);
};




