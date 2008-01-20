// $Id$
//
// Copyright (C) 2008 Greg Landrum
// All Rights Reserved
//

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "RDKFuncs.h"

using namespace RDKit;
ROMol *MolFromSmiles(std::string smi) {
  return SmilesToMol(smi);
};
ROMol *MolFromSmarts(std::string sma) {
  return SmartsToMol(sma);
};
std::string MolToSmiles(ROMol *mol,bool doIsomericSmiles,
                        bool doKekule, int rootedAtAtom) {
  return MolToSmiles(*mol,doIsomericSmiles,doKekule,rootedAtAtom);
};
bool HasSubstructMatch(ROMol &mol,ROMol &query,bool useChirality,
                       bool registerQuery){
  MatchVectType mv;
  return SubstructMatch(mol,query,mv,true,useChirality,registerQuery);
};





