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
#if 0
MatchVectType GetSubstructMatches(ROMol &mol,ROMol &query,bool useChirality,
                                  bool registerQuery){
  MatchVectType mv;
  return SubstructMatch(mol,query,mv,true,useChirality,registerQuery);
};
#endif
std::string PickleMol(RDKit::ROMol *mol){
  std::string res="";
  RDKit::MolPickler::pickleMol(mol,res);
  return res;
}
RDKit::ROMol *MolFromPickle(std::string pkl){
  RDKit::ROMol *res=new RDKit::ROMol(pkl);
  return res;
}




