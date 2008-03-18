//
// Copyright (C) 2008 Greg Landrum
// All Rights Reserved
//
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>


//RDKit::ROMol *MolFromSmiles(std::string smi);
RDKit::ROMOL_SPTR MolFromSmiles(std::string smi);
RDKit::ROMOL_SPTR MolFromSmarts(std::string sma);
std::string MolToSmiles(RDKit::ROMOL_SPTR mol,bool doIsomericSmiles=false,
                        bool doKekule=false, int rootedAtAtom=-1);
RDKit::ChemicalReaction *ReactionFromSmarts(std::string sma);
bool hasSubstructMatch(RDKit::ROMol &mol,RDKit::ROMol &query,bool useChirality=false,
                       bool registerQuery=false);
