//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../Substruct/SubstructMatch.h"
#include "StructChecker.h"
#include "Pattern.h"
#include "Tautomer.h"

namespace RDKit {
namespace StructureCheck {

bool StructCheckTautomer::applyTautomer(unsigned it) {
  if (Options.FromTautomer.size() <= it || Options.ToTautomer.size() <= it) {
    if (Options.Verbose)
      BOOST_LOG(rdInfoLog) << "ERROR: incorrect Tautomer index it=" << it
                           << "\n";
    return false;
  }
  const ROMol &fromTautomer = *Options.FromTautomer[it];
  const ROMol &toTautomer = *Options.ToTautomer[it];
  if (toTautomer.getNumAtoms() != fromTautomer.getNumAtoms()) {
    if (Options.Verbose)
      BOOST_LOG(rdInfoLog) << "ERROR: incorrect data toTautomer.getNumAtoms() "
                              "!= fromTautomer.getNumAtoms()\n";
    // incorrect data
    // throw(.....);
    return false;
  }
  const unsigned nta = toTautomer.getNumAtoms();
  const unsigned ntb = toTautomer.getNumBonds();
  MatchVectType match;  // The format is (queryAtomIdx, molAtomIdx)
  std::vector<unsigned> atomIdxMap(
      Mol.getNumAtoms());  // matched tau atom indeces

  if (!SubstructMatch(Mol, *Options.FromTautomer[it],
                      match))  // SSMatch(mp, from_tautomer, SINGLE_MATCH);
    return false;
  if (Options.Verbose)
    BOOST_LOG(rdInfoLog) << "found match for from_tautomer with " << nta
                         << " atoms\n";
  // init
  for (unsigned i = 0; i < Mol.getNumAtoms(); i++) atomIdxMap[i] = -1;
  for (MatchVectType::const_iterator mit = match.begin(); mit != match.end();
       mit++) {
    unsigned tai = mit->first;   // From and To Tautomer Atom index
    unsigned mai = mit->second;  // Mol Atom index
    atomIdxMap[mai] = tai;
  }
  // scan for completely mapped bonds and replace bond order with mapped bond
  // from to_tautomer
  for (RDKit::BondIterator_ bond = Mol.beginBonds(); bond != Mol.endBonds();
       bond++) {
    unsigned ti = atomIdxMap[(*bond)->getBeginAtomIdx()];
    unsigned tj = atomIdxMap[(*bond)->getEndAtomIdx()];
    if (-1 == ti || -1 == tj) continue;
    const Bond *tb = toTautomer.getBondBetweenAtoms(ti, tj);
    if (tb && (*bond)->getBondType() != tb->getBondType()) {
      (*bond)->setBondType(tb->getBondType());
    }
  }
  // apply charge/radical fixes if any
  for (unsigned i = 0; i < match.size(); i++) {
    Atom &atom = *Mol.getAtomWithIdx(match[i].second);
    const Atom &ta = *toTautomer.getAtomWithIdx(match[i].first);
    atom.setFormalCharge(ta.getFormalCharge());
    atom.setNumRadicalElectrons(ta.getNumRadicalElectrons());
  }

  return true;
}

}  // namespace StructureCheck
}  // namespace RDKit
