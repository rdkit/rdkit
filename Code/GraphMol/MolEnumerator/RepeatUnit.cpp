//
//  Copyright (C) 2021 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>
#include <algorithm>

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
namespace MolEnumerator {

const std::string polymarker = "_polymeratom";
const std::string headmarker = "_headatom";
const std::string headmarker_frame = "_headatom_frame";
const std::string tailmarker = "_tailatom";
const std::string tailmarker_frame = "_tailatom_frame";
const std::string headheadmarker = "_headhead";

void RepeatUnitOp::initFromMol(const ROMol &mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}
void RepeatUnitOp::initFromMol() {
  // we're making an assumption here that each atom has at most one bond to an
  // atom not in the repeat unit
  if (!dp_mol) {
    return;
  }

  d_repeats.clear();
  d_countAtEachPoint.clear();
  dp_frame.reset(new RWMol(*dp_mol));
  dp_frame->beginBatchEdit();
  boost::dynamic_bitset<> atomsSeen(dp_mol->getNumAtoms());
  for (auto &sg : getSubstanceGroups(*dp_mol)) {
    std::string typ;
    if (!sg.getPropIfPresent("TYPE", typ) || typ != "SRU") {
      continue;
    }
    std::string connect;
    sg.getPropIfPresent("CONNECT", connect);
    if (!connect.empty()) {
      if (connect != "HT" && connect != "HH") {
        BOOST_LOG(rdWarningLog)
            << "can only enumerate SRUs with CONNECT=HT or CONNECT=HH"
            << std::endl;
        continue;
      }
    }
    // FIX: need to come back and fix this for ladder polymers
    if (sg.getBonds().size() != 2) {
      BOOST_LOG(rdWarningLog)
          << "can only handle SRUs with two bonds" << std::endl;
      continue;
    }
    std::shared_ptr<RWMol> repeat(new RWMol(*dp_mol));

    // tag the atoms in the repeat unit:
    boost::dynamic_bitset<> sgatoms(dp_mol->getNumAtoms());
    for (auto aidx : sg.getAtoms()) {
      if (atomsSeen[aidx]) {
        throw ValueErrorException("cannot enumerate overlapping SRU groups");
      }
      dp_mol->getAtomWithIdx(aidx)->setProp(polymarker, 1);
      sgatoms.set(aidx);
      // remove it from the frame
      dp_frame->removeAtom(aidx);
    }
    atomsSeen |= sgatoms;

    // tag the head and tail atoms
    const auto headBond = repeat->getBondWithIdx(sg.getBonds()[0]);
    if (sgatoms[headBond->getBeginAtomIdx()]) {
      headBond->getBeginAtom()->setProp(headmarker,
                                        sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        headBond->getBeginAtom()->setProp(headheadmarker, 1);
      }
      auto frameAtom = dp_frame->getAtomWithIdx(headBond->getEndAtomIdx());
      std::vector<unsigned int> vs;
      frameAtom->getPropIfPresent(tailmarker_frame, vs);
      vs.push_back(sg.getProp<unsigned>("index"));
      frameAtom->setProp(tailmarker_frame, vs);
    } else if (sgatoms[headBond->getEndAtomIdx()]) {
      headBond->getEndAtom()->setProp(headmarker,
                                      sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        headBond->getEndAtom()->setProp(headheadmarker, 1);
      }
      auto frameAtom = dp_frame->getAtomWithIdx(headBond->getBeginAtomIdx());
      std::vector<unsigned int> vs;
      frameAtom->getPropIfPresent(tailmarker_frame, vs);
      vs.push_back(sg.getProp<unsigned>("index"));
      frameAtom->setProp(tailmarker_frame, vs);
    } else {
      throw ValueErrorException(
          "neither atom in an SRU bond is in the polymer");
    }
    const auto tailBond = repeat->getBondWithIdx(sg.getBonds()[1]);
    if (sgatoms[tailBond->getBeginAtomIdx()]) {
      tailBond->getBeginAtom()->setProp(tailmarker,
                                        sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        tailBond->getBeginAtom()->setProp(headheadmarker, 1);
      }
      auto frameAtom = dp_frame->getAtomWithIdx(tailBond->getEndAtomIdx());
      std::vector<unsigned int> vs;
      frameAtom->getPropIfPresent(headmarker_frame, vs);
      vs.push_back(sg.getProp<unsigned>("index"));
      frameAtom->setProp(headmarker_frame, vs);
    } else if (sgatoms[tailBond->getEndAtomIdx()]) {
      tailBond->getEndAtom()->setProp(tailmarker,
                                      sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        tailBond->getEndAtom()->setProp(headheadmarker, 1);
      }
      auto frameAtom = dp_frame->getAtomWithIdx(tailBond->getBeginAtomIdx());
      std::vector<unsigned int> vs;
      frameAtom->getPropIfPresent(headmarker_frame, vs);
      vs.push_back(sg.getProp<unsigned>("index"));
      frameAtom->setProp(headmarker_frame, vs);
    } else {
      throw ValueErrorException(
          "neither atom in an SRU bond is in the polymer");
    }
    repeat->beginBatchEdit();
    for (auto aidx = 0u; aidx < repeat->getNumAtoms(); ++aidx) {
      if (!sgatoms[aidx]) {
        repeat->removeAtom(aidx);
      } else {
        repeat->getAtomWithIdx(aidx)->setProp(polymarker, 1);
      }
    }
    repeat->commitBatchEdit();
    d_repeats.push_back(repeat);
    d_countAtEachPoint.push_back(d_defaultRepeatCount);
  }
  dp_frame->commitBatchEdit();

}  // namespace MolEnumerator

std::vector<size_t> RepeatUnitOp::getVariationCounts() const {
  return d_countAtEachPoint;
}

std::unique_ptr<ROMol> RepeatUnitOp::operator()(
    const std::vector<size_t> &which) const {
  PRECONDITION(dp_mol, "no molecule");
  PRECONDITION(dp_frame, "not initialized");
  if (which.size() != d_countAtEachPoint.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }
  // quick error checking before we do any work:
  for (size_t i = 0; i < which.size(); ++i) {
    if (which[i] >= d_countAtEachPoint[i]) {
      throw ValueErrorException("bad element value in enumeration");
    }
  }

  std::unique_ptr<RWMol> res{new RWMol()};
  res->insertMol(*dp_frame);

  // ---------------------
  // we will use these maps of atoms in the frame with head/tail markers later
  std::map<unsigned, Atom *> headMap;
  std::map<unsigned, Atom *> tailMap;
  for (auto atom : res->atoms()) {
    std::vector<unsigned int> vals;
    if (atom->getPropIfPresent(headmarker_frame, vals)) {
      for (auto val : vals) {
        if (headMap.find(val) != headMap.end()) {
          throw ValueErrorException(
              "SRU group present with more than one head atom");
        }
        headMap[val] = atom;
      }
    }
    if (atom->getPropIfPresent(tailmarker_frame, vals)) {
      for (auto val : vals) {
        if (tailMap.find(val) != tailMap.end()) {
          throw ValueErrorException(
              "SRU group present with more than one tail atom");
        }
        tailMap[val] = atom;
      }
    }
  }

  for (size_t i = 0; i < which.size(); ++i) {
    RWMol filling;
    for (size_t iter = 0; iter < which[i]; ++iter) {
      auto origAtomCount = filling.getNumAtoms();
      filling.insertMol(*d_repeats[i]);

      if (iter % 2) {
        // check for any HH groups on odd iterations and invert them
        for (auto aidx2 = origAtomCount; aidx2 < filling.getNumAtoms();
             ++aidx2) {
          auto at2 = filling.getAtomWithIdx(aidx2);
          if (at2->hasProp(headheadmarker)) {
            if (at2->hasProp(headmarker)) {
              at2->setProp(tailmarker, at2->getProp<unsigned>(headmarker));
              at2->clearProp(headmarker);
            } else if (at2->hasProp(tailmarker)) {
              at2->setProp(headmarker, at2->getProp<unsigned>(tailmarker));
              at2->clearProp(tailmarker);
            }
          }
        }
      }
      filling.beginBatchEdit();
      // if we aren't adding the first repeat unit, then we need to connect
      // things
      if (iter) {
        // std::cerr << "ITER " << iter << std::endl;
        for (auto aidx1 = 0u; aidx1 < origAtomCount; ++aidx1) {
          auto at1 = filling.getAtomWithIdx(aidx1);
          unsigned tailIdx;
          if (at1->getPropIfPresent(tailmarker, tailIdx)) {
            // std::cerr << "   tail: " << aidx1 << " " << tailIdx << std::endl;
            bool connected = false;
            for (auto aidx2 = origAtomCount; aidx2 < filling.getNumAtoms();
                 ++aidx2) {
              auto at2 = filling.getAtomWithIdx(aidx2);
              unsigned int headIdx;
              if (at2->getPropIfPresent(headmarker, headIdx) &&
                  tailIdx == headIdx) {
                // std::cerr << "   head: " << aidx2 << " " << headIdx
                //           << std::endl;
                connected = true;
                at2->clearProp(headmarker);
                // remove any atom connected to the head which isn't in the
                // repeat unit
                for (const auto &nbri : boost::make_iterator_range(
                         filling.getAtomNeighbors(at2))) {
                  auto nbr = filling[nbri];
                  if (!nbr->hasProp(polymarker)) {
                    filling.removeAtom(nbr);
                  }
                }
                // FIX: not dealing with multiple bonds from the repeat unit
                filling.addBond(at1, at2, Bond::BondType::SINGLE);
                break;
              }
            }
            if (connected) {
              at1->clearProp(tailmarker);
              // remove any atom connected to the tail which isn't in the
              // repeat unit
              for (const auto &nbri :
                   boost::make_iterator_range(filling.getAtomNeighbors(at1))) {
                auto nbr = filling[nbri];
                if (!nbr->hasProp(polymarker)) {
                  filling.removeAtom(nbr);
                }
              }
            } else {
              throw ValueErrorException("no head found for tail");
            }
          }
        }
      }
      filling.commitBatchEdit();
    }  // end of loop over iteration

    // ok, add the fragment generated by enumerating that SRU to the result:
    unsigned int nOrigAtoms = res->getNumAtoms();
    res->insertMol(filling);

    // and connect it to the frame:
    for (auto aidx = nOrigAtoms; aidx < res->getNumAtoms(); ++aidx) {
      auto sruAtom = res->getAtomWithIdx(aidx);
      unsigned int val;
      if (sruAtom->getPropIfPresent(headmarker, val)) {
        if (tailMap.find(val) != tailMap.end()) {
          // there's an atom in the frame to connect to:
          res->addBond(sruAtom, tailMap[val], Bond::BondType::SINGLE);
          sruAtom->clearProp(headmarker);
          tailMap.erase(val);
        }
      }
      if (sruAtom->getPropIfPresent(tailmarker, val)) {
        if (headMap.find(val) != headMap.end()) {
          // there's an atom in the frame to connect to:
          res->addBond(sruAtom, headMap[val], Bond::BondType::SINGLE);
          sruAtom->clearProp(tailmarker);
          headMap.erase(val);
        }
      }
    }

  }  // end of loop over SRU
  // connect any remaining dangling heads and tails
  for (auto &tpl : headMap) {
    auto iter = tailMap.find(tpl.first);
    if (iter != tailMap.end()) {
      res->addBond(tpl.second, iter->second, Bond::BondType::SINGLE);
      tpl.second->clearProp(headmarker_frame);
      iter->second->clearProp(tailmarker_frame);
    }
  }
  res->commitBatchEdit();

  return std::unique_ptr<ROMol>(new ROMol(*res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
