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
const std::string tailmarker = "_tailatom";
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

  dp_frame.reset(new RWMol(*dp_mol));
  bool canEnumerate = false;
  for (auto &sg : getSubstanceGroups(*dp_frame)) {
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
    // tag the atoms in the repeat unit:
    boost::dynamic_bitset<> sgatoms(dp_frame->getNumAtoms());
    for (auto aidx : sg.getAtoms()) {
      dp_frame->getAtomWithIdx(aidx)->setProp(polymarker, 1);
      sgatoms.set(aidx);
    }
    // tag the head and tail atoms
    // FIX: need to come back and fix this for ladder polymers
    if (sg.getBonds().size() != 2) {
      BOOST_LOG(rdWarningLog)
          << "can only handle SRUs with two bonds" << std::endl;
      continue;
    }

    const auto headBond = dp_frame->getBondWithIdx(sg.getBonds()[0]);
    if (sgatoms[headBond->getBeginAtomIdx()]) {
      headBond->getBeginAtom()->setProp(headmarker,
                                        sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        headBond->getBeginAtom()->setProp(headheadmarker, 1);
      }
    } else if (sgatoms[headBond->getEndAtomIdx()]) {
      headBond->getEndAtom()->setProp(headmarker,
                                      sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        headBond->getEndAtom()->setProp(headheadmarker, 1);
      }
    } else {
      throw ValueErrorException(
          "neither atom in an SRU bond is in the polymer");
    }
    const auto tailBond = dp_frame->getBondWithIdx(sg.getBonds()[1]);
    if (sgatoms[tailBond->getBeginAtomIdx()]) {
      tailBond->getBeginAtom()->setProp(tailmarker,
                                        sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        tailBond->getBeginAtom()->setProp(headheadmarker, 1);
      }
    } else if (sgatoms[tailBond->getEndAtomIdx()]) {
      tailBond->getEndAtom()->setProp(tailmarker,
                                      sg.getProp<unsigned>("index"));
      if (connect == "HH") {
        tailBond->getEndAtom()->setProp(headheadmarker, 1);
      }
    } else {
      throw ValueErrorException(
          "neither atom in an SRU bond is in the polymer");
    }
    canEnumerate = true;
  }
  d_countAtEachPoint.clear();
  if (canEnumerate) {
    d_countAtEachPoint.push_back(d_defaultRepeatCount);
  }
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

  //   for (const auto at : dp_frame->atoms()) {
  //     unsigned headidx, tailidx;
  //     if (at->getPropIfPresent(headmarker, headidx)) {
  //       std::cerr << " HEAD " << at->getIdx() << " : " << headidx <<
  //       std::endl;
  //     }
  //     if (at->getPropIfPresent(tailmarker, tailidx)) {
  //       std::cerr << " TAIL " << at->getIdx() << " : " << tailidx <<
  //       std::endl;
  //     }
  //   }

  std::unique_ptr<RWMol> res{new RWMol()};
  for (size_t i = 0; i < which.size(); ++i) {
    for (size_t iter = 0; iter < which[i]; ++iter) {
      auto origAtomCount = res->getNumAtoms();
      res->insertMol(*dp_frame);

      if (iter % 2) {
        // check for any HH groups on odd iterations and invert them
        for (auto aidx2 = origAtomCount; aidx2 < res->getNumAtoms(); ++aidx2) {
          auto at2 = res->getAtomWithIdx(aidx2);
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
      res->beginBatchEdit();
      // if we aren't adding the first repeat unit, then we need to connect
      // things
      if (iter) {
        // std::cerr << "ITER " << iter << std::endl;
        for (auto aidx1 = 0u; aidx1 < origAtomCount; ++aidx1) {
          auto at1 = res->getAtomWithIdx(aidx1);
          unsigned tailIdx;
          if (at1->getPropIfPresent(tailmarker, tailIdx)) {
            // std::cerr << "   tail: " << aidx1 << " " << tailIdx << std::endl;
            bool connected = false;
            for (auto aidx2 = origAtomCount; aidx2 < res->getNumAtoms();
                 ++aidx2) {
              auto at2 = res->getAtomWithIdx(aidx2);
              unsigned int headIdx;
              if (at2->getPropIfPresent(headmarker, headIdx) &&
                  tailIdx == headIdx) {
                // std::cerr << "   head: " << aidx2 << " " << headIdx
                //           << std::endl;
                connected = true;
                at2->clearProp(headmarker);
                // remove any atom connected to the head which isn't in the
                // repeat unit
                for (const auto &nbri :
                     boost::make_iterator_range(res->getAtomNeighbors(at2))) {
                  auto nbr = (*res)[nbri];
                  if (!nbr->hasProp(polymarker)) {
                    res->removeAtom(nbr);
                  }
                }
                // FIX: not dealing with multiple bonds from the repeat unit
                res->addBond(at1, at2, Bond::BondType::SINGLE);
                break;
              }
            }
            if (connected) {
              at1->clearProp(tailmarker);
              // remove any atom connected to the tail which isn't in the
              // repeat unit
              for (const auto &nbri :
                   boost::make_iterator_range(res->getAtomNeighbors(at1))) {
                auto nbr = (*res)[nbri];
                if (!nbr->hasProp(polymarker)) {
                  res->removeAtom(nbr);
                }
              }
            } else {
              throw ValueErrorException("no head found for tail");
            }
          }
        }
      }
      res->commitBatchEdit();
    }
  }

  return std::unique_ptr<ROMol>(new ROMol(*res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
