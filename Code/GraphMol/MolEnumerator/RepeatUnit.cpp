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
#include <charconv>

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
namespace MolEnumerator {

const std::string polymarker = "_polymeratom";
const std::string headmarker = "_headatom";
const std::string headmarker_frame = "_headatom_frame";
const std::string tailmarker = "_tailatom";
const std::string tailmarker_frame = "_tailatom_frame";
const std::string headheadmarker = "_headhead";
const unsigned ladderoffset = 100000;

namespace {
void tagAtoms(std::shared_ptr<ROMol> mol, const Bond *bond,
              const boost::dynamic_bitset<> &sgatoms, unsigned int index,
              const std::string &marker, const std::string &framemarker,
              const std::string &connect) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(bond, "bad bond");
  if (sgatoms[bond->getBeginAtomIdx()]) {
    bond->getBeginAtom()->setProp(marker, index);
    if (connect == "HH") {
      bond->getBeginAtom()->setProp(headheadmarker, 1);
    }
    auto frameAtom = mol->getAtomWithIdx(bond->getEndAtomIdx());
    std::vector<unsigned int> vs;
    frameAtom->getPropIfPresent(framemarker, vs);
    vs.push_back(index);
    frameAtom->setProp(framemarker, vs);
  } else if (sgatoms[bond->getEndAtomIdx()]) {
    bond->getEndAtom()->setProp(marker, index);
    if (connect == "HH") {
      bond->getEndAtom()->setProp(headheadmarker, 1);
    }
    auto frameAtom = mol->getAtomWithIdx(bond->getBeginAtomIdx());
    std::vector<unsigned int> vs;
    frameAtom->getPropIfPresent(framemarker, vs);
    vs.push_back(index);
    frameAtom->setProp(framemarker, vs);
  } else {
    throw ValueErrorException("neither atom in an SRU bond is in the polymer");
  }
}
}  // namespace

const size_t RepeatUnitOp::DEFAULT_REPEAT_COUNT = 4;

void RepeatUnitOp::initFromMol(const ROMol &mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}

// if the SRU label has information about the repetition range, we should use
// that information for the enumerator. There are these cases:
//      * if label is a number, N, assume desired range is 0-N
//      * if label is M-N, where M and N are non-negative integers
//        and M <= N, we'll assume the desired range is {M, N}.
//      * well assume the desired range is {0,
//      RepeatUnitOp::DEFAULT_REPEAT_COUNT}
//        for anything else
[[nodiscard]] static std::pair<size_t, size_t> parse_repeat_counts(
    const std::string &sru_label) {
  static constexpr std::pair<size_t, size_t> default_repeat_counts{
      0, RepeatUnitOp::DEFAULT_REPEAT_COUNT};

  if (sru_label.empty()) {
    return default_repeat_counts;
  }

  auto end = sru_label.data() + sru_label.size();

  size_t min_repeats = 0;
  auto status = std::from_chars(sru_label.data(), end, min_repeats);

  // if label is something like 'N', max count should be N
  if (status.ec == std::errc() && status.ptr == end) {
    // offset by one since range is inclusive
    return std::make_pair(0u, min_repeats + 1);
  }

  // if we failed to parse the first number for some reason
  if (status.ec != std::errc() || (status.ptr != end && status.ptr[0] != '-')) {
    return default_repeat_counts;
  }

  size_t max_repeats = 0;
  status = std::from_chars(status.ptr + 1, end, max_repeats);
  if (status.ec != std::errc() || status.ptr != end) {
    return default_repeat_counts;
  }

  // assume it's nonsense
  if (max_repeats < min_repeats) {
    return default_repeat_counts;
  }

  // offset by one since range is inclusive
  return std::make_pair(min_repeats, max_repeats + 1);
}

void RepeatUnitOp::initFromMol() {
  // we're making an assumption here that each atom has at most one bond to an
  // atom not in the repeat unit
  if (!dp_mol) {
    return;
  }

  d_repeats.clear();
  d_countAtEachPoint.clear();

  // start by figuring out which atoms are in sgroups which will be enumerated
  boost::dynamic_bitset<> atomsInSRUs(dp_mol->getNumAtoms());
  std::vector<boost::dynamic_bitset<>> atomsPerSRU;
  std::vector<const SubstanceGroup *> enumerated_SGroups;
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

    std::string label;
    sg.getPropIfPresent("LABEL", label);
    auto [min_repeats, max_repeats] = parse_repeat_counts(label);
    d_minRepeatCounts.push_back(min_repeats);
    d_countAtEachPoint.push_back(
        std::min(max_repeats - min_repeats, d_maxNumRounds));

    // tag the atoms in the repeat unit:
    boost::dynamic_bitset<> sgatoms(dp_mol->getNumAtoms());
    for (auto aidx : sg.getAtoms()) {
      if (atomsInSRUs[aidx]) {
        throw ValueErrorException("cannot enumerate overlapping SRU groups");
      }
      sgatoms.set(aidx);
      atomsInSRUs.set(aidx);
    }
    atomsPerSRU.push_back(sgatoms);

    auto sgIdx = sg.getProp<unsigned>("index");

    // tag the head and tail atoms
    const auto &bnds = sg.getBonds();
    if (bnds.size() == 2) {
      // simple case with only two bonds: here we just have head and tail
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(bnds[0]), sgatoms, sgIdx,
               headmarker, tailmarker_frame, connect);
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(bnds[1]), sgatoms, sgIdx,
               tailmarker, headmarker_frame, connect);
    } else if (bnds.size() == 4) {
      // four bonds are what we see for a ladder polymer, here we need two
      // different heads and two different tails. We mark the second set with
      // a large offset so that they don't get confused, but otherwise the
      // rest of the code handles this automatically.

      // NOTE: theoretically we could support larger numbers of head/tail
      // pairs, but I don't believe these show up in reality

      // We may have XBCORR to indicate which bonds correspond to which
      std::vector<unsigned int> xbcorr;
      sg.getPropIfPresent("XBCORR", xbcorr);
      // if it's not there, or if it wasn't the right size, just use
      // the bonds:
      if (xbcorr.size() != 4) {
        xbcorr = {bnds[0], bnds[2], bnds[1], bnds[3]};
      }
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(xbcorr[0]), sgatoms, sgIdx,
               headmarker, tailmarker_frame, connect);
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(xbcorr[2]), sgatoms,
               sgIdx + ladderoffset, headmarker, tailmarker_frame, connect);
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(xbcorr[1]), sgatoms, sgIdx,
               tailmarker, headmarker_frame, connect);
      tagAtoms(dp_mol, dp_mol->getBondWithIdx(xbcorr[3]), sgatoms,
               sgIdx + ladderoffset, tailmarker, headmarker_frame, connect);
    } else {
      throw ValueErrorException("can only handle SRUs with two or four bonds");
    }

    enumerated_SGroups.push_back(&sg);
  }

  if (!enumerated_SGroups.empty() &&
      dp_mol->hasProp(common_properties::molFileLinkNodes)) {
    throw ValueErrorException(
        "cannot enumerate molecules which include both SRUs and LINKNODEs");
  }

  // copy the molecule over as the frame. We'll remove atoms in SRUs from this
  // below
  dp_frame.reset(new RWMol(*dp_mol));
  dp_frame->beginBatchEdit();

  // now set up the repeat units for each of the SRUs
  for (const auto *sgp : enumerated_SGroups) {
    const auto &sg = *sgp;
    std::shared_ptr<RWMol> repeat(new RWMol(*dp_mol));

    // remove the atoms in the repeat unit from the frame:
    boost::dynamic_bitset<> sgatoms(dp_mol->getNumAtoms());
    for (auto aidx : sg.getAtoms()) {
      sgatoms.set(aidx);
      // remove it from the frame
      dp_frame->removeAtom(aidx);
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
  }
  dp_frame->commitBatchEdit();

}  // namespace MolEnumerator

std::vector<size_t> RepeatUnitOp::getVariationCounts() const {
  return d_countAtEachPoint;
}

namespace {
void flipHeadHeadGroups(unsigned int origAtomCount, RWMol &mol) {
  for (auto aidx = origAtomCount; aidx < mol.getNumAtoms(); ++aidx) {
    auto atom = mol.getAtomWithIdx(aidx);
    if (atom->hasProp(headheadmarker)) {
      if (atom->hasProp(headmarker)) {
        atom->setProp(tailmarker, atom->getProp<unsigned>(headmarker));
        atom->clearProp(headmarker);
      } else if (atom->hasProp(tailmarker)) {
        atom->setProp(headmarker, atom->getProp<unsigned>(tailmarker));
        atom->clearProp(tailmarker);
      }
      if (atom->hasProp(headmarker_frame)) {
        atom->setProp(
            tailmarker_frame,
            atom->getProp<std::vector<unsigned int>>(headmarker_frame));
        atom->clearProp(headmarker_frame);
      } else if (atom->hasProp(tailmarker_frame)) {
        atom->setProp(
            headmarker_frame,
            atom->getProp<std::vector<unsigned int>>(tailmarker_frame));
        atom->clearProp(tailmarker_frame);
      }
    }
  }
}

void connectRepeatAtomsAndRemoveExtras(unsigned int origAtomCount, RWMol &mol) {
  mol.beginBatchEdit();
  for (auto aidx1 = 0u; aidx1 < origAtomCount; ++aidx1) {
    auto at1 = mol.getAtomWithIdx(aidx1);
    unsigned tailIdx;
    if (at1->getPropIfPresent(tailmarker, tailIdx)) {
      bool connected = false;
      for (auto aidx2 = origAtomCount; aidx2 < mol.getNumAtoms(); ++aidx2) {
        auto at2 = mol.getAtomWithIdx(aidx2);
        unsigned int headIdx;
        if (at2->getPropIfPresent(headmarker, headIdx) && tailIdx == headIdx) {
          connected = true;
          at2->clearProp(headmarker);
          // remove any atom connected to the head which isn't in the
          // repeat unit
          for (const auto &nbri :
               boost::make_iterator_range(mol.getAtomNeighbors(at2))) {
            auto nbr = mol[nbri];
            if (!nbr->hasProp(polymarker)) {
              mol.removeAtom(nbr);
            }
          }
          // FIX: not dealing with multiple bonds from the repeat unit
          mol.addBond(at1, at2, Bond::BondType::SINGLE);
          break;
        }
      }
      if (connected) {
        at1->clearProp(tailmarker);
        // remove any atom connected to the tail which isn't in the
        // repeat unit
        for (const auto &nbri :
             boost::make_iterator_range(mol.getAtomNeighbors(at1))) {
          auto nbr = mol[nbri];
          if (!nbr->hasProp(polymarker)) {
            mol.removeAtom(nbr);
          }
        }
      } else {
        throw ValueErrorException("no head found for tail");
      }
    }
  }
  mol.commitBatchEdit();
}

void connectRepeatToFrame(unsigned int nOrigAtoms, RWMol &mol,
                          std::map<unsigned, Atom *> &headMap,
                          std::map<unsigned, Atom *> &tailMap) {
  for (auto aidx = nOrigAtoms; aidx < mol.getNumAtoms(); ++aidx) {
    auto sruAtom = mol.getAtomWithIdx(aidx);
    unsigned int val;
    if (sruAtom->getPropIfPresent(headmarker, val)) {
      if (tailMap.find(val) != tailMap.end()) {
        // there's an atom in the frame to connect to:
        mol.addBond(sruAtom, tailMap[val], Bond::BondType::SINGLE);
        sruAtom->clearProp(headmarker);
        tailMap.erase(val);
      }
    }
    if (sruAtom->getPropIfPresent(tailmarker, val)) {
      if (headMap.find(val) != headMap.end()) {
        // there's an atom in the frame to connect to:
        mol.addBond(sruAtom, headMap[val], Bond::BondType::SINGLE);
        sruAtom->clearProp(tailmarker);
        headMap.erase(val);
      }
    }
    std::vector<unsigned int> vals;
    if (sruAtom->getPropIfPresent(headmarker_frame, vals)) {
      for (const auto val : vals) {
        headMap[val] = sruAtom;
      }
      sruAtom->clearProp(headmarker_frame);
    }
    if (sruAtom->getPropIfPresent(tailmarker_frame, vals)) {
      for (const auto val : vals) {
        tailMap[val] = sruAtom;
      }
      sruAtom->clearProp(tailmarker_frame);
    }
  }
}

void constructHeadAndTailMaps(RWMol &mol, std::map<unsigned, Atom *> &headMap,
                              std::map<unsigned, Atom *> &tailMap) {
  for (auto atom : mol.atoms()) {
    std::vector<unsigned int> vals;
    if (atom->getPropIfPresent(headmarker_frame, vals)) {
      for (auto val : vals) {
        if (headMap.find(val) != headMap.end()) {
          throw ValueErrorException(
              "SRU group present with multiple head atoms with the same index");
        }
        headMap[val] = atom;
      }
    }
    if (atom->getPropIfPresent(tailmarker_frame, vals)) {
      for (auto val : vals) {
        if (tailMap.find(val) != tailMap.end()) {
          throw ValueErrorException(
              "SRU group present with multiple tail atoms with the same index");
        }
        tailMap[val] = atom;
      }
    }
  }
}
}  // namespace

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
  // we will use these maps from head/tail markerss to atoms in the frame
  // later
  std::map<unsigned, Atom *> headMap;
  std::map<unsigned, Atom *> tailMap;
  constructHeadAndTailMaps(*res, headMap, tailMap);

  for (size_t i = 0; i < which.size(); ++i) {
    RWMol filling;
    // if there SRU sgroups has a repetition range like 3-5, we should make
    // sure each output has at least 3 repetitions
    auto offset = d_minRepeatCounts[i];
    for (size_t iter = 0; iter < offset + which[i]; ++iter) {
      auto origAtomCount = filling.getNumAtoms();
      filling.insertMol(*d_repeats[i]);

      if (iter % 2) {
        // check for any HH groups on odd iterations and invert them
        flipHeadHeadGroups(origAtomCount, filling);
      }
      // if we aren't adding the first repeat unit, then we need to connect
      // things
      if (iter) {
        connectRepeatAtomsAndRemoveExtras(origAtomCount, filling);
      }
    }  // end of loop over iteration

    // ok, add the fragment generated by enumerating that SRU to the result:
    unsigned int nOrigAtoms = res->getNumAtoms();
    res->insertMol(filling);

    // and connect it to the frame:
    connectRepeatToFrame(nOrigAtoms, *res, headMap, tailMap);

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
