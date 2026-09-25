//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

#include "../MolOps.h"
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../Substruct/SubstructMatch.h"
#include <GraphMol/new_canon.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include "MMPA.h"

// #define MMPA_DEBUG // enable debug info output

namespace RDKit {
namespace MMPA {

typedef std::vector<std::pair<unsigned, unsigned>>
    BondVector_t;  // pair of BeginAtomIdx, EndAtomIdx

namespace detail {
unsigned long long computeMorganCodeHash(const ROMol &mol) {
  auto nv = mol.getNumAtoms();
  std::vector<unsigned long> currCodes;
  currCodes.reserve(nv);
  std::vector<unsigned long> prevCodes;
  auto nIterations = mol.getNumBonds();
  if (nIterations > 5) {
    nIterations = 5;
  }

  for (const auto a : mol.atoms()) {
    auto atomCode = a->getAtomicNum();
    atomCode |= a->getIsotope() << 8;

    auto charge = a->getFormalCharge();
    atomCode |= std::abs(charge) << 16;
    if (charge < 0) {
      atomCode |= 1 << 29;
    }
    atomCode |= (a->getIsAromatic() ? 1 : 0) << 30;
    currCodes.push_back(atomCode);
  }

  for (size_t iter = 0; iter < nIterations; iter++) {
    prevCodes = currCodes;

    for (const auto bond : mol.bonds()) {
      unsigned order = bond->getBondType();
      unsigned atom1 = bond->getBeginAtomIdx();
      unsigned atom2 = bond->getEndAtomIdx();
      unsigned v1 = prevCodes[atom1];
      unsigned v2 = prevCodes[atom2];

      currCodes[atom1] += v2 * v2 + (v2 + 23) * (order + 1721);
      currCodes[atom2] += v1 * v1 + (v1 + 23) * (order + 1721);
    }
  }

  return std::accumulate(currCodes.begin(), currCodes.end(), 0ULL,
                         [](unsigned long long acc, unsigned long code) {
                           return acc + code * (code + 6849) + 29;
                         });
}

// if skipDoubleBonds is true, double bonds are skipped
// if doubleBondFlag is false, non-double bonds are skipped
void addBondsFromTemplate(const RWMol &templateMol, RWMol &newMol,
                          const std::map<unsigned, unsigned> &newAtomMap,
                          const boost::dynamic_bitset<> &isAtomInFragment,
                          bool skipDoubleBonds) {
  for (const auto templateBond : templateMol.bonds()) {
    bool isDoubleBond = (templateBond->getBondType() == Bond::DOUBLE);
    bool shouldProcessBond = (isDoubleBond ^ skipDoubleBonds);
    if (!shouldProcessBond ||
        !isAtomInFragment.test(templateBond->getBeginAtomIdx()) ||
        !isAtomInFragment.test(templateBond->getEndAtomIdx())) {
      continue;
    }
    auto ai1 = newAtomMap.at(templateBond->getBeginAtomIdx());
    auto ai2 = newAtomMap.at(templateBond->getEndAtomIdx());
    auto newBondIdx = newMol.addBond(ai1, ai2, templateBond->getBondType()) - 1;
    auto newBond = newMol.getBondWithIdx(newBondIdx);
    newBond->setBondDir(templateBond->getBondDir());
    if (isDoubleBond) {
      const auto &stereoAtoms = templateBond->getStereoAtoms();
      if (stereoAtoms.size() == 2) {
        newBond->setStereoAtoms(newAtomMap.at(stereoAtoms[0]),
                                newAtomMap.at(stereoAtoms[1]));
      }
    }
  }
}

void extractAtoms(const RWMol &templateMol, RWMol &newMol,
                  const INT_VECT &atomIndices) {
  boost::dynamic_bitset<> isAtomInFragment(templateMol.getNumAtoms());
  // key is atom index in template molecule
  std::map<unsigned int, unsigned int> newAtomMap;
  for (int ai : atomIndices) {
    isAtomInFragment.set(ai);
    auto a = templateMol.getAtomWithIdx(ai);
    newAtomMap[ai] = newMol.addAtom(a->copy(), true, true);
  }
  // add bonds from this fragment skipping double bonds
  addBondsFromTemplate(templateMol, newMol, newAtomMap, isAtomInFragment, true);
  // add bonds from this fragment skipping non-double bonds
  addBondsFromTemplate(templateMol, newMol, newAtomMap, isAtomInFragment,
                       false);
}
}  // namespace detail

static inline void convertMatchingToBondVect(
    std::vector<BondVector_t> &matching_bonds,
    const std::vector<MatchVectType> &matching_atoms) {
  for (const auto &matching_atom : matching_atoms) {
    matching_bonds.emplace_back();
    auto &mb = matching_bonds.back();  // current match
    // assume pattern is only one bond pattern
    auto a1 =
        static_cast<unsigned>(matching_atom[0].second);  // mol atom 1 index
    auto a2 =
        static_cast<unsigned>(matching_atom[1].second);  // mol atom 2 index
    mb.emplace_back(a1, a2);
  }
}

static void addResult(std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>
                          &res,  // const SignatureVector& resSignature,
                      const ROMol &mol, const BondVector_t &bonds_selected,
                      size_t maxCuts) {
#ifdef MMPA_DEBUG
  std::cout << res.size() + 1 << ": ";
#endif
  // loop through the bonds to delete. == deleteBonds()
  unsigned isotope = 0;
  std::map<unsigned, unsigned> isotope_track;
  std::vector<unsigned int> bondIndices;
  bondIndices.reserve(bonds_selected.size());
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  dummyLabels.reserve(bonds_selected.size());
  for (const auto &bi : bonds_selected) {
#ifdef MMPA_DEBUG
    {
      std::string symbol =
          mol.getAtomWithIdx(bonds_selected[bi].first)->getSymbol();
      int label = 0;
      mol.getAtomWithIdx(bonds_selected[bi].first)
          ->getPropIfPresent(common_properties::molAtomMapNumber, label);
      char a1[32];
      if (0 == label) {
        sprintf(a1, "\'%s\'", symbol.c_str(), label);
      } else {
        sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
      }
      symbol = mol.getAtomWithIdx(bonds_selected[bi].second)->getSymbol();
      label = 0;
      mol.getAtomWithIdx(bonds_selected[bi].second)
          ->getPropIfPresent(common_properties::molAtomMapNumber, label);
      char a2[32];
      if (0 == label) {
        sprintf(a2, "\'%s\'", symbol.c_str(), label);
      } else {
        sprintf(a2, "\'%s:%u\'", symbol.c_str(), label);
      }

      std::cout << "(" << bonds_selected[bi].first << a1 << ","
                << bonds_selected[bi].second << a2 << ") ";
    }
#endif
    ++isotope;
    const auto oldBond = mol.getBondBetweenAtoms(bi.first, bi.second);
    CHECK_INVARIANT(oldBond, "bond not found");
    bondIndices.push_back(oldBond->getIdx());
    dummyLabels.emplace_back(isotope, isotope);
  }
  std::unique_ptr<ROMol> fragmentedMol(
      MolFragmenter::fragmentOnBonds(mol, bondIndices, true, &dummyLabels));
  for (auto ai = mol.getNumAtoms(); ai < fragmentedMol->getNumAtoms(); ++ai) {
    auto a = fragmentedMol->getAtomWithIdx(ai);
    isotope = a->getIsotope();
    CHECK_INVARIANT(isotope, "isotope should be >0");
    a->setProp(common_properties::molAtomMapNumber, static_cast<int>(isotope));
    a->setIsotope(0);
    isotope_track[ai] = isotope;
  }
#ifdef MMPA_DEBUG
  std::cout << "\n";
#endif
  RWMOL_SPTR core, side_chains;  // core & side_chains output molecules

  if (bondIndices.size() == 1) {
    side_chains =
        RWMOL_SPTR(new RWMol(*fragmentedMol));  // output = '%s,%s,,%s.%s'
// DEBUG PRINT
#ifdef MMPA_DEBUG
// OK: std::cout<<res.size()+1<<" isotope="<< isotope <<","<<
// MolToSmiles(*side_chains, true) <<"\n";
#endif
  } else if (bondIndices.size() >= 2) {
    VECT_INT_VECT frags;
    MolOps::getMolFrags(*fragmentedMol, frags);

    // #check if its a valid triple or bigger cut.  matchObj = re.search(
    //'\*.*\*.*\*', f)
    //  check if exists a fragment with maxCut connection points (*.. *.. *)
    if (isotope >= 3) {
      bool valid = std::any_of(
          frags.begin(), frags.end(),
          [&isotope_track, &maxCuts](const INT_VECT &frag) {
            size_t nLabels = std::count_if(
                frag.begin(), frag.end(), [&isotope_track](int ai) {
                  return isotope_track.end() != isotope_track.find(ai);
                });
            return (nLabels >= maxCuts);
          });
      if (!valid) {
#ifdef MMPA_DEBUG
        std::cout << "isotope>=3: invalid fragments. fragment with maxCut "
                     "connection points not found"
                  << "\n";
#endif
        return;
      }
    }

    auto iCore = std::numeric_limits<size_t>::max();
    side_chains = RWMOL_SPTR(new RWMol);
    size_t maxAttachments = 0;
    for (size_t i = 0; i < frags.size(); i++) {
      size_t nAttachments = std::count_if(
          frags[i].begin(), frags[i].end(), [&isotope_track](int ai) {
            return isotope_track.end() != isotope_track.find(ai);
          });
      if (maxAttachments < nAttachments) {
        maxAttachments = nAttachments;
      }
      if (1 == nAttachments) {  // build side-chain set of molecules from
                                // selected fragment
        detail::extractAtoms(*fragmentedMol, *side_chains, frags[i]);
      } else {  // select the core fragment
// DEBUG PRINT
#ifdef MMPA_DEBUG
        if (iCore != -1) {
          std::cout << "Next CORE found. iCore=" << iCore << " New i=" << i
                    << " nAttachments=" << nAttachments << "\n";
        }
#endif
        if (nAttachments >= maxAttachments) {  // Choose a fragment with maximal
                                               // number of connection points as
                                               // a core
          iCore = i;
        }
      }
    }
    // build core molecule from selected fragment
    if (iCore != std::numeric_limits<size_t>::max()) {
      core = RWMOL_SPTR(new RWMol);
      detail::extractAtoms(*fragmentedMol, *core, frags[iCore]);
// DEBUG PRINT
#ifdef MMPA_DEBUG
// std::cout<<res.size()+1<<" isotope="<< isotope <<" "<< MolToSmiles(*core,
// true)<<", "<<MolToSmiles(*side_chains, true)<<"\n";
#endif
    }  // iCore != -1
  }
  // check for duplicates:
  bool resFound = false;
  for (const auto &r : res) {
    if (side_chains->getNumAtoms() == r.second->getNumAtoms() &&
        side_chains->getNumBonds() == r.second->getNumBonds() &&
        ((nullptr == core.get() && nullptr == r.first.get()) ||
         (nullptr != core.get() && nullptr != r.first.get() &&
          core->getNumAtoms() == r.first->getNumAtoms() &&
          core->getNumBonds() == r.first->getNumBonds()))) {
      // ToDo accurate check:
      // 1. compare hash code
      if (detail::computeMorganCodeHash(*side_chains) ==
              detail::computeMorganCodeHash(*r.second) &&
          (nullptr == core || detail::computeMorganCodeHash(*core) ==
                                  detail::computeMorganCodeHash(*r.first))) {
        // 2. final check to exclude hash collisions
        // We decided that it is not necessary to implement
        resFound = true;
        break;
      }
    }
  }
  if (!resFound) {
    // std::cerr << "**********************" << std::endl;
    // From rfrag.py
    // now change the labels on sidechains and core
    // to get the new labels, cansmi the dot-disconnected side chains
    // the first fragment in the side chains has attachment label 1, 2nd: 2,
    // 3rd: 3
    // then change the labels accordingly in the core
    std::map<unsigned int, int> canonicalAtomMaps;
    if (side_chains.get()) {
      RWMol tmp_side_chain(*(side_chains.get()));
      std::vector<int> oldMaps(tmp_side_chain.getNumAtoms(), 0);

      // clear atom labels (they are used in canonicalization)
      //  and move them to dummy storage
      for (auto at : tmp_side_chain.atoms()) {
        int label = 0;
        if (at->getAtomicNum() == 0 &&
            at->getPropIfPresent(common_properties::molAtomMapNumber, label)) {
          at->clearProp(common_properties::molAtomMapNumber);
          oldMaps[at->getIdx()] = label;
        }
      }

      const bool doIsomericSmiles = true;  // should this be false???
      auto smiles = MolToSmiles(tmp_side_chain, doIsomericSmiles);
      // std::cerr << "smiles: " << smiles << std::endl;

      // Get the canonical output order and use it to remap
      // the atom maps in the side chains
      // these will get reapplied to the core (if there is a core)
      const auto &ranks = tmp_side_chain.getProp<std::vector<unsigned int>>(
          common_properties::_smilesAtomOutputOrder);

      std::vector<std::pair<unsigned int, int>> rankedAtoms;

      for (size_t idx = 0; idx < ranks.size(); ++idx) {
        unsigned int atom_idx = ranks[idx];
        if (oldMaps[atom_idx] > 0) {
          const int label = oldMaps[atom_idx];
          // std::cerr << "atom_idx: " << atom_idx << " rank: " <<
          // ranks[atom_idx] <<
          //    " molAtomMapNumber: " << label << std::endl;
          rankedAtoms.emplace_back(idx, label);
        }
      }
      std::sort(rankedAtoms.begin(), rankedAtoms.end());
      int nextMap = 0;
      for (auto &rankedAtom : rankedAtoms) {
        if (canonicalAtomMaps.find(rankedAtom.second) ==
            canonicalAtomMaps.end()) {
          canonicalAtomMaps[rankedAtom.second] = ++nextMap;
        }
      }
    }

    // std::cerr << "======== Remap core " << std::endl;
    if (core.get()) {  // remap core if it exists
      for (auto at : core->atoms()) {
        int label = 0;
        if (at->getAtomicNum() == 0 &&
            at->getPropIfPresent(common_properties::molAtomMapNumber, label)) {
          // std::cerr << "remapping core: " << label << " :" <<
          // canonicalAtomMaps[label] <<
          //    std::endl;
          at->setProp(common_properties::molAtomMapNumber,
                      canonicalAtomMaps.at(label));
        }
      }
    }

    // std::cerr << "======== Remap side-chain " << std::endl;
    for (auto at : side_chains->atoms()) {
      int label = 0;
      if (at->getAtomicNum() == 0 &&
          at->getPropIfPresent(common_properties::molAtomMapNumber, label)) {
        // std::cerr << "remapping side chain: " << label << " :" <<
        // canonicalAtomMaps[label] << std::endl;
        at->setProp(common_properties::molAtomMapNumber,
                    canonicalAtomMaps.at(label));
      }
    }

    res.emplace_back(core, side_chains);  //
  }
#ifdef MMPA_DEBUG
  else {
    std::cout << res.size() + 1 << " --- DUPLICATE Result FOUND --- ri=" << ri
              << "\n";
  }
#endif
}

//=====================================================================
static inline void appendBonds(BondVector_t &bonds,
                               const BondVector_t &matching_bonds) {
  bonds.reserve(bonds.size() + matching_bonds.size());
  bonds.insert(bonds.end(), matching_bonds.begin(), matching_bonds.end());
}

static inline void processCuts(
    size_t i, size_t minCuts, size_t maxCuts,
    BondVector_t &bonds_selected,
    const std::vector<BondVector_t> &matching_bonds, const ROMol &mol,
    std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &res) {
  if (maxCuts < minCuts) {
    throw ValueErrorException("supplied maxCuts is less than minCuts");
  }

  if (minCuts == 0) {
    throw ValueErrorException("minCuts must be greater than 0");
  }

  for (size_t x = i; x < matching_bonds.size(); x++) {
    appendBonds(bonds_selected, matching_bonds[x]);
    if (bonds_selected.size() >= minCuts) {
      addResult(res, mol, bonds_selected, maxCuts);
    }
    if (bonds_selected.size() < maxCuts) {
      processCuts(x + 1, minCuts, maxCuts, bonds_selected, matching_bonds, mol,
                  res);
    }

    bonds_selected.pop_back();
  }
}

//=====================================================================
// Public API implementation:
//=====================================================================

bool fragmentMol(const ROMol &mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &res,
                 unsigned int maxCuts, unsigned int maxCutBonds,
                 const std::string &pattern) {
  return fragmentMol(mol, res, 1, maxCuts, maxCutBonds, pattern);
}

bool fragmentMol(const ROMol &mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &res,
                 unsigned int minCuts, unsigned int maxCuts,
                 unsigned int maxCutBonds, const std::string &pattern) {
#ifdef MMPA_DEBUG
  for (size_t i = 0; i < mol.getNumAtoms(); i++) {
    std::string symbol = mol.getAtomWithIdx(i)->getSymbol();
    int label = 0;
    mol.getAtomWithIdx(i)->getPropIfPresent(common_properties::molAtomMapNumber,
                                            label);
    char a1[32];
    if (0 == label) {
      sprintf(a1, "\'%s\'", symbol.c_str(), label);
    } else {
      sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
    }
    std::cout << "Atom " << i << ": " << a1;  //<<" Bonds:";
    std::cout << "\n";
  }
#endif

  res.clear();
  std::unique_ptr<const ROMol> smarts(
      static_cast<ROMol *>(SmartsToMol(pattern)));
  std::vector<MatchVectType>
      matching_atoms;  // one bond per match ! with default pattern
  auto total = SubstructMatch(mol, *smarts, matching_atoms);
#ifdef MMPA_DEBUG
  std::cout << "total substructs =" << total
            << "\nmatching bonds (atom1, atom2):\n";
#endif
  if (0 == total) {  // Not found.  Return empty set of molecules
    return false;
  }
#ifdef MMPA_DEBUG
  for (size_t i = 0; i < matching_atoms.size(); i++) {
    std::string symbol =
        mol.getAtomWithIdx(matching_atoms[i][0].second)->getSymbol();
    int label = 0;
    mol.getAtomWithIdx(matching_atoms[i][0].second)
        ->getPropIfPresent(common_properties::molAtomMapNumber, label);
    char a1[32];
    if (0 == label) {
      sprintf(a1, "\'%s\'", symbol.c_str(), label);
    } else {
      sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
    }
    symbol = mol.getAtomWithIdx(matching_atoms[i][1].second)->getSymbol();
    label = 0;
    mol.getAtomWithIdx(matching_atoms[i][1].second)
        ->getPropIfPresent(common_properties::molAtomMapNumber, label);
    char a2[32];
    if (0 == label) {
      sprintf(a2, "\'%s\'", symbol.c_str(), label);
    } else {
      sprintf(a2, "\'%s:%u\'", symbol.c_str(), label);
    }

    std::cout << i << ": (" << matching_atoms[i][0].second << a1 << ","
              << matching_atoms[i][1].second << a2 << ") \n";
  }
#endif

  std::vector<BondVector_t> matching_bonds;  // List of matched query's bonds
  convertMatchingToBondVect(matching_bonds, matching_atoms);
  if (matching_bonds.size() > maxCutBonds) {
    return false;
  }
#ifdef MMPA_DEBUG
  std::cout << "total matching_bonds = " << matching_bonds.size() << "\n";
#endif

  // loop to generate every cut in the molecule
  BondVector_t bonds_selected;
  processCuts(0, minCuts, maxCuts, bonds_selected, matching_bonds, mol, res);
  return true;
}

bool fragmentMol(const ROMol &mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &res,
                 const std::vector<unsigned int> &bondsToCut,
                 unsigned int minCuts, unsigned int maxCuts) {
  std::vector<BondVector_t> matching_bonds;  // List of matched query's bonds

  for (auto i : bondsToCut) {
    const auto bond = mol.getBondWithIdx(i);
    BondVector_t bonds;
    auto a1 = bond->getBeginAtomIdx();
    auto a2 = bond->getEndAtomIdx();
    bonds.emplace_back(a1, a2);
    matching_bonds.push_back(std::move(bonds));
  }

  // loop to generate every cut in the molecule
  BondVector_t bonds_selected;
  processCuts(0, minCuts, maxCuts, bonds_selected, matching_bonds, mol, res);
  return true;
}
}  // namespace MMPA
}  // namespace RDKit
