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
#include <math.h>

#include "../MolOps.h"
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../Substruct/SubstructMatch.h"
#include <GraphMol/new_canon.h>
#include <GraphMol/MolOps.h>
#include "MMPA.h"

//#define _DEBUG // enable debug info output

namespace RDKit {
namespace MMPA {

typedef std::vector<std::pair<unsigned, unsigned> >
    BondVector_t;  // pair of BeginAtomIdx, EndAtomIdx

static inline unsigned long long computeMorganCodeHash(const ROMol& mol) {
  size_t nv = mol.getNumAtoms();
  size_t ne = mol.getNumBonds();
  std::vector<unsigned long> currCodes(nv);
  std::vector<unsigned long> prevCodes(nv);
  size_t nIterations = mol.getNumBonds();
  if (nIterations > 5) nIterations = 5;

  for (unsigned ai = 0; ai < nv; ai++) {
    const Atom& a = *mol.getAtomWithIdx(ai);
    unsigned atomCode = a.getAtomicNum();
    atomCode |= a.getIsotope() << 8;
    atomCode |= a.getFormalCharge() << 16;
    atomCode |= (a.getIsAromatic() ? 1 : 0) << 30;
    currCodes[ai] = atomCode;
  }

  for (size_t iter = 0; iter < nIterations; iter++) {
    for (size_t i = 0; i < nv; i++) prevCodes[i] = currCodes[i];

    for (size_t bi = 0; bi < ne; bi++) {
      const Bond* bond = mol.getBondWithIdx(bi);
      unsigned order = bond->getBondType();
      unsigned atom1 = bond->getBeginAtomIdx();
      unsigned atom2 = bond->getEndAtomIdx();
      unsigned v1 = prevCodes[atom1];
      unsigned v2 = prevCodes[atom2];

      currCodes[atom1] += v2 * v2 + (v2 + 23) * (order + 1721);
      currCodes[atom2] += v1 * v1 + (v1 + 23) * (order + 1721);
    }
  }

  unsigned long long result = 0;
  for (unsigned ai = 0; ai < nv; ai++) {
    unsigned long code = currCodes[ai];
    result += code * (code + 6849) + 29;
  }
  return result;
}

static inline void convertMatchingToBondVect(
    std::vector<BondVector_t>& matching_bonds,
    const std::vector<MatchVectType>& matching_atoms, const ROMol& mol) {
  RDUNUSED_PARAM(mol);
  for (std::vector<MatchVectType>::const_iterator m = matching_atoms.begin();
       m != matching_atoms.end(); ++m) {
    matching_bonds.push_back(BondVector_t());
    BondVector_t& mb = matching_bonds.back();  // current match
    // assume patern is only one bond pattern
    unsigned a1 = (unsigned)(*m)[0].second;  // mol atom 1 index
    unsigned a2 = (unsigned)(*m)[1].second;  // mol atom 2 index
    mb.push_back(std::pair<unsigned, unsigned>(a1, a2));
  }
}

static void addResult(std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR> >&
                          res,  // const SignatureVector& resSignature,
                      const ROMol& mol,
                      const BondVector_t& bonds_selected, size_t maxCuts) {
#ifdef _DEBUG
  std::cout << res.size() + 1 << ": ";
#endif
  RWMol em(mol);
  // loop through the bonds to delete. == deleteBonds()
  unsigned isotope = 0;
  std::map<unsigned, unsigned> isotope_track;
  for (size_t i = 0; i < bonds_selected.size(); i++) {
#ifdef _DEBUG
    {
      std::string symbol =
          em.getAtomWithIdx(bonds_selected[i].first)->getSymbol();
      int label = 0;
      em.getAtomWithIdx(bonds_selected[i].first)
          ->getPropIfPresent(common_properties::molAtomMapNumber, label);
      char a1[32];
      if (0 == label)
        sprintf(a1, "\'%s\'", symbol.c_str(), label);
      else
        sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
      symbol = em.getAtomWithIdx(bonds_selected[i].second)->getSymbol();
      label = 0;
      em.getAtomWithIdx(bonds_selected[i].second)
          ->getPropIfPresent(common_properties::molAtomMapNumber, label);
      char a2[32];
      if (0 == label)
        sprintf(a2, "\'%s\'", symbol.c_str(), label);
      else
        sprintf(a2, "\'%s:%u\'", symbol.c_str(), label);

      std::cout << "(" << bonds_selected[i].first << a1 << ","
                << bonds_selected[i].second << a2 << ") ";
    }
#endif
    isotope += 1;
    // remove the bond
    em.removeBond(bonds_selected[i].first, bonds_selected[i].second);

    // now add attachement points and set attachment point lables
    Atom* a = new Atom(0);
    a->setProp(common_properties::molAtomMapNumber, (int)isotope);
    unsigned newAtomA = em.addAtom(a, true, true);
    em.addBond(bonds_selected[i].first, newAtomA, Bond::SINGLE);
    a = new Atom(0);
    a->setProp(common_properties::molAtomMapNumber, (int)isotope);
    unsigned newAtomB = em.addAtom(a, true, true);
    em.addBond(bonds_selected[i].second, newAtomB, Bond::SINGLE);

    // keep track of where to put isotopes
    isotope_track[newAtomA] = isotope;
    isotope_track[newAtomB] = isotope;
  }
#ifdef _DEBUG
  std::cout << "\n";
#endif
  RWMOL_SPTR core, side_chains;  // core & side_chains output molecules

  if (isotope == 1) {
    side_chains = RWMOL_SPTR(new RWMol(em));  // output = '%s,%s,,%s.%s'
// DEBUG PRINT
#ifdef _DEBUG
// OK: std::cout<<res.size()+1<<" isotope="<< isotope <<","<<
// MolToSmiles(*side_chains, true) <<"\n";
#endif
  } else if (isotope >= 2) {
    std::vector<std::vector<int> > frags;
    unsigned int nFrags = MolOps::getMolFrags(em, frags);

    //#check if its a valid triple or bigger cut.  matchObj = re.search(
    //'\*.*\*.*\*', f)
    // check if exists a fragment with maxCut connection points (*.. *.. *)
    if (isotope >= 3) {
      bool valid = false;
      for (size_t i = 0; i < nFrags; i++) {
        unsigned nLabels = 0;
        for (size_t ai = 0; ai < frags[i].size(); ai++) {
          if (isotope_track.end() !=
              isotope_track.find(frags[i][ai]))  // new added atom
            ++nLabels;                           // found connection point
        }
        if (nLabels >=
            maxCuts) {  // looks like it should be selected as core !  ??????
          valid = true;
          break;
        }
      }
      if (!valid) {
#ifdef _DEBUG
        std::cout << "isotope>=3: invalid fragments. fragment with maxCut "
                     "connection points not found"
                  << "\n";
#endif
        return;
      }
    }

    size_t iCore = std::numeric_limits<size_t>::max();
    side_chains = RWMOL_SPTR(new RWMol);
    std::map<unsigned, unsigned>
        visitedBonds;  // key is bond index in source molecule
    unsigned maxAttachments = 0;
    for (size_t i = 0; i < frags.size(); i++) {
      unsigned nAttachments = 0;
      for (size_t ai = 0; ai < frags[i].size(); ai++) {
        if (isotope_track.end() !=
            isotope_track.find(
                frags[i][ai]))  // == if(a->hasProp("molAtomMapNumber"))
          ++nAttachments;
      }
      if (maxAttachments < nAttachments) maxAttachments = nAttachments;
      if (1 == nAttachments) {  // build side-chain set of molecules from
                                // selected fragment
        std::map<unsigned, unsigned>
            newAtomMap;  // key is atom index in source molecule
        for (size_t ai = 0; ai < frags[i].size(); ai++) {
          Atom* a = em.getAtomWithIdx(frags[i][ai]);
          newAtomMap[frags[i][ai]] =
              side_chains->addAtom(a->copy(), true, true);
        }
        // add all bonds from this fragment
        for (size_t ai = 0; ai < frags[i].size(); ai++) {
          Atom* a = em.getAtomWithIdx(frags[i][ai]);
          ROMol::OEDGE_ITER beg, end;
          for (boost::tie(beg, end) = em.getAtomBonds(a); beg != end; ++beg) {
            const BOND_SPTR bond = em[*beg];
            if (newAtomMap.end() == newAtomMap.find(bond->getBeginAtomIdx()) ||
                newAtomMap.end() == newAtomMap.find(bond->getEndAtomIdx()) ||
                visitedBonds.end() != visitedBonds.find(bond->getIdx()))
              continue;
            unsigned ai1 = newAtomMap[bond->getBeginAtomIdx()];
            unsigned ai2 = newAtomMap[bond->getEndAtomIdx()];
            unsigned bi = side_chains->addBond(ai1, ai2, bond->getBondType());
            visitedBonds[bond->getIdx()] = bi;
          }
        }
      } else {  // select the core fragment
// DEBUG PRINT
#ifdef _DEBUG
        if (iCore != -1)
          std::cout << "Next CORE found. iCore=" << iCore << " New i=" << i
                    << " nAttachments=" << nAttachments << "\n";
#endif
        if (nAttachments >= maxAttachments)  // Choose a fragment with maximal
                                             // number of connection points as a
                                             // core
          iCore = i;
      }
    }
    // build core molecule from selected fragment
    if (iCore != std::numeric_limits<size_t>::max()) {
      core = RWMOL_SPTR(new RWMol);
      visitedBonds.clear();
      std::map<unsigned, unsigned>
          newAtomMap;  // key is atom index in source molecule
      for (size_t i = 0; i < frags[iCore].size(); i++) {
        unsigned ai = frags[iCore][i];
        Atom* a = em.getAtomWithIdx(ai);
        newAtomMap[ai] = core->addAtom(a->copy(), true, true);
      }
      // add all bonds from this fragment
      for (size_t ai = 0; ai < frags[iCore].size(); ai++) {
        Atom* a = em.getAtomWithIdx(frags[iCore][ai]);
        ROMol::OEDGE_ITER beg, end;
        for (boost::tie(beg, end) = em.getAtomBonds(a); beg != end; ++beg) {
          const BOND_SPTR bond = em[*beg];
          if (newAtomMap.end() == newAtomMap.find(bond->getBeginAtomIdx()) ||
              newAtomMap.end() == newAtomMap.find(bond->getEndAtomIdx()) ||
              visitedBonds.end() != visitedBonds.find(bond->getIdx()))
            continue;
          unsigned ai1 = newAtomMap[bond->getBeginAtomIdx()];
          unsigned ai2 = newAtomMap[bond->getEndAtomIdx()];
          unsigned bi = core->addBond(ai1, ai2, bond->getBondType());
          visitedBonds[bond->getIdx()] = bi;
        }
      }
// DEBUG PRINT
#ifdef _DEBUG
// std::cout<<res.size()+1<<" isotope="<< isotope <<" "<< MolToSmiles(*core,
// true)<<", "<<MolToSmiles(*side_chains, true)<<"\n";
#endif
    }  // iCore != -1
  }
  // check for duplicates:
  bool resFound = false;
  size_t ri = 0;
  for (ri = 0; ri < res.size(); ri++) {
    const std::pair<ROMOL_SPTR, ROMOL_SPTR>& r = res[ri];
    if (side_chains->getNumAtoms() == r.second->getNumAtoms() &&
        side_chains->getNumBonds() == r.second->getNumBonds() &&
        ((NULL == core.get() && NULL == r.first.get()) ||
         (NULL != core.get() && NULL != r.first.get() &&
          core->getNumAtoms() == r.first->getNumAtoms() &&
          core->getNumBonds() == r.first->getNumBonds()))) {
      // ToDo accurate check:
      // 1. compare hash code
      if (computeMorganCodeHash(*side_chains) ==
              computeMorganCodeHash(*r.second) &&
          (NULL == core ||
           computeMorganCodeHash(*core) == computeMorganCodeHash(*r.first))) {
        // 2. final check to exclude hash collisions
        // We decided that it does not neccessary to implement
        resFound = true;
        break;
      }
    }
  }
  if (!resFound) {
    //std::cerr << "**********************" << std::endl;
    // From rfrag.py
    // now change the labels on sidechains and core
    // to get the new labels, cansmi the dot-disconnected side chains
    // the first fragment in the side chains has attachment label 1, 2nd: 2, 3rd: 3
    // then change the labels accordingly in the core
    std::map<unsigned int, int> canonicalAtomMaps;
    if( side_chains.get() ) {
      RWMol tmp_side_chain(*(side_chains.get()));
      std::vector<int> oldMaps(tmp_side_chain.getNumAtoms(), 0);
      
      // clear atom labels (they are used in canonicalization)
      //  and move them to dummy storage
      for (ROMol::AtomIterator at = tmp_side_chain.beginAtoms(); at != tmp_side_chain.endAtoms();
           ++at) {
        int label = 0;
        if ((*at)->getPropIfPresent(common_properties::molAtomMapNumber, label) ) {
          (*at)->clearProp(common_properties::molAtomMapNumber);
          oldMaps[(*at)->getIdx()] = label;
        }
      }

      const bool doIsomericSmiles = true; // should this be false???
      std::string smiles = MolToSmiles(tmp_side_chain, doIsomericSmiles);
      //std::cerr << "smiles: " << smiles << std::endl;

      // Get the canonical output order and use it to remap
      //  the atom maps int the side chains
      //  these will get reapplied to the core (if there is a core)
      const std::vector<unsigned int> &ranks = tmp_side_chain.getProp<
        std::vector<unsigned int> >(
          common_properties::_smilesAtomOutputOrder);
      
      std::vector<std::pair<unsigned int, int> > rankedAtoms;

      for(size_t idx=0;idx<ranks.size();++idx) {
        unsigned int atom_idx = ranks[idx];
        if(oldMaps[atom_idx] >0) {
          const int label = oldMaps[atom_idx];
          //std::cerr << "atom_idx: " << atom_idx << " rank: " << ranks[atom_idx] <<
          //    " molAtomMapNumber: " << label << std::endl;
          rankedAtoms.push_back(std::make_pair(idx, label));
        }
      }
      std::sort(rankedAtoms.begin(), rankedAtoms.end());
      int nextMap = 0;
      for(size_t i=0;i<rankedAtoms.size();++i) {
        if(canonicalAtomMaps.find(rankedAtoms[i].second) == canonicalAtomMaps.end()) {
          //std::cerr << "Remapping: " << rankedAtoms[i].second << " " << " to " << (i+1) <<
          //    std::endl;
          canonicalAtomMaps[rankedAtoms[i].second] = ++nextMap;
        }
      }
    }
    
    //std::cerr << "======== Remap core " << std::endl;
    if( core.get() ) { // remap core if it exists
      for (ROMol::AtomIterator at = core->beginAtoms(); at != core->endAtoms();
           ++at) {
        int label = 0;
        if ((*at)->getPropIfPresent(common_properties::molAtomMapNumber, label) ) {
          //std::cerr << "remapping core: " << label << " :" << canonicalAtomMaps[label] <<
          //    std::endl;
          (*at)->setProp(common_properties::molAtomMapNumber, canonicalAtomMaps[label]);
        }
      }
    }
    
    //std::cerr << "======== Remap side-chain " << std::endl;
    for (ROMol::AtomIterator at = side_chains->beginAtoms(); at != side_chains->endAtoms();
         ++at) {
      int label = 0;
      if ((*at)->getPropIfPresent(common_properties::molAtomMapNumber, label) ) {
        //std::cerr << "remapping side chain: " << label << " :" << 
        // canonicalAtomMaps[label] << std::endl;
        (*at)->setProp(common_properties::molAtomMapNumber, canonicalAtomMaps[label]);
      }
    }

    
    res.push_back(std::pair<ROMOL_SPTR, ROMOL_SPTR>(core, side_chains));  //
  }
#ifdef _DEBUG
  else
    std::cout << res.size() + 1 << " --- DUPLICATE Result FOUND --- ri=" << ri
              << "\n";
#endif
}

//=====================================================================
static inline void appendBonds(BondVector_t& bonds,
                               const BondVector_t& matching_bonds) {
  for (BondVector_t::const_iterator b = matching_bonds.begin();
       b != matching_bonds.end(); ++b)
    bonds.push_back(*b);
}

static inline void processCuts(
    size_t i, size_t minCuts, size_t maxCuts, BondVector_t& bonds_selected,
    const std::vector<BondVector_t>& matching_bonds, const ROMol& mol,
    std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR> >& res) {
  if(maxCuts < minCuts)
    throw ValueErrorException("supplied maxCuts is less than minCuts");
        
  if(minCuts==0)
    throw ValueErrorException("minCuts must be greater than 0");

  for (size_t x = i; x < matching_bonds.size(); x++) {
    appendBonds(bonds_selected, matching_bonds[x]);
    if(bonds_selected.size() >= minCuts) {
      addResult(res, mol, bonds_selected, maxCuts);
    }
    if (bonds_selected.size() < maxCuts) {
      processCuts(x + 1, minCuts, maxCuts, bonds_selected, matching_bonds, mol, res);
    }
    
    bonds_selected.pop_back();
  }
}

//=====================================================================
// Public API implementation:
//=====================================================================

bool fragmentMol(const ROMol& mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR> >& res,
                 unsigned int maxCuts,
                 unsigned int maxCutBonds,
                 const std::string& pattern) {
  return fragmentMol(mol, res, 1, maxCuts, maxCutBonds, pattern);
}

bool fragmentMol(const ROMol& mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR> >& res,
                 unsigned int minCuts,
                 unsigned int maxCuts,
                 unsigned int maxCutBonds,
                 const std::string& pattern) {
#ifdef _DEBUG
  for (size_t i = 0; i < mol.getNumAtoms(); i++) {
    std::string symbol = mol.getAtomWithIdx(i)->getSymbol();
    int label = 0;
    mol.getAtomWithIdx(i)->getPropIfPresent(common_properties::molAtomMapNumber,
                                            label);
    char a1[32];
    if (0 == label)
      sprintf(a1, "\'%s\'", symbol.c_str(), label);
    else
      sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
    std::cout << "Atom " << i << ": " << a1;  //<<" Bonds:";
    std::cout << "\n";
  }
#endif

  res.clear();
  std::auto_ptr<const ROMol> smarts((const ROMol*)SmartsToMol(pattern));
  std::vector<MatchVectType>
      matching_atoms;  // one bond per match ! with default pattern
  unsigned int total = SubstructMatch(mol, *smarts, matching_atoms);
#ifdef _DEBUG
  std::cout << "total substructs =" << total
            << "\nmatching bonds (atom1, atom2):\n";
#endif
  if (0 == total)  // Not found.  Return empty set of molecules
    return false;
#ifdef _DEBUG
  for (size_t i = 0; i < matching_atoms.size(); i++) {
    std::string symbol =
        mol.getAtomWithIdx(matching_atoms[i][0].second)->getSymbol();
    int label = 0;
    mol.getAtomWithIdx(matching_atoms[i][0].second)
        ->getPropIfPresent(common_properties::molAtomMapNumber, label);
    char a1[32];
    if (0 == label)
      sprintf(a1, "\'%s\'", symbol.c_str(), label);
    else
      sprintf(a1, "\'%s:%u\'", symbol.c_str(), label);
    symbol = mol.getAtomWithIdx(matching_atoms[i][1].second)->getSymbol();
    label = 0;
    mol.getAtomWithIdx(matching_atoms[i][1].second)
        ->getPropIfPresent(common_properties::molAtomMapNumber, label);
    char a2[32];
    if (0 == label)
      sprintf(a2, "\'%s\'", symbol.c_str(), label);
    else
      sprintf(a2, "\'%s:%u\'", symbol.c_str(), label);

    std::cout << i << ": (" << matching_atoms[i][0].second << a1 << ","
              << matching_atoms[i][1].second << a2 << ") \n";
  }
#endif

  std::vector<BondVector_t> matching_bonds;  // List of matched query's bonds
  convertMatchingToBondVect(matching_bonds, matching_atoms, mol);
  if (matching_bonds.size() > maxCutBonds) return false;
#ifdef _DEBUG
  std::cout << "total matching_bonds = " << matching_bonds.size() << "\n";
#endif

  // loop to generate every cut in the molecule
  BondVector_t bonds_selected;
  processCuts(0, minCuts, maxCuts, bonds_selected, matching_bonds, mol, res);
  return true;
}

bool fragmentMol(const ROMol& mol,
                 std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR> >& res,
                 const std::vector<unsigned int>& bondsToCut,
                 unsigned int minCuts,
                 unsigned int maxCuts) {
  std::vector<BondVector_t> matching_bonds;  // List of matched query's bonds

  BOOST_FOREACH(unsigned int i, bondsToCut) {
    const Bond *bond = mol.getBondWithIdx(i);
    BondVector_t bonds;
    unsigned int a1 = bond->getBeginAtomIdx();
    unsigned int a2 = bond->getEndAtomIdx();
    bonds.push_back( std::make_pair(a1, a2) );
    matching_bonds.push_back(bonds);
  }

  // loop to generate every cut in the molecule
  BondVector_t bonds_selected;
  processCuts(0, minCuts, maxCuts, bonds_selected, matching_bonds, mol, res);
  return true;
}

}
}  // namespace RDKit
