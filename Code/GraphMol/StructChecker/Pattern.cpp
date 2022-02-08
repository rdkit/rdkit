//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#if defined(__CYGWIN__) && !defined(_GNU_SOURCE)
// -std=c++11 doesn't declare strtok_r
#define _GNU_SOURCE
#endif

#include <map>
#include "../QueryAtom.h"
#include "Pattern.h"

// define snprintf for msvc
#if _MSC_VER
#if _MSC_VER < 1900
#define snprintf _snprintf
#endif
#endif

namespace RDKit {
namespace StructureCheck {

static void verboseAtom(const ROMol &mol, int atom_idx,
                        std::vector<Neighbourhood> &neighbours) {
  char nstr[512] = "";
  for (unsigned int j : neighbours[atom_idx].Atoms) {
    const Atom &a = *mol.getAtomWithIdx(j);
    snprintf(nstr + strlen(nstr), sizeof(nstr) - strlen(nstr), ", %s",
             a.getSymbol().c_str());
    if (0 != a.getFormalCharge())
      snprintf(nstr + strlen(nstr), sizeof(nstr) - strlen(nstr), "%+d",
               a.getFormalCharge());
  }
  BOOST_LOG(rdInfoLog) << "* CheckAtom idx=" << atom_idx << ", "
                       << mol.getAtomWithIdx(atom_idx)->getSymbol()
                       << " neighbours=" << neighbours[atom_idx].Atoms.size()
                       << nstr << "\n";
}

RDKit::Bond::BondType convertBondType(AABondType bt) {
  static const RDKit::Bond::BondType rdbt[] = {
      RDKit::Bond::UNSPECIFIED,  // BT_NONE = 0,
      RDKit::Bond::SINGLE,      RDKit::Bond::DOUBLE,
      RDKit::Bond::TRIPLE,      RDKit::Bond::AROMATIC,
      RDKit::Bond::UNSPECIFIED,  // ??? SINGLE_DOUBLE = 5,
      RDKit::Bond::ONEANDAHALF,  // SINGLE_AROMATIC = 6,
      RDKit::Bond::TWOANDAHALF,  // DOUBLE_AROMATIC = 7,
  };
  return (bt <= DOUBLE_AROMATIC) ? rdbt[bt] : RDKit::Bond::OTHER;  // ??
}

AABondType convertBondType(RDKit::Bond::BondType rdbt) {
  static const AABondType bt[] = {
      BT_NONE,         SINGLE, DOUBLE, TRIPLE,
      BT_NONE,          // QUADRUPLE,
      BT_NONE,          // QUINTUPLE,
      BT_NONE,          // HEXTUPLE,
      SINGLE_AROMATIC,  // ONEANDAHALF
      DOUBLE_AROMATIC,  // TWOANDAHALF
      BT_NONE,          // THREEANDAHALF,
      BT_NONE,          // FOURANDAHALF,
      BT_NONE,          // FIVEANDAHALF,
      AROMATIC,
  };
  return (rdbt <= RDKit::Bond::AROMATIC) ? bt[rdbt] : BT_NONE;  //??
}

/*
 * Checks if the ligand is compatible with the bond and the atom.
 */
bool LigandMatches(const Atom &a, const Bond &b, const Ligand &l,
                   bool use_charge) {
  if (l.BondType != ANY_BOND && !isBondTypeMatch(b, l.BondType)) return false;
  if (l.Radical != ANY_RADICAL && a.getNumRadicalElectrons() != l.Radical)
    return false;
  if ((l.Charge != ANY_CHARGE || use_charge) && l.Charge != ANY_CHARGE &&
      a.getFormalCharge() != l.Charge)
    return false;
  return (AtomSymbolMatch(a.getSymbol(), l.AtomSymbol));
}

static void applyAtomSymbolList(RWMol &mol, const std::string symbol,
                                Atom *&atom) {
  // replace mol atom with Ligrad.AtomSymbol
  // single atom symbol:
  if (std::string::npos == symbol.find(',')) {
    atom->setAtomicNum(getAtomicNumber(symbol));
    return;
  }
  // comma separated list of atom symbols:
  char *context;
#ifdef WIN32
#define strtok_r strtok_s  // thread safe strtok()
#endif
  const char *atsym = symbol.c_str();
  char buf[512];
  char *tokp;
  strcpy(buf, atsym);
  tokp = strtok_r(buf, ",", &context);
  atom->setAtomicNum(getAtomicNumber(tokp));
  QueryAtom a(*atom);
  a.setQuery(makeAtomNumQuery(getAtomicNumber(tokp)));
  //  QueryAtom a(atom);  // copy all mol's atom properties to keep its
  for (tokp = strtok_r((char *)nullptr, ",", &context); tokp;
       tokp = strtok_r((char *)nullptr, ",", &context)) {
    a.expandQuery(makeAtomNumQuery(getAtomicNumber(tokp)),
                  Queries::COMPOSITE_OR);
  }
  unsigned int idx = atom->getIdx();
  mol.replaceAtom(idx, &a);
  atom = mol.getAtomWithIdx(idx);
}

struct AtomNeighbor {
  unsigned AtomIdx;  // the molecule's atom
  unsigned BondIdx;  // the molecule's bond
  AtomNeighbor(unsigned a, unsigned b) : AtomIdx(a), BondIdx(b) {}
};

bool TransformAugmentedAtoms(
    RWMol &mol,
    const std::vector<std::pair<AugmentedAtom, AugmentedAtom>> &aapair,
    bool verbose) {
  bool transformed = false;
  // find neighbours:
  std::vector<std::vector<AtomNeighbor>> neighbours(mol.getNumAtoms());
  std::vector<bool> bondInRing(mol.getNumBonds());
  std::vector<bool> atomInRing(mol.getNumAtoms());
  std::vector<std::pair<unsigned, unsigned>> bondsToRemove;

  std::vector<Neighbourhood> neighbours_d(mol.getNumAtoms());
  SetupNeighbourhood(mol, neighbours_d);

  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    unsigned a1 = bond->getEndAtomIdx();
    unsigned a2 = bond->getBeginAtomIdx();
    neighbours[a1].push_back(AtomNeighbor(a2, i));
    neighbours[a2].push_back(AtomNeighbor(a1, i));
    bondInRing[i] = false;                    // init
    atomInRing[a1] = atomInRing[a2] = false;  // init
  }

  const RingInfo::VECT_INT_VECT &rings = mol.getRingInfo()->bondRings();
  for (const auto &ring : rings) {
    for (int j : ring) {
      bondInRing[j] = true;
      const Bond *bond = mol.getBondWithIdx(j);
      unsigned a1 = bond->getEndAtomIdx();
      unsigned a2 = bond->getBeginAtomIdx();
      atomInRing[a1] = atomInRing[a2] = true;
    }
  }
  // find matched atoms and replace
  for (size_t i = 0; i < aapair.size(); i++) {
    // AAFix():
    for (unsigned j = 0; j < mol.getNumAtoms(); j++) {
      const Atom *atom = mol.getAtomWithIdx(j);
      const AugmentedAtom &aa1 = aapair[i].first;
      /* a lot of prints
            if(verbose) {
               verboseAtom(mol, j, neighbours_d);
            }
      */
      // check if atom match to aapair[i].first atom and topology
      if (neighbours[j].size() == aa1.Ligands.size() &&
          (ANY_CHARGE == aa1.Charge || atom->getFormalCharge() == aa1.Charge) &&
          (ANY_RADICAL == aa1.Radical ||
           atom->getNumRadicalElectrons() == aa1.Radical) &&
          (TP_NONE == aa1.Topology || (RING == aa1.Topology && atomInRing[j]) ||
           (CHAIN == aa1.Topology && !atomInRing[j])) &&
          AtomSymbolMatch(atom->getSymbol(), aa1.AtomSymbol)) {
        // RecMatch():
        //        if ( ! RecMatch(mol, j, aa1, neighbours_d, verbose))
        //          continue;  // go to next molecule's atom
        unsigned matched = 0;
        const AugmentedAtom &aa2 = aapair[i].second;
        std::vector<bool> visited(aa2.Ligands.size());
        std::vector<unsigned> match_a(aa2.Ligands.size());
        std::vector<unsigned> match_b(aa2.Ligands.size());
        for (size_t l = 0; l < aa2.Ligands.size(); l++) {
          visited[l] = false;
          match_a[l] = -1;
          match_b[l] = -1;
        }

        for (size_t k = 0; k < neighbours[j].size(); k++) {
          const Atom *nbrAtom = mol.getAtomWithIdx(neighbours[j][k].AtomIdx);
          const Bond *nbrBond = mol.getBondWithIdx(neighbours[j][k].BondIdx);
          for (size_t l = 0; l < aa1.Ligands.size(); l++)
            if (!visited[l]) {
              const Ligand &ligand = aa1.Ligands[l];
              if ((0 == ligand.SubstitutionCount ||
                   neighbours[j].size() == ligand.SubstitutionCount) &&
                  (ANY_CHARGE == ligand.Charge ||
                   nbrAtom->getFormalCharge() == ligand.Charge) &&
                  (ANY_RADICAL == ligand.Radical ||
                   nbrAtom->getNumRadicalElectrons() == ligand.Radical) &&
                  ((ANY_BOND == ligand.BondType ||
                    isBondTypeMatch(*nbrBond, ligand.BondType)) ||
                   (SINGLE_DOUBLE == ligand.BondType &&
                    (RDKit::Bond::SINGLE == nbrBond->getBondType() ||
                     RDKit::Bond::DOUBLE == nbrBond->getBondType()))) &&
                  AtomSymbolMatch(nbrAtom->getSymbol(), ligand.AtomSymbol)) {
                matched++;
                match_a[l] = neighbours[j][k].AtomIdx;
                match_b[l] = neighbours[j][k].BondIdx;
                visited[l] = true;
                break;
              }
            }
        }
        if (matched != neighbours[j].size())
          continue;  // go to next molecule's atom

        // Replace Augmented Atom. TransformAA()
        if (verbose)
          BOOST_LOG(rdInfoLog)
              << "Replace " << LogNeighbourhood(mol, j, neighbours_d)
              << " with AAPair i=" << i << " " << aapair[i].first.AtomSymbol
              << " -> " << aapair[i].second.AtomSymbol << "\n";
        // change central atom
        Atom *a = mol.getAtomWithIdx(j);
        if (aa1.AtomSymbol != aa2.AtomSymbol)
          applyAtomSymbolList(mol, aa2.AtomSymbol, a);
        if (a->getFormalCharge() != aa2.Charge) a->setFormalCharge(aa2.Charge);
        if (a->getNumRadicalElectrons() != aa2.Radical)
          a->setNumRadicalElectrons(aa2.Radical);
        // change ligand atoms and their bonds with central atom
        for (size_t l = 0; l < aa2.Ligands.size(); l++) {
          Atom *al = mol.getAtomWithIdx(match_a[l]);
          const Ligand &ligand = aa2.Ligands[l];
          {
            // AVALON pattern.c:154
            if (aa1.Ligands[l].AtomSymbol != ligand.AtomSymbol)
              applyAtomSymbolList(mol, ligand.AtomSymbol, al);
            //            }
          }
          if (al->getFormalCharge() != ligand.Charge)
            al->setFormalCharge(ligand.Charge);
          if (al->getNumRadicalElectrons() != ligand.Radical)
            al->setNumRadicalElectrons(ligand.Radical);

          if (BT_NONE == ligand.BondType)  // remove bond
            bondsToRemove.push_back(
                std::pair<unsigned, unsigned>(j, neighbours[j][l].AtomIdx));
          else if (aa1.Ligands[l].BondType != ligand.BondType)
            mol.getBondWithIdx(match_b[l])
                ->setBondType(convertBondType(ligand.BondType));
        }
        transformed = true;
      }
    }
  }
  // remove bonds marked to delete // BT_NONE
  if (!bondsToRemove.empty()) {
    for (auto &i : bondsToRemove) mol.removeBond(i.first, i.second);
  }
  return transformed;
}

//-----------------------------------------------------------------------------
// CheckAtoms():
//-----------------------------------------------------------------------------

bool isBondTypeMatch(const RDKit::Bond &b, AABondType lbt) {
  RDKit::Bond::BondType bt = b.getBondType();
  bool aromatic = b.getIsAromatic();

  if (ANY_BOND == lbt)
    return true;
  else if (lbt == AROMATIC)
    return aromatic || bt == RDKit::Bond::AROMATIC ||
           bt == RDKit::Bond::ONEANDAHALF || bt == RDKit::Bond::TWOANDAHALF ||
           bt == RDKit::Bond::THREEANDAHALF;
  else if (lbt == SINGLE_AROMATIC)
    return (bt == RDKit::Bond::ONEANDAHALF ||
            (aromatic && bt == RDKit::Bond::SINGLE));
  else if (lbt == DOUBLE_AROMATIC)
    return (bt == RDKit::Bond::TWOANDAHALF ||
            (aromatic && bt == RDKit::Bond::DOUBLE));
  else if (lbt == SINGLE)
    return (bt == RDKit::Bond::SINGLE || bt == RDKit::Bond::ONEANDAHALF ||
            bt == RDKit::Bond::AROMATIC || aromatic);
  else if (lbt == DOUBLE)
    return (bt == RDKit::Bond::DOUBLE || bt == RDKit::Bond::TWOANDAHALF ||
            bt == RDKit::Bond::AROMATIC || aromatic);
  else if (lbt == TRIPLE)
    return (bt == RDKit::Bond::TRIPLE || bt == RDKit::Bond::THREEANDAHALF ||
            bt == RDKit::Bond::AROMATIC);
  else
    return false;
}

static bool RecMatchNeigh(std::vector<Neighbourhood> &matched_pairs,
                          std::vector<bool> &visited, std::vector<bool> &used) {
  for (size_t i = 0; i < matched_pairs.size(); i++) {
    if (visited[i]) {
      continue;
    }
    visited[i] = true;
    std::vector<unsigned> &pairs = matched_pairs[i].Atoms;
    for (unsigned int pair : pairs) {
      if (used[pair]) continue;
      used[pair] = true;
      if (RecMatchNeigh(matched_pairs, visited, used)) {
        return true;
      }
      used[pair] = false;
    }
    visited[i] = false;
  }
  for (auto &&j : visited)
    if (!j) return false;
  for (auto &&j : used)
    if (!j) return false;
  return true;
}

/*
 * Recursively searches *mp for the next atom to be appended to
 * match[0..level]. Returns TRUE if all ligands in aap are mapped
 * successfully. nbp[i] describes the neighbour bonds and atoms
 * of atom i.
 */
bool RecMatch(const ROMol &mol, unsigned atomIdx, const AugmentedAtom &aa,
              const std::vector<Neighbourhood> &nbp, bool verbose) {
  const Neighbourhood &nbph = nbp[atomIdx];
  std::vector<Neighbourhood> matched_pairs(nbph.Atoms.size());
  bool found;
  for (unsigned i = 0; i < nbph.Atoms.size(); i++) {
    const Atom &atom = *mol.getAtomWithIdx(nbph.Atoms[i]);
    const Bond &bond = *mol.getBondWithIdx(nbph.Bonds[i]);
    found = false;

    for (unsigned j = 0; j < aa.Ligands.size(); j++) {
      const Ligand &au_ligand = aa.Ligands[j];

      if (au_ligand.Charge != ANY_CHARGE &&
          atom.getFormalCharge() != au_ligand.Charge)
        continue;
      if (au_ligand.Radical != ANY_RADICAL &&
          atom.getNumRadicalElectrons() != au_ligand.Radical)
        continue;
      if (au_ligand.SubstitutionCount != 0 &&
          nbp[nbph.Atoms[i]].Atoms.size() != au_ligand.SubstitutionCount)
        continue;
      if (!isBondTypeMatch(bond, au_ligand.BondType)) continue;
      if (!AtomSymbolMatch(atom.getSymbol(), au_ligand.AtomSymbol)) continue;
      // found matched item:
      found = true;
      matched_pairs[i].Atoms.push_back(j);
    }

    if (!found) {
      return false;
    }
  }
  std::vector<bool> visited(nbph.Atoms.size());
  std::vector<bool> used(aa.Ligands.size());

  for (auto &&j : visited) j = false;
  for (auto &&j : used) j = false;

  if (RecMatchNeigh(matched_pairs, visited, used)) {
    if (verbose)
      BOOST_LOG(rdInfoLog) << "RecMatch ret TRUE " << aa.ShortName << "\n";
    return true;
  }
  //  if (verbose)
  //    BOOST_LOG(rdInfoLog) << "RecMatch ret FALSE " << aa.Ligands.size() << "
  //    "
  //                         << aa.ShortName << "\n";
  return false;
}

/*OLD version (==avalon code)
bool RecMatch(const ROMol &mol, std::vector<unsigned> &match, unsigned int
level,
              const AugmentedAtom &aa,  const std::vector<Neighbourhood> &nbp,
bool verbose) {
    bool is_new;

    if (level == aa.Ligands.size())
        return true;
    const Neighbourhood &nbph = nbp[match[0]];

    for (unsigned i = 0; i < nbph.Atoms.size(); i++) {
        is_new = true;
        for (unsigned j = 1; j < match.size(); j++)
            if (nbph.Atoms[i] == match[j])
                is_new = false;
        if (!is_new)
            continue;

        const Atom &atom = *mol.getAtomWithIdx(nbph.Atoms[i]);
        const Bond &bond = *mol.getBondWithIdx(nbph.Bonds[i]);
        const Ligand& au_ligand = aa.Ligands[level];

        if ((atom.getFormalCharge() == au_ligand.Charge || au_ligand.Charge ==
ANY_CHARGE)
            && (atom.getNumRadicalElectrons() == au_ligand.Radical ||
au_ligand.Radical == ANY_RADICAL)
            && (convertBondType(bond.getBondType()) == au_ligand.BondType
                || (bond.getBondType() == RDKit::Bond::AROMATIC &&
au_ligand.BondType != TRIPLE))
            && (au_ligand.SubstitutionCount == 0 ||
                nbp[nbph.Atoms[i]].Atoms.size() == au_ligand.SubstitutionCount)
            && AtomSymbolMatch(atom.getSymbol(), au_ligand.AtomSymbol)) {

            match.push_back(nbph.Atoms[i]);
            if (RecMatch(mol, match, level + 1, aa, nbp, verbose)) {
                if (verbose)
                    BOOST_LOG(rdInfoLog) << "RecMatch level=" << level << " ret
TRUE\n";
                return true;
            }
        }
    }
    if (verbose)
        BOOST_LOG(rdInfoLog) << "RecMatch level="<< level << " ret FALSE "<<
aa.Ligands.size()<<" "<< aa.ShortName <<"\n";
    return false;
}
*/
/*
 * Tests if atom i in *mp matches the augmented atom description
 * *aap. nbp[] is used to speed up access to neighbour atoms and
 * bonds. The first matching atom mapping is placed into match[1..].
 * i is stored in match[0].
 */
bool AAMatch(const ROMol &mol, unsigned i,
             //    std::vector<unsigned> &match,
             const AugmentedAtom &aa,
             const std::vector<unsigned> &atom_ring_status,
             const std::vector<Neighbourhood> &nbp, bool verbose) {
  const Atom &atom = *mol.getAtomWithIdx(i);

  if (nbp[i].Atoms.size() == aa.Ligands.size() &&
      (aa.Charge == ANY_CHARGE || atom.getFormalCharge() == aa.Charge) &&
      (aa.Radical == ANY_RADICAL ||
       atom.getNumRadicalElectrons() == aa.Radical) &&
      AtomSymbolMatch(atom.getSymbol(), aa.AtomSymbol)) {
    if (!atom_ring_status.empty() && aa.Topology == RING &&
        atom_ring_status[i] == 0) {
      // DEBUG:
      // if (verbose)
      //  BOOST_LOG(rdInfoLog) << "AAMatch i=" << i
      //                       << " aa.Topology RING. ret FALSE\n";
      return false;
    }
    if (!atom_ring_status.empty() && aa.Topology == CHAIN &&
        atom_ring_status[i] != 0) {
      // DEBUG:
      // if (verbose)
      //  BOOST_LOG(rdInfoLog) << "AAMatch i=" << i
      //                       << " aa.Topology CHAIN. ret FALSE\n";
      return false;
    }
    //        if(match.size() == 0) {           match.push_back(i);        }
    //        match[0] = i;
    return RecMatch(mol, i, aa, nbp, verbose);
    //        return RecMatch(mol, match, 0, aa, nbp, verbose);
  } else {
    return false;
  }
}

/*
 * Computes how many basis rings each bond shares and how many
 * ring bonds are attached to an atom. The results are stored in
 * atom_status[] and bond_status[] respectively.
 */
static void RingState(const ROMol &mol, std::vector<unsigned> &atom_status,
                      std::vector<unsigned> &bond_status) {
  atom_status.resize(mol.getNumAtoms());
  bond_status.resize(mol.getNumBonds());
  for (unsigned int &i : atom_status) i = 0;
  for (unsigned int &i : bond_status) i = 0;
  if (mol.getNumBonds() == 0) return;
  // for each bond compute amount of rings that contains the bond
  const RDKit::RingInfo &ringInfo = *mol.getRingInfo();
  const VECT_INT_VECT &bondRings = ringInfo.bondRings();
  for (const auto &ring : bondRings) {
    for (unsigned i = 0; i < ring.size(); i++) {
      bond_status[i]++;
    }
  }
  for (unsigned i = 0; i < bond_status.size(); i++)
    if (bond_status[i] > 0) {
      const Bond &bond = *mol.getBondWithIdx(i);
      atom_status[bond.getBeginAtomIdx()]++;
      atom_status[bond.getEndAtomIdx()]++;
    }
}

/*
 * Checks if every atom in *mp matches one of the augmented atoms
 * in good_atoms[0..ngood-1]. It returns TRUE if all atoms gave a match or FALSE
 * otherwise.
 */
bool CheckAtoms(const ROMol &mol, const std::vector<AugmentedAtom> &good_atoms,
                bool verbose) {
  if (good_atoms.empty()) return true;
  std::vector<Neighbourhood> neighbours(mol.getNumAtoms());
  std::vector<unsigned> atom_status(mol.getNumAtoms());
  std::vector<unsigned> bond_status(mol.getNumBonds());
  //    if(verbose)
  //        BOOST_LOG(rdInfoLog) << "CheckAtoms: " << mol.getNumAtoms() << "
  //        atoms in molecule\n";

  RingState(mol, atom_status, bond_status);
  SetupNeighbourhood(mol, neighbours);

  unsigned nmatch = 0;
  //    int len = good_atoms.size();
  //    for (int xx = 300; xx < 305; xx++) {
  //      BOOST_LOG(rdInfoLog) << xx < "*******\n";
  //      for (int j = 0; j < good_atoms[xx].Ligands.size(); j++) {
  //         BOOST_LOG(rdInfoLog) << good_atoms[xx].Ligands[j].AtomSymbol <<
  //         "\n";
  //      }
  //   }

  for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
    unsigned prevn = nmatch;
    // DEBUG:
    if (verbose) {
      verboseAtom(mol, i, neighbours);
    }
    unsigned j = 0;
    for (; j < good_atoms.size(); j++) {
      // check for ring state of central atom. Ring matches to ring only
      if (good_atoms[j].Topology == RING && 0 == atom_status[i]) {
        // DEBUG:
        if (verbose &&  // j >= 12 && j <= 21 &&
            mol.getAtomWithIdx(i)->getAtomicNum() == 6)  // 'C'
          BOOST_LOG(rdInfoLog) << "UNMATCHED ring state RING of atom idx=" << i
                               << " " << mol.getAtomWithIdx(i)->getSymbol()
                               << " status=" << atom_status[i] << "\n";
        continue;
      }
      if (good_atoms[j].Topology == CHAIN && 0 != atom_status[i]) {
        // DEBUG:
        if (verbose &&  // j >= 12 && j <= 21 &&
            mol.getAtomWithIdx(i)->getAtomicNum() == 6)  // 'C'
          BOOST_LOG(rdInfoLog) << "UNMATCHED ring state CHAIN of atom idx=" << i
                               << " " << mol.getAtomWithIdx(i)->getSymbol()
                               << " status=" << atom_status[i] << "\n";
        continue;
      }

      // DEBUG:
      //            if (verbose && j >= 12 && j <= 21 &&
      //            mol.getAtomWithIdx(i)->getAtomicNum() == 6) // 'C'
      //                BOOST_LOG(rdInfoLog) << "AAMatch i="<<i<<" j="<<j;

      std::vector<unsigned> match;  // unused [MAXNEIGHBOURS + 1];
      if (neighbours[i].Atoms.size() == good_atoms[j].Ligands.size() &&
          AAMatch(mol, i, good_atoms[j], atom_status, neighbours, verbose)) {
        if (verbose) {
          BOOST_LOG(rdInfoLog)
              << "AAMatch idx=" << i << " aa_idx=" << j << " ret TRUE !\n";
        }
        // DEBUG:
        //        if (verbose && j >= 12 && j <= 21 &&
        //            mol.getAtomWithIdx(i)->getAtomicNum() == 6)  // 'C'
        //          BOOST_LOG(rdInfoLog) << "AAMatch i=" << i << " j=" << j
        //                               << " ret TRUE !\n";
        nmatch++;
        break;
      }
      // DEBUG:
      //      else if (verbose && j >= 12 && j <= 21 &&
      //               mol.getAtomWithIdx(i)->getAtomicNum() == 6)  // 'C'
      //        BOOST_LOG(rdInfoLog) << "AAMatch i=" << i << " j=" << j
      //                             << " ret FALSE\n";
    }
    if (verbose && nmatch == prevn)  // UNMATCHED atom
      BOOST_LOG(rdInfoLog) << "UNMATCHED atom idx=" << i << " "
                           << mol.getAtomWithIdx(i)->getSymbol()
                           << " status=" << atom_status[i]
                           << " nmatch=" << nmatch << std::endl
                           << LogNeighbourhood(mol, i, neighbours) << std::endl;

    if (verbose && j == good_atoms.size()) {  // failed this atom .. log it
      BOOST_LOG(rdWarningLog)
          << LogNeighbourhood(mol, i, neighbours) << std::endl;
    }
  }
  return (nmatch == mol.getNumAtoms());
}

}  // namespace StructureCheck
}  // namespace RDKit
