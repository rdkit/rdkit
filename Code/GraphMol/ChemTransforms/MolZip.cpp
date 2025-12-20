//
//  Copyright (C) 2020-2025 Brian Kelley and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MolFragmenter.h"

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "ChemTransforms.h"
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <Geometry/Transform3D.h>

#include <cstdint>
#include <map>
#include <optional>
#include <vector>

namespace RDKit {
namespace {
const unsigned int NOLABEL = std::numeric_limits<unsigned int>::max();
// Get the atom label - this might be useful as a util class
unsigned int get_label(const Atom *a, const MolzipParams &p) {
  PRECONDITION(a, "bad atom in MolZip::get_label")
  unsigned int idx = NOLABEL;
  switch (p.label) {
    case MolzipLabel::AtomMapNumber:
      if (a->getAtomicNum() == 0) {
        auto mapno = a->getAtomMapNum();
        return mapno ? mapno : NOLABEL;
      }
      break;

    case MolzipLabel::Isotope:
      if (a->getAtomicNum() == 0) {
        auto iso = a->getIsotope();
        return iso ? iso : NOLABEL;
      }
      break;

    case MolzipLabel::AtomType:
      idx = std::distance(p.atomSymbols.begin(),
                          std::find(p.atomSymbols.begin(), p.atomSymbols.end(),
                                    a->getSymbol()));
      if (idx == p.atomSymbols.size()) {
        idx = NOLABEL;
      }
      break;

    case MolzipLabel::FragmentOnBonds:
      // shouldn't ever get here
      CHECK_INVARIANT(
          0, "FragmentOnBonds is not an atom label, it is an atom index");
      break;

    case MolzipLabel::AtomProperty:
      a->getPropIfPresent<unsigned int>(p.atomProperty, idx);
      break;

    default:
      CHECK_INVARIANT(0, "bogus MolZipLabel value in MolZip::get_label");
  }
  return idx;
}

// Return the connected atom
//  n.b. There can be only one connection from a mapped atom
Atom *get_other_atom(Atom *a) {
  PRECONDITION(a, "null atom in MolZip::get_other_atom");
  auto &m = a->getOwningMol();
  if (m.getAtomDegree(a) != 1) {
    return nullptr;
  }

  return m[*m.getAtomNeighbors(a).first];
}

int num_swaps_to_interconvert(std::vector<unsigned int> &orders) {
  int nswaps = 0;
  std::vector<bool> seen(orders.size());
  for (size_t i = 0; i < orders.size(); ++i) {
    if (!seen[i]) {
      auto j = i;
      while (orders[j] != i) {
        j = orders[j];
        CHECK_INVARIANT(
            j < orders.size(),
            "molzip: bond index outside of number of bonds for atom")
        seen[j] = true;
        nswaps++;
      }
    }
  }
  return nswaps;
}

// Simple bookkeeping class to bond attachments and handle stereo
struct ZipBond {
  Atom *a = nullptr;  // atom being bonded
  Atom *a_dummy =
      nullptr;  // Labelled atom, i.e. [*:1]-C  will bond the C to something
  Atom *b = nullptr;  // atom being bonded
  Atom *b_dummy =
      nullptr;  // Labelled atom, i.e. [*:1]-O will bond the O to something
  Atom *a_link =
      nullptr;  // Link bonds have six atoms, [*:a][*b].[*:a]C.[*:b]D, four
                // dummies and two atoms
  Atom *b_link = nullptr;
  bool isLinker = false;          // is this a straight linker bond
  Bond::BondType linkerBondType;  // The linker bond type

  // Backup the original chirality mark_chirality must be called first
  //  as it checks the datastructure for validity;
  void mark_chirality() const {
    PRECONDITION(a, "Must have a begin atom to bond");
    PRECONDITION(b, "Must have an end atom to bond");
    PRECONDITION(a_dummy, "Must have a begin dummy atom");
    PRECONDITION(b_dummy, "Must have an end dummy atom");

    mark(a, a_dummy, b);
    mark(b, b_dummy, a);
  }

  // bond a<->b for now only use single bonds
  //  XXX FIX ME take the highest bond order.
  // RETURNS true if a new bond was made
  bool bond(RWMol &newmol, const MolzipParams &params) const {
    if (!a || !b || !a_dummy || !b_dummy) {
      BOOST_LOG(rdWarningLog)
          << "Incomplete atom labelling, cannot make bond" << std::endl;
      return false;
    }

    // Fragment on bonds allows multiple links to the same atom
    // i.e. C.[1C].[1C]
    //  otherwise throw an invariant error

    CHECK_INVARIANT(
        params.label == MolzipLabel::FragmentOnBonds ||
            !a->getOwningMol().getBondBetweenAtoms(a->getIdx(), b->getIdx()),
        "molzip: zipped Bond already exists, perhaps labels are duplicated");
    if (!a->getOwningMol().getBondBetweenAtoms(a->getIdx(), b->getIdx())) {
      CHECK_INVARIANT(&a->getOwningMol() == &newmol,
                      "Owning mol is not the combined molecule!!");
      if (isLinker) {
        // This is the easy bit, just link a and b, and schedule a_dummy and
        // b_dummy
        //  for deletion
        CHECK_INVARIANT(
            a && b && a_dummy && b_dummy && a_link && b_link,
            "molzip: Link Bond is missing one or more labelled atoms");
        newmol.addBond(a, b, linkerBondType);
        a_link->setProp("__molzip_used", true);
        b_link->setProp("__molzip_used", true);
      } else {
        auto bnd = newmol.getBondBetweenAtoms(a->getIdx(), a_dummy->getIdx());
        CHECK_INVARIANT(
            bnd != nullptr,
            "molzip: begin atom and specified dummy atom connection "
            "are not bonded.")
        auto bond_type_a = bnd->getBondType();
        auto bond_dir_a = bnd->getBondDir();
        auto a_is_start = bnd->getBeginAtom() == a;

        bnd = newmol.getBondBetweenAtoms(b->getIdx(), b_dummy->getIdx());
        CHECK_INVARIANT(bnd != nullptr,
                        "molzip: end atom and specified dummy connection atom "
                        "are not bonded.")
        auto bond_type_b = bnd->getBondType();

        auto bond_dir_b = bnd->getBondDir();
        auto b_is_start = bnd->getBeginAtom() == b;

        unsigned int bnd_idx = 0;
        // Fusion bond-dir logic table
        // a-* b-* => a-b
        //  < = wedge

        //  a<* b-* => a<b
        //  a>* b-* => a>b
        //  a-* b>* => a<b
        //  a-* b<* => a>b
        Bond::BondDir bond_dir{Bond::BondDir::NONE};
        auto start = a;
        auto end = b;
        if (bond_dir_a != Bond::BondDir::NONE &&
            bond_dir_b != Bond::BondDir::NONE) {
          // are we consistent between the two bond orders check for the case of
          // fragment on bonds where a<* and b>* or a>* and b<* when < is either
          // a hash or wedge bond but not both.
          bool consistent_directions = false;
          if (bond_dir_a == bond_dir_b) {
            if ((a_is_start != b_is_start)) {
              consistent_directions = true;
            }
          }
          if (!consistent_directions) {
            BOOST_LOG(rdWarningLog)
                << "inconsistent bond directions when merging fragments, ignoring..."
                << std::endl;
            bond_dir_a = bond_dir_b = Bond::BondDir::NONE;
          } else {
            bond_dir_b = Bond::BondDir::NONE;
          }
        }

        if (bond_dir_a != Bond::BondDir::NONE) {
          if (!a_is_start) {
            start = b;
            end = a;
          }
          bond_dir = bond_dir_a;
        } else if (bond_dir_b != Bond::BondDir::NONE) {
          if (b_is_start) {
            start = b;
            end = a;
          }
          bond_dir = bond_dir_b;
        }

        if (bond_type_a != Bond::BondType::SINGLE) {
          bnd_idx = newmol.addBond(start, end, bond_type_a);
        } else if (bond_type_b != Bond::BondType::SINGLE) {
          bnd_idx = newmol.addBond(start, end, bond_type_b);
        } else {
          bnd_idx = newmol.addBond(start, end, Bond::BondType::SINGLE);
        }

        newmol.getBondWithIdx(bnd_idx - 1)->setBondDir(bond_dir);
      }
    }
    a_dummy->setProp("__molzip_used", true);
    b_dummy->setProp("__molzip_used", true);

    return true;
  }

  // Restore the marked chirality (mark_chirality must be called first)
  void restore_chirality(std::set<Atom *> &already_checked) const {
    PRECONDITION(a, "Must have a begin atom to bond");
    PRECONDITION(b, "Must have an end atom to bond");
    PRECONDITION(a_dummy, "Must have a begin dummy atom");
    PRECONDITION(b_dummy, "Must have an end dummy atom");
    if (already_checked.find(a) == already_checked.end()) {
      restore(a);
      already_checked.insert(a);
    }
    if (already_checked.find(b) == already_checked.end()) {
      restore(b);
      already_checked.insert(b);
    }

    // now do bond stereo
    std::string mark = "__molzip_bond_stereo_mark";
    for (auto *bond : a->getOwningMol().bonds()) {
      if (bond->hasProp(mark)) {
        std::vector<int> atoms;
        for (auto *atom : bond->getProp<std::vector<Atom *>>(mark)) {
          atoms.push_back(rdcast<int>(atom->getIdx()));
        }
        bond->getStereoAtoms().swap(atoms);
        bond->setStereo(
            bond->getProp<Bond::BondStereo>("__molzip_bond_stereo"));
      }
    }
  }

 private:
  // Mark the original order of the nbr atoms including the dummy
  //  The goal is to copy the dummy chiral order over to the
  //  atom being bonded
  void mark(Atom *chiral_atom, Atom *dummy_atom, Atom *new_atom) const {
    if (chiral_atom->getChiralTag()) {
      std::string mark =
          "__molzip_mark_" + std::to_string(chiral_atom->getIdx());
      chiral_atom->setProp("__molzip_chiral_mark", mark);
      int order = 0;
      auto &m = chiral_atom->getOwningMol();
      for (auto nbrIdx :
           boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
        m[nbrIdx]->setProp(mark, order);
        ++order;
      }
      new_atom->setProp(mark, dummy_atom->getProp<int>(mark));
    }

    // check bond stereo
    auto &m = chiral_atom->getOwningMol();
    for (auto bond : m.atomBonds(chiral_atom)) {
      if (bond->getStereo() != Bond::BondStereo::STEREONONE) {
        std::string mark = "__molzip_bond_stereo_mark";
        std::vector<Atom *> atoms;
        bool has_dummy = false;
        for (auto idx : bond->getStereoAtoms()) {
          if (static_cast<unsigned>(idx) == dummy_atom->getIdx()) {
            atoms.push_back(new_atom);
            has_dummy = true;
          } else {
            atoms.push_back(m.getAtomWithIdx(idx));
          }
        }
        if (has_dummy) {
          bond->setProp(mark, atoms);
          bond->setProp<Bond::BondStereo>("__molzip_bond_stereo",
                                          bond->getStereo());
        }
      }
    }
  }

  // Restore the atom's chirality by comparing the original order
  //  to the current
  void restore(Atom *chiral_atom) const {
    if (!chiral_atom->getChiralTag()) {
      return;
    }
    std::string mark =
        chiral_atom->getProp<std::string>("__molzip_chiral_mark");
    // std::vector<unsigned int> orders1;
    std::vector<unsigned int> orders2;
    auto &m = chiral_atom->getOwningMol();
    for (auto nbrIdx :
         boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
      orders2.push_back(m[nbrIdx]->getProp<int>(mark));
    }
    if (num_swaps_to_interconvert(orders2) % 2 == 1) {
      chiral_atom->invertChirality();
    }
  }
};

void rotateFragmentToBondVector(
    ROMol &mol, const Atom &a, const Atom &b, const Atom &a_dummy,
    const Atom &b_dummy, const std::vector<int> &fragmentForAtom,
    const std::map<int, std::vector<int>> &atomsInFragment, int confId = -1) {
  if (!mol.getNumConformers()) {
    return;
  }
  auto &conf = mol.getConformer(confId);

  //-----------------------------------------
  //  Notation:
  //    Pmc: molecule connection point (the atom that will be
  //     removed from the molecule).
  //    Pma: molecule attachment point (the atom to which we'll form
  //     the bond).
  //    Psc: sidechain connection point
  //    Psa: sidechain attachment point
  //    Vm: Pmc-Pma (molecular attachment vector)
  //    Vs: Psc-Psa (sidechain attachment vector)
  //
  //-----------------------------------------
  const auto Pma = conf.getAtomPos(a.getIdx());
  const auto Psa = conf.getAtomPos(b.getIdx());
  const auto Pmc = conf.getAtomPos(a_dummy.getIdx());
  const auto Psc = conf.getAtomPos(b_dummy.getIdx());

  auto Um = Pma.directionVector(Pmc);
  // note the opposite direction here:
  auto Us = Psc.directionVector(Psa);

  RDGeom::Transform3D templateTform;
  templateTform.SetTranslation(Pma);

  double sinT, cosT;
  cosT = Us.dotProduct(Um);
  if (cosT > 1.0) cosT = 1.0;
  if (fabs(cosT) < 1.0) {
    sinT = sqrt(1.0 - cosT * cosT);
    RDGeom::Point3D rotnAxis = Us.crossProduct(Um);
    rotnAxis.normalize();
    templateTform.SetRotation(cosT, sinT, rotnAxis);
  } else if (cosT == 1.0) {
    RDGeom::Point3D normal(1, 0, 0);
    if (fabs(Us.dotProduct(normal)) == 1.0) {
      normal = RDGeom::Point3D(0, 1, 0);
    }
    RDGeom::Point3D rotnAxis = Us.crossProduct(normal);
    templateTform.SetRotation(-1, 0, rotnAxis);
  }

  // we use the second attachment vector to set the bond distance
  RDGeom::Transform3D tmpTform;
  tmpTform.SetTranslation(Psc * -1.0);
  templateTform *= tmpTform;

  // ---------
  // transform the atomic positions in the sidechain:
  // ---------
  auto fragId = fragmentForAtom[b.getIdx()];
  for (const auto atomIdx : atomsInFragment.at(fragId)) {
    RDGeom::Point3D pos = conf.getAtomPos(atomIdx);
    templateTform.TransformPoint(pos);
    conf.setAtomPos(atomIdx, pos);
  }
}

}  // namespace

static const std::string indexPropName("__zipIndex");

std::unique_ptr<ROMol> molzip(
    const ROMol &a, const ROMol &b, const MolzipParams &params,
    std::optional<std::map<int, int>> &attachmentMapping) {
  if (attachmentMapping) {
    attachmentMapping->clear();
  }
  std::unique_ptr<RWMol> newmol;
  if (b.getNumAtoms()) {
    newmol.reset(static_cast<RWMol *>(combineMols(a, b)));
  } else {
    newmol.reset(new RWMol(a));
  }

  // doing the coordinate alignment is quicker if we know which atoms are in
  // which fragments
  std::vector<int> fragmentForAtom;
  std::map<int, std::vector<int>> atomsInFragment;
  if (params.alignCoordinates) {
    MolOps::getMolFrags(*newmol, fragmentForAtom);
    for (size_t i = 0; i < fragmentForAtom.size(); ++i) {
      atomsInFragment[fragmentForAtom[i]].push_back(static_cast<int>(i));
    }
  }

  std::map<unsigned int, ZipBond> mappings;
  std::map<Atom *, std::vector<const ZipBond *>> mappings_by_atom;

  // Linker bonds resolve to the same zip bond by using the
  //  lowest label.  I.e.
  //   [*:1][*:2] sets the link bond to label 1
  //    so this sets linkerBonds[1] == linkerBonds[2] = the same ZipBond
  std::map<unsigned int, ZipBond *> linkerBonds;

  std::vector<Atom *> deletions;
  if (params.label == MolzipLabel::FragmentOnBonds) {
    for (auto *atom : newmol->atoms()) {
      if (atom->getAtomicNum() == 0) {
        auto molno = atom->getIsotope();
        auto attached_atom = get_other_atom(atom);
        auto &bond = mappings[molno];
        bond.a = attached_atom;
        bond.a_dummy = atom;
        bond.b = newmol->getAtomWithIdx(molno);
        for (auto nbrIdx :
             boost::make_iterator_range(newmol->getAtomNeighbors(bond.b))) {
          auto *nbr = (*newmol)[nbrIdx];
          if (nbr->getAtomicNum() == 0 &&
              nbr->getIsotope() == attached_atom->getIdx()) {
            bond.b_dummy = nbr;
            break;
          }
        }
        if (!bond.b_dummy) {
          BOOST_LOG(rdErrorLog)
              << "Cannot find atom to bond using FragmentOnBond labelling"
              << std::endl;
          return std::unique_ptr<ROMol>();
        }
        mappings_by_atom[atom].push_back(&bond);
        deletions.push_back(atom);
        if (attachmentMapping) {
          if (int otherIndex, dummyIndex;
              atom->getPropIfPresent(indexPropName, dummyIndex) &&
              bond.b->getPropIfPresent(indexPropName, otherIndex)) {
            (*attachmentMapping)[dummyIndex] = otherIndex;
          }
        }
      }
    }
  } else {  // Non Fragment By Bonds attaching
    for (auto *atom : newmol->atoms()) {
      auto molno = get_label(atom, params);
      if (molno != NOLABEL) {
        auto attached_atom = get_other_atom(atom);
        auto attached_molno =
            attached_atom ? get_label(attached_atom, params) : NOLABEL;
        if (attached_molno != NOLABEL) {
          // we have a linker bond
          //  [*:1][*:2].[*:1]C.[*:2]S links C and S and drops all dummies
          //   Note:  the linker bond MUST come first here
          // Get the min molno and use this to assign the bonds to link
          if (molno > attached_molno) {
            std::swap(molno, attached_molno);
            std::swap(atom, attached_atom);
          }

          auto link_bond = atom->getOwningMol().getBondBetweenAtoms(
              atom->getIdx(), attached_atom->getIdx());
          CHECK_INVARIANT(
              link_bond,
              ("molzip: link bond with labels: " + std::to_string(molno) + "," +
               std::to_string(attached_molno) + " is missing"));
          auto bondType = link_bond->getBondType();
          if (mappings.find(molno) == mappings.end()) {
            auto &bond = mappings[molno];
            CHECK_INVARIANT(
                linkerBonds.find(attached_molno) == linkerBonds.end(),
                ("molzip: Linker attachment point with label: " +
                 std::to_string(attached_molno) + " found before linker bond"));
            linkerBonds[attached_molno] = &bond;
            bond.isLinker = true;
            bond.linkerBondType = bondType;
            bond.a_link = atom;
            bond.b_link = attached_atom;
            deletions.push_back(atom);
            deletions.push_back(attached_atom);
          } else {
            // we'll find this bond twice, so let's make sure it is setup
            // correctly
            auto &bond = mappings[molno];
            CHECK_INVARIANT(
                bondType = bond.linkerBondType,
                ("molzip: LINKER bond with labels: " + std::to_string(molno) +
                 "," + std::to_string(attached_molno) +
                 " has inconsistent bond types"));
            CHECK_INVARIANT(
                bond.isLinker,
                ("molzip: LINKER bond with labels: " + std::to_string(molno) +
                 "," + std::to_string(attached_molno) +
                 " found but not encountered first in the molecules to be zipped."));
            CHECK_INVARIANT(
                (bond.a_link == atom && bond.b_link == attached_atom) ||
                    (bond.b_link == atom && bond.a_link == attached_atom),
                ("molzip: Linker Bond with labels " + std::to_string(molno) +
                 "," + std::to_string(attached_molno) +
                 " not setup correctly"));
          }
        } else if (mappings.find(molno) == mappings.end() &&
                   linkerBonds.find(molno) == linkerBonds.end()) {
          // Normal linkage C[*:1].S[*:1] links C and S
          // LinkBond [*:1][*:2].C[*:1].S[*:2] links C and S
          auto &bond = mappings[molno];
          CHECK_INVARIANT(
              !bond.a,
              "molzip: bond info already setup for bgn atom with label:" +
                  std::to_string(molno));
          bond.a = attached_atom;
          bond.a_dummy = atom;
        } else {
          auto &bond = linkerBonds.find(molno) == linkerBonds.end()
                           ? mappings[molno]
                           : *linkerBonds[molno];
          if (bond.isLinker) {
            CHECK_INVARIANT(
                !bond.a || !bond.b,
                "molzip: Linker bond has multiple attachments for label: " +
                    std::to_string(molno));
            if (bond.a) {
              bond.b = attached_atom;
              bond.b_dummy = atom;
              deletions.push_back(bond.b_dummy);
            } else {
              bond.a = attached_atom;
              bond.a_dummy = atom;
              deletions.push_back(bond.a_dummy);
            }
          } else {
            CHECK_INVARIANT(
                bond.a,
                "molzip: bond info not properly setup for bgn atom with label:" +
                    std::to_string(molno));
            CHECK_INVARIANT(
                !bond.b,
                "molzip: bond info already exists for end atom with label:" +
                    std::to_string(molno));

            bond.b = attached_atom;
            bond.b_dummy = atom;
          }

          mappings_by_atom[bond.a].push_back(&bond);
          if (attachmentMapping) {
            if (int otherIndex, dummyIndex;
                bond.a_dummy->getPropIfPresent(indexPropName, dummyIndex) &&
                bond.b->getPropIfPresent(indexPropName, otherIndex)) {
              (*attachmentMapping)[dummyIndex] = otherIndex;
            }
            if (int otherIndex, dummyIndex;
                bond.b_dummy->getPropIfPresent(indexPropName, dummyIndex) &&
                bond.a->getPropIfPresent(indexPropName, otherIndex)) {
              (*attachmentMapping)[dummyIndex] = otherIndex;
            }
          }
        }
        deletions.push_back(atom);
      }
    }
  }

  // Mark the existing chirality so we can try and restore it later
  for (auto &kv : mappings_by_atom) {
    for (auto &bond : kv.second) {
      bond->mark_chirality();
    }
  }

  // Make all the bonds
  boost::dynamic_bitset<> aligned(atomsInFragment.size());
  for (auto &kv : mappings) {
    if (kv.second.bond(*newmol, params) && params.alignCoordinates) {
      // handle rotating coordinates here if we made the bond
      if (params.alignCoordinates) {
        // make sure we haven't aligned this fragment already:
        if (!aligned[fragmentForAtom[kv.second.b->getIdx()]]) {
          std::cerr << "rotate: " << kv.second.a->getIdx() << "-"
                    << kv.second.b->getIdx() << " "
                    << kv.second.a_dummy->getIdx() << "-"
                    << kv.second.b_dummy->getIdx() << std::endl;
          rotateFragmentToBondVector(*newmol, *(kv.second.a), *(kv.second.b),
                                     *(kv.second.a_dummy), *(kv.second.b_dummy),
                                     fragmentForAtom, atomsInFragment);
          aligned.set(fragmentForAtom[kv.second.b->getIdx()]);
        }
      }
    }
    // handle rotating coordinates here if we made the bond
    // be sure to detect cases where the same fragment ends up attaching
    // multiple times. Can probably do that using the __molzip_used tag
  }
  newmol->beginBatchEdit();

  // Remove the used atoms
  for (auto &atom : deletions) {
    if (atom->hasProp("__molzip_used")) {
      newmol->removeAtom(atom);
    }
  }
  newmol->commitBatchEdit();

  // Try and restore the chirality now that we have new bonds
  std::set<Atom *> already_checked;
  for (auto &kv : mappings_by_atom) {
    for (auto &bond : kv.second) {
      bond->restore_chirality(already_checked);
    }
  }

  // remove all molzip tags
  for (auto *atom : newmol->atoms()) {
    auto propnames = atom->getPropList();
    for (auto &prop : propnames) {
      if (prop.find("__molzip") == 0) {
        atom->clearProp(prop);
      }
    }
  }
  for (auto *bond : newmol->bonds()) {
    auto propnames = bond->getPropList();
    for (auto &prop : propnames) {
      if (prop.find("__molzip") == 0) {
        bond->clearProp(prop);
      }
    }
  }
  newmol->updatePropertyCache(params.enforceValenceRules);
  newmol->setProp(common_properties::_StereochemDone, true);
  return newmol;
}

RDKIT_CHEMTRANSFORMS_EXPORT std::unique_ptr<ROMol> molzip(
    const ROMol &a, const ROMol &b, const MolzipParams &params) {
  std::optional<std::map<int, int>> opt(std::nullopt);
  return molzip(a, b, params, opt);
}

std::unique_ptr<ROMol> molzip(const ROMol &a, const MolzipParams &params) {
  const static ROMol b;
  return molzip(a, b, params);
}

std::unique_ptr<ROMol> molzip(std::vector<ROMOL_SPTR> &decomposition,
                              const MolzipParams &params) {
  if (params.generateCoordinates) {
    int index = 0;
    for (const auto &mol : decomposition) {
      for (const auto atom : mol->atoms()) {
        atom->setProp(indexPropName, index++);
      }
    }
  }

  if (decomposition.empty()) {
    return nullptr;
  }

  // When the rgroup decomposition splits a ring, it puts it in both
  //  rgroups, so remove these
  std::vector<ROMOL_SPTR> mols;
  if (params.label != MolzipLabel::FragmentOnBonds &&
      decomposition.size() > 1) {
    std::vector<std::string> existing_smiles;
    for (size_t idx = 1; idx < decomposition.size(); ++idx) {
      auto &mol = decomposition[idx];
      auto smiles = MolToSmiles(*mol);
      if (std::find(existing_smiles.begin(), existing_smiles.end(), smiles) ==
          existing_smiles.end()) {
        mols.push_back(mol);
        existing_smiles.push_back(smiles);
      }
    }
  }

  auto combinedMol = decomposition[0];
  if (!mols.empty()) {
    combinedMol = std::accumulate(
        mols.begin(), mols.end(), decomposition[0],
        [](const auto &combined, const auto &mol) {
          return boost::shared_ptr<ROMol>(combineMols(*combined, *mol));
        });
  }
  const static ROMol b;
  std::optional attachmentMappingOption = std::map<int, int>();
  auto zippedMol = molzip(*combinedMol, b, params, attachmentMappingOption);

  if (params.generateCoordinates && zippedMol->getNumAtoms() > 0) {
    const auto confId = RDDepict::compute2DCoords(*zippedMol);
    const auto zippedConf = zippedMol->getConformer(confId);
    auto attachmentMapping = *attachmentMappingOption;
    for (auto &mol : decomposition) {
      const auto newConf = new Conformer(mol->getNumAtoms());
      newConf->set3D(false);
      for (const auto atom : mol->atoms()) {
        int zippedIndex = atom->getProp<int>(indexPropName);
        atom->clearProp(indexPropName);
        if (const auto attachment = attachmentMapping.find(zippedIndex);
            attachment != attachmentMapping.end()) {
          zippedIndex = (*attachment).second;
        }
        auto zipppedAtoms = zippedMol->atoms();
        auto zippedAtom = std::find_if(
            zipppedAtoms.begin(), zipppedAtoms.end(),
            [zippedIndex](const Atom *zippedAtom) {
              const auto index = zippedAtom->getProp<int>(indexPropName);
              return index == zippedIndex;
            });

        newConf->setAtomPos(atom->getIdx(),
                            zippedConf.getAtomPos((*zippedAtom)->getIdx()));
      }
      mol->addConformer(newConf, true);
    }
    for (const auto atom : zippedMol->atoms()) {
      atom->clearProp(indexPropName);
    }
  }

  return zippedMol;
}

std::unique_ptr<ROMol> molzip(const std::map<std::string, ROMOL_SPTR> &row,
                              const MolzipParams &params) {
  auto core = row.find("Core");
  PRECONDITION(core != row.end(), "RGroup has no Core, cannot molzip");
  std::vector<ROMOL_SPTR> mols;
  mols.push_back(core->second);
  for (auto it : row) {
    if (it.first != "Core") {
      mols.push_back(it.second);
    }
  }

  return molzip(mols, params);
}

}  // end of namespace RDKit
