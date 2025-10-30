#include "utils.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileWriters.h>
namespace RDKit {
namespace ChemDraw {
std::string NodeType(CDXNodeType nodetype) {
  switch (nodetype) {
    case kCDXNodeType_Unspecified:
      return "Unspecified";
    case kCDXNodeType_Element:
      return "Element";
    case kCDXNodeType_ElementList:
      return "ElementList";
    case kCDXNodeType_ElementListNickname:
      return "ElementListNickname";
    case kCDXNodeType_Nickname:
      return "Nickname";
    case kCDXNodeType_Fragment:
      return "Fragment";
    case kCDXNodeType_Formula:
      return "Formula";
    case kCDXNodeType_GenericNickname:
      return "GenericNickname";
    case kCDXNodeType_AnonymousAlternativeGroup:
      return "Anonymous Alternative Group";
    case kCDXNodeType_NamedAlternativeGroup:
      return "Named Alternative Group";
    case kCDXNodeType_MultiAttachment:
      return "MultiAttachment";
    case kCDXNodeType_VariableAttachment:
      return "Variable Attachment";
    case kCDXNodeType_ExternalConnectionPoint:
      return "ExternalConnectionPoint";
    case kCDXNodeType_LinkNode:
      return "LinkNode";
    case kCDXNodeType_Monomer:
      return "Monomer";
    default:
      return "?";
  }
}

void scaleBonds(const ROMol &mol, Conformer &conf, double targetBondLength,
                double bondLength) {
  double avg_bond_length = 0.0;
  if (bondLength < 0) {
    // If we don't have a bond length for any reason, just scale the average
    // bond length
    for (auto &bond : mol.bonds()) {
      avg_bond_length += (conf.getAtomPos(bond->getBeginAtomIdx()) -
                          conf.getAtomPos(bond->getEndAtomIdx()))
                             .length();
    }
    avg_bond_length /= mol.getNumBonds();
  } else {
    avg_bond_length = bondLength;
  }

  if (avg_bond_length > 0) {
    double scale = targetBondLength / avg_bond_length;
    for (auto &pos : conf.getPositions()) {
      pos *= scale;
    }
  }
}

unsigned int get_fuse_label(Atom *atm) {
  // return atm->getAtomMapNum(); easier debugging
  unsigned int label = 0;  // default is no label
  atm->getPropIfPresent<unsigned int>(FUSE_LABEL, label);
  return label;
}

void set_fuse_label(Atom *atm, unsigned int idx) {
  // atm->setAtomMapNum(idx); //for debugging
  if (idx) {
    atm->setProp<unsigned int>(FUSE_LABEL, idx);
  } else {
    atm->clearProp(FUSE_LABEL);
  }
}

struct FragmentReplacement {
  // R = Replacement
  // F = Fragment
  // C = Connection
  //                    C R C F     F
  //                    N=*=C.*=CCC=*
  //  label               1   1     1
  //  has bond ordering
  //
  //  goal replace the atom R with the connections
  unsigned int label = 0;
  Atom *replacement_atom = nullptr;

  std::vector<Atom *> replacement_connection_atoms;
  std::vector<Atom *> fragment_atoms;

  bool replace(RWMol &mol) {
    if (!replacement_atom) {
      return true;
    }

    auto bond_ordering =
        replacement_atom->getProp<std::vector<int>>(CDX_BOND_ORDERING);

    // The "addBond" lower in the loop modifies the atomBonds iterator, which
    // made this loop re-run the 1st iteration. So we copy the bonds first to
    // ensure we don't modify while iterating.
    std::vector<Bond *> replacement_bonds;
    for (auto bond : mol.atomBonds(replacement_atom)) {
      replacement_bonds.push_back(bond);
    }

    // Find the connecting atoms and and do the replacement
    for (auto bond : replacement_bonds) {
      // find the position of the attachment bonds in the bond ordering
      unsigned bond_id = 0;
      if (!bond->getPropIfPresent<unsigned int>(CDX_BOND_ID, bond_id)) {
        BOOST_LOG(rdWarningLog)
            << "bond missing internal CDX BOND id, can't attach fragment at bond:"
            << std::endl;
        return false;
      }
      auto it = std::find(bond_ordering.begin(), bond_ordering.end(), bond_id);
      if (it == bond_ordering.end()) {
        return false;
      }

      auto pos = std::distance(bond_ordering.begin(), it);
      if (pos < 0 || (size_t)pos >= fragment_atoms.size()) {
        BOOST_LOG(rdWarningLog)
            << "bond ordering and number of atoms in fragment mismatch, can't attach fragment at bond:"
            << bond_id << std::endl;

        return false;
      }

      // The "addBond" lower in the loop modifies the atomBonds iterator, which
      // gave trouble with the replacement bonds. To ensure safety, we copy the bonds first.
      auto &xatom = fragment_atoms[pos];
      std::vector<Bond *> xbonds;
      for (auto &xbond : mol.atomBonds(xatom)) {
        xbonds.push_back(xbond);
      }
      for (auto &xbond : xbonds) {
        // xatom is the fragment dummy atom
        // xbond is the fragment bond
        if (bond->getBeginAtom() == replacement_atom) {
          // bond from fragement to after replacement atom
          mol.addBond(xbond->getOtherAtom(xatom), bond->getEndAtom(),
                      bond->getBondType());
        } else {
          // bond from before replacement atom to fragment
          mol.addBond(bond->getBeginAtom(), xbond->getOtherAtom(xatom),
                      bond->getBondType());
        }
      }
    }

    mol.removeAtom(replacement_atom);
    for (auto &atom : fragment_atoms) {
      mol.removeAtom(atom);
    }
    return true;
  }
};

// Replace fragments that are not possible with molzip
bool replaceFragments(RWMol &mol) {
  // Anything with a single atom that is supposed to be replaced via a fragment
  // is here
  std::map<int, FragmentReplacement> replacements;

  for (auto &atom : mol.atoms()) {
    auto label = get_fuse_label(atom);
    if (label) {
      if (atom->hasProp(CDX_BOND_ORDERING)) {
        auto &frag = replacements[label];
        frag.label = label;
        frag.replacement_atom = atom;
      } else {
        // The is the fragment attachment atoms that need to
        //  be attached to the ones connected to the atom being replaced
        auto &frag = replacements[label];
        frag.fragment_atoms.push_back(atom);
      }
    }
  }
  mol.beginBatchEdit();
  for (auto &replacement : replacements) {
    replacement.second.replace(mol);
  }
  mol.commitBatchEdit();
  return true;
}
namespace {
Atom::ChiralType getChirality(ROMol &mol, Atom *center_atom, Conformer &conf) {
  if (center_atom->hasProp(CDX_BOND_ORDERING)) {
    auto bond_ordering =
        center_atom->getProp<std::vector<int>>(CDX_BOND_ORDERING);
    if (bond_ordering.size() < 3) {
      return Atom::ChiralType::CHI_UNSPECIFIED;
    }
    std::vector<Atom *> atoms;

    std::vector<std::pair<double, unsigned int>> angles;
    auto center = conf.getAtomPos(center_atom->getIdx());

    for (auto cdx_id : bond_ordering) {
      if (cdx_id == 0) {
        continue;
      }

      for (auto bond : mol.atomBonds(center_atom)) {
        int bond_id;
        if (bond->getPropIfPresent<int>(CDX_BOND_ID, bond_id)) {
        } else {
          return Atom::ChiralType::CHI_UNSPECIFIED;
        }
        if (bond_id == cdx_id) {
          auto atom = bond->getOtherAtom(center_atom);
          if (!atom) {
            // something went really wrong
            return Atom::ChiralType::CHI_UNSPECIFIED;
          }
          auto pos = conf.getAtomPos(atom->getIdx()) - center;
          double angle = atan2(pos.x, pos.y);
          angles.push_back(std::make_pair(angle, bond->getIdx()));
        }
      }
    }

    std::sort(angles.begin(), angles.end());

    // angles are now sorted in a clockwise rotation
    INT_LIST bonds;
    for (auto &angle : angles) {
      bonds.push_back(angle.second);
    }

    if (bonds.size() < 3) {
      return Atom::ChiralType::CHI_UNSPECIFIED;
    }

    auto nswaps = center_atom->getPerturbationOrder(bonds);
    if (bonds.size() == 3 && center_atom->getTotalNumHs() == 1) {
      ++nswaps;
    }
    // This is supports the HDot and HDash available in chemdraw
    //  one is an implicit wedged hydrogen and one is a dashed hydrogen
    if (center_atom->hasProp(CDX_IMPLICIT_HYDROGEN_STEREO) &&
        center_atom->getProp<char>(CDX_IMPLICIT_HYDROGEN_STEREO) == 'w') {
      nswaps++;
    }

    if (nswaps % 2) {
      return Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    }
    return Atom::ChiralType::CHI_TETRAHEDRAL_CW;
  }

  return Atom::ChiralType::CHI_UNSPECIFIED;
}
}  // namespace
void checkChemDrawTetrahedralGeometries(RWMol &mol) {
  std::vector<std::pair<char, Atom *>> unsetTetrahedralAtoms;
  Conformer *conf = nullptr;
  if (mol.getNumConformers()) {
    conf = &mol.getConformer();
  }
  bool chiralityChanged = false;

  for (auto atom : mol.atoms()) {
    // only deal with unspecified chiralities
    if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
      atom->clearProp(CDX_CIP);
      continue;
    }
    if (conf && !conf->is3D()) {
      atom->setChiralTag(getChirality(mol, atom, *conf));
      if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
        chiralityChanged = true;
      }
    }
    // If we have a cip code, might as well check it too
    CDXAtomCIPType cip;
    if (atom->getPropIfPresent<CDXAtomCIPType>(CDX_CIP, cip)) {
      // assign, possibly wrong, initial stereo.
      // note: we can probably deduce this through CDX_BOND_ORDERING, but
      //  I currently don't understand that well enough.
      switch (cip) {
        case kCDXCIPAtom_R:
          if (!chiralityChanged) {
            atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
          }
          unsetTetrahedralAtoms.push_back(std::make_pair('R', atom));
          break;
        case kCDXCIPAtom_r:
          if (!chiralityChanged) {
            atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
          }
          unsetTetrahedralAtoms.push_back(std::make_pair('r', atom));
          break;
        case kCDXCIPAtom_S:
          if (!chiralityChanged) {
            atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
          }
          unsetTetrahedralAtoms.push_back(std::make_pair('S', atom));
          break;
        case kCDXCIPAtom_s:
          if (!chiralityChanged) {
            atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
          }
          unsetTetrahedralAtoms.push_back(std::make_pair('s', atom));
          break;
        default:
          break;
      }
    }
  }

  // Now that we have missing chiralities, let's check the CIP codes and reset
  // if necessary.
  //  This is an expensive way of doing this, but we only have stereo->cip not
  //  cip->stereo implemented currently

  for (auto cipatom : unsetTetrahedralAtoms) {
    try {
      CIPLabeler::assignCIPLabels(mol);
    } catch (...) {
      // can throw std::runtime error?
      break;
    }
    std::string cipcode;
    if (cipatom.second->getPropIfPresent<std::string>(
            common_properties::_CIPCode, cipcode)) {
      if (cipcode.size() && cipcode[0] != cipatom.first) {
        // need to swap
        if (cipatom.second->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
          cipatom.second->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
          cipatom.second->updatePropertyCache();
          chiralityChanged = true;
        } else if (cipatom.second->getChiralTag() ==
                   Atom::ChiralType::CHI_TETRAHEDRAL_CCW) {
          cipatom.second->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
          cipatom.second->updatePropertyCache();
          chiralityChanged = true;
        }
      }
    }
  }
  if (chiralityChanged) {
    const bool cleanIt = true;
    const bool force = true;
    MolOps::assignStereochemistry(mol, cleanIt, force);
  }
}
}  // namespace ChemDraw
}  // namespace RDKit
