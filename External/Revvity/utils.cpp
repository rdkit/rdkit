#include "utils.h"

namespace RDKit {

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
      return "Forumla";
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
    // If we don't have a bond length for any reason, just scale the avgerage
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
    // C = Conneciton
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

  bool replace(RWMol& mol) {
    if (!replacement_atom) return true;

    auto &bond_ordering =
        replacement_atom->getProp<std::vector<int>>(CDX_BOND_ORDERING);
    
    // Find the connecting atoms and and do the replacement
    size_t i = 0;  
    for (auto bond : mol.atomBonds(replacement_atom)) {
        // find the position of the attachement bonds in the bond ordering
      auto bond_id = bond->getProp<unsigned int>(CDX_BOND_ID);
      auto it = std::find(bond_ordering.begin(), bond_ordering.end(), bond_id);
      if (it == bond_ordering.end()) return false;

      auto pos = std::distance(bond_ordering.begin(), it);

      auto &xatom = fragment_atoms[pos];

      for (auto &xbond : mol.atomBonds(xatom)) {
          // xatom is the fragment dummy atom
          // xbond is the fragment bond
        if (bond->getBeginAtom() == replacement_atom) {
          mol.addBond(xbond->getOtherAtom(xatom), bond->getEndAtom(),
                      bond->getBondType());
        } else {
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

}  // namespace RDKit
