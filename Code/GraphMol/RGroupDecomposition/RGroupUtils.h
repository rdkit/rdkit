//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_UTILS
#define RGROUP_UTILS

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include "RGroupDecomp.h"

#include <map>
namespace RDKit {

RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string RLABEL;
RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string RLABEL_TYPE;
RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string RLABEL_CORE_INDEX;
RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string SIDECHAIN_RLABELS;
RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string done;
RDKIT_RGROUPDECOMPOSITION_EXPORT extern const std::string
    UNLABELED_CORE_ATTACHMENT;

const unsigned int EMPTY_CORE_LABEL = -100000;

// Various places where rgroups can be labeled
//  the order of precedence
enum class Labelling {
  RGROUP_LABELS,
  ISOTOPE_LABELS,
  ATOMMAP_LABELS,
  INDEX_LABELS,
  DUMMY_LABELS,
  INTERNAL_LABELS
};

//! return the user friendly name for the given labelling
std::string labellingToString(Labelling type);

//! Get the RLabels,atom mapping for the current molecule
std::map<int, Atom *> getRlabels(const RWMol &mol);

//! Remove the  user labels from the atom
void clearInputLabels(Atom *atom);

//! Set the rgroup label for the current atom, this also sets the
/// appropriate  MDL or other label
bool setLabel(Atom *atom, int label, std::set<int> &labels, int &maxLabel,
              bool relabel, Labelling type);

//! Returns true if the core has a dummy atom
bool hasDummy(const RWMol &core);

//! Returns true if the core atom is either an atom with multiple
/// connections or an atom with a single connection that has no user
/// defined rgroup label
bool isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(const Atom &atom);

//! Return true if the atom has a user-defined R group label
bool isUserRLabel(const Atom &atom);

// ! Return true if the atom is a terminal dummy R group (user labelled or
// unlabelled)
bool isDummyRGroupAttachment(const Atom &atom);

//! Returns true if the core atom is either a dummy atom with multiple
/// connections or a dummy atom with a single connection that has no user
/// defined rgroup label
inline bool isAnyAtomWithMultipleNeighborsOrNotUserRLabel(const Atom &atom) {
  if (atom.getAtomicNum()) {
    return false;
  }
  return isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(atom);
}

//! Returns a JSON form
/// The prefix argument is added to each line in the output
RDKIT_RGROUPDECOMPOSITION_EXPORT std::string toJSON(
    const RGroupRow &rgr, const std::string &prefix = "");
//! Returns a JSON form
/// The prefix argument is added to each line in the output
RDKIT_RGROUPDECOMPOSITION_EXPORT std::string toJSON(
    const RGroupRows &rgr, const std::string &prefix = "");
//! Returns a JSON form
/// The prefix argument is added to each line in the output
RDKIT_RGROUPDECOMPOSITION_EXPORT std::string toJSON(
    const RGroupColumn &rgr, const std::string &prefix = "");
//! Returns a JSON form
/// The prefix argument is added to each line in the output
RDKIT_RGROUPDECOMPOSITION_EXPORT std::string toJSON(
    const RGroupColumns &rgr, const std::string &prefix = "");

//! Relabel dummy atoms bearing an R-group mapping (as
/// atom map number, isotope or MDLRGroup label) such that
/// they will be displayed by the rendering code as R# rather
/// than #*, *:#, #*:#, etc. By default, only the MDLRGroup label
/// is retained on output; this may be configured through the
/// outputLabels parameter.
/// In case there are multiple potential R-group mappings,
/// the priority on input is Atom map number > Isotope > MDLRGroup.
/// The inputLabels parameter allows to configure which mappings
/// are taken into consideration.
RDKIT_RGROUPDECOMPOSITION_EXPORT void relabelMappedDummies(
    ROMol &mol, unsigned int inputLabels = AtomMap | Isotope | MDLRGroup,
    unsigned int outputLabels = MDLRGroup);

}  // namespace RDKit

#endif
