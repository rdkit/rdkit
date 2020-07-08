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
#include <map>
namespace RDKit {
extern const std::string RLABEL;
extern const std::string SIDECHAIN_RLABELS;
extern const std::string done;

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
//  appropriate  MDL or other label
bool setLabel(Atom *atom, int label, std::set<int> &labels, int &maxLabel,
              bool relabel, Labelling type);

//! Returns true if the core has a dummy atom
bool hasDummy(const RWMol &core);
}

#endif
