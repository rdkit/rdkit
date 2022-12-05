//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupUtils.h"
#include <boost/format.hpp>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

namespace RDKit {
std::string labellingToString(Labelling type) {
  switch (type) {
    case Labelling::RGROUP_LABELS:
      return "RGroupLabels";
    case Labelling::ISOTOPE_LABELS:
      return "IsotopeLabels";
    case Labelling::ATOMMAP_LABELS:
      return "AtomMapLabels";
    case Labelling::INDEX_LABELS:
      return "IndexLabels";
    case Labelling::DUMMY_LABELS:
      return "DummyLabels";
    case Labelling::INTERNAL_LABELS:
      return "InternalLabels";
  }
  return "unknown";
}

// Return the current set of rlabels for the molecule
//  Negative RLabels are ones not assigned by the user and are
//   either set to the index of the atom OR set to the index of
//   the atom from the core with the best MCS.
//  Positive rlabels are user defined (i.e. the user has specified
//   specific R1, R2 etc...
std::map<int, Atom *> getRlabels(const RWMol &mol) {
  std::map<int, Atom *> atoms;

  for (auto atom : mol.atoms()) {
    if (atom->hasProp(RLABEL)) {
      int rlabel = atom->getProp<int>(RLABEL);  // user label
      CHECK_INVARIANT(atoms.find(rlabel) == atoms.end(),
                      "Duplicate labels in rgroup core!");
      atoms[rlabel] = atom;
    }
  }
  return atoms;
}

void clearInputLabels(Atom *atom) {
  // atom->setIsotope(0); Don't want to clear deuterium and things like that if
  // they aren't labels
  atom->setAtomMapNum(0);
  if (atom->hasProp(common_properties::_MolFileRLabel)) {
    atom->clearProp(common_properties::_MolFileRLabel);
  }
}

bool setLabel(Atom *atom, int label, std::set<int> &labels, int &maxLabel,
              bool relabel, Labelling type) {
  if (type == Labelling::ISOTOPE_LABELS) {
    atom->setIsotope(0);
  } else if (type == Labelling::ATOMMAP_LABELS) {
    atom->setAtomMapNum(0);
  } else if (type == Labelling::RGROUP_LABELS) {
    if (atom->hasProp(common_properties::_MolFileRLabel)) {
      atom->clearProp(common_properties::_MolFileRLabel);
      atom->setIsotope(0);
    }
  }

  if (label) {
    if (labels.find(label) != labels.end()) {
      if (relabel) {
        if (type == Labelling::INTERNAL_LABELS) {
          BOOST_LOG(rdWarningLog) << "Relabelling existing label" << std::endl;
        }
        label = maxLabel;
      } else {
        // XXX FIX me - get label id
        throw ValueErrorException(
            std::string("Duplicate label in input, current type is:") +
            labellingToString(type));
      }
    }

    atom->setProp<int>(RLABEL, label);
    atom->setProp<int>(RLABEL_TYPE, static_cast<int>(type));
    labels.insert(label);
    maxLabel = (std::max)(maxLabel, label + 1);
    return true;
  }
  return false;
}

bool isUserRLabel(const Atom &atom) {
  return atom.hasProp(RLABEL) && atom.hasProp(RLABEL_TYPE) &&
         static_cast<Labelling>(atom.getProp<int>(RLABEL_TYPE)) !=
             Labelling::INDEX_LABELS;
}

bool isDummyRGroupAttachment(const Atom &atom) {
  if (atom.getAtomicNum() != 0 || atom.getDegree() != 1) {
    return false;
  }
  if (isUserRLabel(atom)) {
    return true;
  }
  int rlabel_type = 0;
  bool unlabelled_core_attachment = false;
  if (atom.hasProp(RLABEL) && atom.getPropIfPresent(RLABEL_TYPE, rlabel_type) &&
      static_cast<Labelling>(rlabel_type) == Labelling::INDEX_LABELS &&
      atom.getPropIfPresent(UNLABELLED_CORE_ATTACHMENT,
                            unlabelled_core_attachment) &&
      unlabelled_core_attachment) {
    return true;
  }
  return false;
}

bool isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(const Atom &atom) {
  if (atom.getDegree() > 1) {
    return true;
  }
  return !isDummyRGroupAttachment(atom);
}

bool hasDummy(const RWMol &core) {
  for (RWMol::ConstAtomIterator atIt = core.beginAtoms();
       atIt != core.endAtoms(); ++atIt) {
    if ((*atIt)->getAtomicNum() == 0) {
      return true;
    }
  }
  return false;
}

namespace {
std::string MolToText(const ROMol &mol) {
  bool hasQuery = false;
  for (const auto atom : mol.atoms()) {
    if (atom->hasQuery() && atom->getQuery()->getDescription() != "AtomNull") {
      hasQuery = true;
      break;
    }
  }
  if (!hasQuery) {
    for (const auto bond : mol.bonds()) {
      if (bond->hasQuery()) {
        hasQuery = true;
        break;
      }
    }
  }
  if (!hasQuery) {
    return MolToSmiles(mol);
  } else {
    return MolToSmarts(mol);
  }
}
}  // namespace

std::string toJSON(const RGroupRow &rgr, const std::string &prefix) {
  std::string res = prefix + "{\n";
  for (const auto &elem : rgr) {
    auto fmt = boost::format{"  \"%1%\":\"%2%\""} % (elem.first) %
               (MolToText(*elem.second));
    res += prefix + fmt.str() + ",\n";
  }
  res.erase(res.end() - 2, res.end());
  res += "\n" + prefix + "}";
  return res;
}

std::string toJSON(const RGroupRows &rows, const std::string &prefix) {
  std::string res = prefix + "[\n";
  auto rowPrefix = prefix + "  ";
  for (const auto &row : rows) {
    res += toJSON(row, rowPrefix) + ",\n";
  }
  res.erase(res.end() - 2, res.end());
  res += "\n" + prefix + "]";
  return res;
}

std::string toJSON(const RGroupColumn &rgr, const std::string &prefix) {
  std::string res = "[\n";
  for (const auto &elem : rgr) {
    auto fmt = boost::format{"  \"%1%\""} % (MolToText(*elem));
    res += prefix + fmt.str() + ",\n";
  }
  res.erase(res.end() - 2, res.end());
  res += "\n" + prefix + "]";
  return res;
}

std::string toJSON(const RGroupColumns &cols, const std::string &prefix) {
  std::string res = prefix + "[\n";
  auto colPrefix = prefix + "  ";
  for (const auto &col : cols) {
    auto fmt = boost::format{"  \"%1%\": %2%"} % (col.first) %
               (toJSON(col.second, colPrefix));
    res += prefix + fmt.str() + ",\n";
  }
  res.erase(res.end() - 2, res.end());
  res += "\n" + prefix + "]";
  return res;
}

}  // namespace RDKit
