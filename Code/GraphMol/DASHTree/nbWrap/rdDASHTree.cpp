//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <RDBoost/Wrap_nb.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/DASHTree/DASHTree.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RDKit;
using namespace nb::literals;

namespace {

nb::tuple featureToPython(const DASH::AtomFeature &feature) {
  return nb::make_tuple(feature.atomicNum, feature.degree, feature.formalCharge,
                        static_cast<bool>(feature.conjugated), feature.numHs);
}

//! borrows the molecules a python list holds; None becomes a null entry
std::vector<const ROMol *> molsFromPython(const nb::list &mols) {
  std::vector<const ROMol *> res;
  for (nb::handle item : mols) {
    res.push_back(item.is_none() ? nullptr : nb::cast<const ROMol *>(item));
  }
  return res;
}

const char *moduleDoc =
    "Module containing DASH: per-atom properties and partial charges from a\n"
    "Dynamic Attention-based Substructure Hierarchy.\n\n"
    "DASHTree(path) opens a .dash container; rdkit.Chem.DASHTree has\n"
    "GetDASHTree() and the convert tool.\n\n"
    "NOTE: This functionality is experimental and the API and/or results may "
    "change in future releases.\n\n"
    "References:\n"
    "  M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6014\n"
    "  M. Lehner et al., J. Chem. Phys. 2024, 161, 074103\n";

}  // namespace

NB_MODULE(rdDASHTree, m) {
  m.doc() = moduleDoc;

  nb::enum_<DASH::Normalization>(m, "Normalization")
      .value("NONE", DASH::Normalization::NONE)
      .value("SYMMETRIC", DASH::Normalization::SYMMETRIC)
      .value("STD_WEIGHTED", DASH::Normalization::STD_WEIGHTED);

  nb::class_<DASH::DASHParams>(
      m, "DASHParams", "Controls how a subgraph match descends into the tree.")
      .def(nb::init<>())
      .def_rw("maxDepth", &DASH::DASHParams::maxDepth,
              "maximum number of tree levels to descend (default 16)")
      .def_rw("attentionThreshold", &DASH::DASHParams::attentionThreshold,
              "stop once the attention accumulated over the descent exceeds "
              "this (default 10.0)")
      .def_rw("attentionIncrementThreshold",
              &DASH::DASHParams::attentionIncrementThreshold,
              "stop once a single step contributes less attention than this "
              "(default 0.0)")
      .def_rw("legacyTieBreaking", &DASH::DASHParams::legacyTieBreaking,
              "where several atoms fit the next node, take the first in atom "
              "order, as the DASH-tree python package does, instead of "
              "following each until the tree tells them apart; the result "
              "then depends on the atom order (default False)")
      .def("__setattr__", &safeSetattr);

  nb::class_<DASH::NormalizationOptions>(
      m, "NormalizationOptions",
      "Controls a normalised assignment, partial charges being the usual "
      "case.")
      .def(nb::init<>())
      .def_rw("valueProperty", &DASH::NormalizationOptions::valueProperty,
              "property column holding the values (default 'result')")
      .def_rw("stdProperty", &DASH::NormalizationOptions::stdProperty,
              "property column holding their standard deviations (default "
              "'std'); needed by STD_WEIGHTED, the others only report it")
      .def_rw("normalization", &DASH::NormalizationOptions::normalization,
              "a Normalization (default STD_WEIGHTED)")
      .def_rw("defaultStdValue", &DASH::NormalizationOptions::defaultStdValue,
              "substituted for a stored deviation that is not positive or "
              "absent (default 0.1)")
      .def_rw("params", &DASH::NormalizationOptions::params,
              "the DASHParams controlling the match")
      .def("__setattr__", &safeSetattr);

  nb::class_<DASH::DASHTreeNode>(
      m, "DASHTreeNode",
      "One node of a DASH tree.\n\n"
      "A node describes itself relative to its parent: the atom-feature "
      "class of\n"
      "the atom it adds to the matched substructure (GetAtomFeatureIndex),\n"
      "the position in that substructure of the atom it attaches to "
      "(GetConAtom),\n"
      "and the bond descriptor between them (GetConType: 1, 2, 3 for the "
      "bond\n"
      "order, 4 for a conjugated bond). A branch root, and the heavy-atom "
      "child of\n"
      "a hydrogen root, attach to nothing and report -1 for both.\n\n"
      "A node is a sequence of its children: len(node), node[i] and\n"
      "'for child in node' all work.\n\n"
      "Handles read straight from the mapped file and are cheap to copy; the "
      "tree\n"
      "they came from is kept alive for as long as one exists.\n")
      .def("GetBranch", &DASH::DASHTreeNode::getBranch,
           "Returns the branch the node belongs to.")
      .def("GetId", &DASH::DASHTreeNode::getId,
           "Returns the node's id within its branch, as GetAtomNodePath "
           "reports it.")
      .def("GetAtomFeatureIndex", &DASH::DASHTreeNode::getAtomFeatureIndex,
           "Returns the atom-feature class of the atom this node adds.")
      .def("GetConAtom", &DASH::DASHTreeNode::getConAtom,
           "Returns the position of the substructure atom this node attaches "
           "to,\n-1 for none.")
      .def("GetConType", &DASH::DASHTreeNode::getConType,
           "Returns the bond descriptor of that attachment, -1 for none.")
      .def("GetNumChildren", &DASH::DASHTreeNode::getNumChildren,
           "Returns how many children the node has.")
      .def("GetChild", &DASH::DASHTreeNode::getChild, "i"_a,
           nb::keep_alive<0, 1>(),
           "Returns the i-th child. Raises IndexError unless i is below\n"
           "GetNumChildren().")
      .def("__len__", &DASH::DASHTreeNode::getNumChildren,
           "The number of children.")
      .def(
          "__getitem__",
          [](const DASH::DASHTreeNode &self, int i) {
            const int n = static_cast<int>(self.getNumChildren());
            const int index = i < 0 ? i + n : i;
            if (index < 0 || index >= n) {
              // std::out_of_range is what nanobind turns into an IndexError,
              // which is what ends a for-loop over the node
              throw std::out_of_range(std::to_string(i));
            }
            return self.getChild(static_cast<unsigned int>(index));
          },
          "i"_a, nb::keep_alive<0, 1>(),
          "The i-th child, counting from the end for negative i. IndexError "
          "past\nthe end, which is what ends a for-loop over the node.")
      .def("GetAttention", &DASH::DASHTreeNode::getAttention,
           "Returns the attention weight the tree assigned to this node.")
      .def("IsStopNode", &DASH::DASHTreeNode::isStopNode,
           "Returns whether a descent at the default threshold stops here.")
      .def("GetValue", &DASH::DASHTreeNode::getValue, "property"_a,
           "Returns the value of a property at this node, NaN if it carries "
           "none.\nRaises ValueError if the property was not resolved.");

  nb::class_<DASH::DASHTree>(m, "DASHTree", "A memory-mapped DASH tree.")
      .def(
          "__init__",
          [](DASH::DASHTree *self, const std::string &filename,
             nb::object properties) {
            std::vector<std::string> names;
            pythonObjectToVect<std::string>(properties, names);
            new (self) DASH::DASHTree(filename, names);
          },
          "filename"_a, "properties"_a = nb::none(),
          "Maps a '.dash' container.\n\n"
          "  ARGUMENTS:\n"
          "    - filename: path to the container\n"
          "    - properties: the property columns to resolve; None or an "
          "empty\n"
          "      sequence means every column in the file\n")

      .def("GetNumBranches", &DASH::DASHTree::getNumBranches,
           "Returns the number of atom-feature branches in the file.")
      .def("GetNumNodes", &DASH::DASHTree::getNumNodes,
           "Returns the total number of tree nodes in the file.")
      .def("GetFilename", &DASH::DASHTree::getFilename,
           "Returns the path the tree was mapped from.")
      .def("GetPropertyNames", &DASH::DASHTree::getPropertyNames,
           "Returns the property columns that were resolved.")
      .def("GetAvailablePropertyNames",
           &DASH::DASHTree::getAvailablePropertyNames,
           "Returns every property column the file holds, whether or not it "
           "was\nresolved at construction.")

      .def(
          "GetAtomNodePath",
          [](const DASH::DASHTree &self, const ROMol &mol, unsigned int atomIdx,
             const DASH::DASHParams &params) {
            std::vector<std::uint32_t> path;
            self.getAtomNodePath(mol, atomIdx, path, params);
            return path;
          },
          "mol"_a, "atomIdx"_a, "params"_a = DASH::DASHParams(),
          "Returns the nodes the match for one atom descended through: the "
          "branch\n"
          "index followed by the node ids of the descent.\n")
      .def(
          "GetMatchedSubstructure",
          [](const DASH::DASHTree &self, const ROMol &mol, unsigned int atomIdx,
             const DASH::DASHParams &params) {
            std::vector<unsigned int> atoms;
            self.getMatchedSubstructure(mol, atomIdx, atoms, params);
            return nb::tuple(nb::cast(atoms));
          },
          "mol"_a, "atomIdx"_a, "params"_a = DASH::DASHParams(),
          "Returns the atom indices a match covers, as a tuple in the order "
          "the\n"
          "descent added them: entry i is the atom that node i of the path\n"
          "matched, and the position a node's GetConAtom refers to. A "
          "hydrogen\n"
          "is matched through its heavy neighbour, so for one the tuple "
          "starts at\n"
          "that neighbour.\n")

      .def("GetRoot", &DASH::DASHTree::getRoot, "branch"_a,
           nb::keep_alive<0, 1>(), "Returns the root node of a branch.")
      .def("GetNode", &DASH::DASHTree::getNode, "branch"_a, "nodeId"_a,
           nb::keep_alive<0, 1>(),
           "Returns a node by branch and id, in the numbering GetAtomNodePath\n"
           "reports: GetNode(path[0], path[k]) is the k-th node of a match.\n")

      .def("GetAtomProperty", &DASH::DASHTree::getAtomProperty, "mol"_a,
           "atomIdx"_a, "property"_a, "params"_a = DASH::DASHParams(),
           "Returns the value of a property for one atom.\n\n"
           "  The deepest node of the match that carries a value wins, so a\n"
           "  sparsely populated column falls back to a more general\n"
           "  substructure. NaN if no node on the path carries one.\n")
      .def(
          "GetMolProperty",
          [](const DASH::DASHTree &self, const ROMol &mol,
             const std::string &property, const DASH::DASHParams &params) {
            std::vector<double> values;
            self.getMolProperty(mol, property, values, params);
            return values;
          },
          "mol"_a, "property"_a, "params"_a = DASH::DASHParams(),
          "Returns the value of a property for every atom of a molecule.\n")
      .def(
          "GetNormalizedMolProperty",
          [](const DASH::DASHTree &self, const ROMol &mol, double target,
             const DASH::NormalizationOptions &options) {
            std::vector<double> values;
            self.getNormalizedMolProperty(mol, target, values, options);
            return values;
          },
          "mol"_a, "target"_a, "options"_a = DASH::NormalizationOptions(),
          "Returns the value of options.valueProperty for every atom, "
          "adjusted\n"
          "so the values sum to target. GetPartialCharges is this with the\n"
          "molecule's formal charge as the target.\n")
      .def(
          "GetNormalizedMolPropertyBatch",
          [](const DASH::DASHTree &self, const nb::object &mols,
             const std::vector<double> &targets,
             const DASH::NormalizationOptions &options, int numThreads) {
            const nb::list items(mols);
            const std::vector<const ROMol *> molPtrs = molsFromPython(items);
            std::vector<std::vector<double>> values;
            {
              NOGIL gil;
              self.getNormalizedMolPropertyBatch(molPtrs, targets, values,
                                                 options, numThreads);
            }
            return values;
          },
          "mols"_a, "targets"_a, "options"_a = DASH::NormalizationOptions(),
          "numThreads"_a = 1,
          "GetNormalizedMolProperty for an iterable of molecules, one target "
          "each.\n"
          "See GetPartialChargesBatch for the threading semantics.\n")
      .def(
          "GetPartialCharges",
          [](const DASH::DASHTree &self, const ROMol &mol,
             const DASH::NormalizationOptions &options) {
            std::vector<double> charges;
            self.getPartialCharges(mol, charges, options);
            return charges;
          },
          "mol"_a, "options"_a = DASH::NormalizationOptions(),
          "Returns one partial charge per atom, normalised so they sum to "
          "the\n"
          "molecule's formal charge.\n")
      .def(
          "GetPartialChargesDetails",
          [](const DASH::DASHTree &self, const ROMol &mol,
             const DASH::NormalizationOptions &options) {
            std::vector<double> charges, rawValues, stds;
            std::vector<unsigned int> matchDepths;
            self.getPartialCharges(mol, charges, rawValues, stds, matchDepths,
                                   options);
            nb::dict res;
            res["charges"] = nb::cast(charges);
            res["raw"] = nb::cast(rawValues);
            res["std"] = nb::cast(stds);
            res["matchDepth"] = nb::cast(matchDepths);
            return res;
          },
          "mol"_a, "options"_a = DASH::NormalizationOptions(),
          "Like GetPartialCharges, but returns a dict also holding the raw "
          "tree\n"
          "values ('raw'), the deviations used ('std') and how deep each "
          "atom's\n"
          "match went ('matchDepth').\n")

      .def(
          "GetPartialChargesBatch",
          [](const DASH::DASHTree &self, const nb::object &mols,
             const DASH::NormalizationOptions &options, int numThreads) {
            // the list keeps the molecules alive while the GIL is released
            const nb::list items(mols);
            const std::vector<const ROMol *> molPtrs = molsFromPython(items);
            std::vector<std::vector<double>> charges;
            {
              // the batch touches no python objects, so the GIL goes; python
              // raises a pending KeyboardInterrupt once it is back
              NOGIL gil;
              self.getPartialChargesBatch(molPtrs, charges, options,
                                          numThreads);
            }
            return charges;
          },
          "mols"_a, "options"_a = DASH::NormalizationOptions(),
          "numThreads"_a = 1,
          "Returns partial charges for an iterable of molecules.\n\n"
          "  ARGUMENTS:\n"
          "    - mols: an iterable of molecules; None entries give empty "
          "results\n"
          "    - numThreads: how many threads to use. 0 means every hardware\n"
          "      thread; a negative value leaves that many of them unused.\n\n"
          "  Molecules are handed out one at a time, so a batch of unequal\n"
          "  molecules still balances across the threads.\n")
      .def(
          "GetMolPropertyBatch",
          [](const DASH::DASHTree &self, const nb::object &mols,
             const std::string &property, const DASH::DASHParams &params,
             int numThreads) {
            const nb::list items(mols);
            const std::vector<const ROMol *> molPtrs = molsFromPython(items);
            std::vector<std::vector<double>> values;
            {
              NOGIL gil;
              self.getMolPropertyBatch(molPtrs, property, values, params,
                                       numThreads);
            }
            return values;
          },
          "mols"_a, "property"_a, "params"_a = DASH::DASHParams(),
          "numThreads"_a = 1,
          "Returns the value of a property for every atom of every molecule "
          "in\n"
          "an iterable. See GetPartialChargesBatch for the threading "
          "semantics.\n")

      .def(
          "GetAtomFeature",
          [](const DASH::DASHTree &self, unsigned int index) {
            return featureToPython(self.getAtomFeature(index));
          },
          "index"_a,
          "Returns the atom-feature class with the given branch index as\n"
          "(atomicNum, degree, formalCharge, conjugated, numHs).\n")
      .def("GetAtomFeatureIndex", &DASH::DASHTree::getAtomFeatureIndex,
           "atom"_a,
           "Returns the branch index of an atom, or -1 if its atom type is "
           "not\n"
           "one of the classes this tree covers.\n");
}
