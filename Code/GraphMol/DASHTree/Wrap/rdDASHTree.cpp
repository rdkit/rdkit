//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <string>
#include <vector>

#include <RDBoost/Wrap.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/DASHTree/DASHTree.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

DASH::DASHTree *makeDASHTree(const std::string &filename,
                             const python::object &properties) {
  std::vector<std::string> names;
  pythonObjectToVect<std::string>(properties, names);
  return new DASH::DASHTree(filename, names);
}

template <typename T>
python::list listFrom(const std::vector<T> &values) {
  python::list res;
  for (const auto &value : values) {
    res.append(value);
  }
  return res;
}

python::list propertyNames(const DASH::DASHTree &self) {
  return listFrom(self.getPropertyNames());
}

python::list availablePropertyNames(const DASH::DASHTree &self) {
  return listFrom(self.getAvailablePropertyNames());
}

python::list getAtomNodePath(const DASH::DASHTree &self, const ROMol &mol,
                             unsigned int atomIdx,
                             const DASH::DASHParams &params) {
  return listFrom(self.getAtomNodePath(mol, atomIdx, params));
}

python::tuple getMatchedSubstructure(const DASH::DASHTree &self,
                                     const ROMol &mol, unsigned int atomIdx,
                                     const DASH::DASHParams &params) {
  return python::tuple(
      listFrom(self.getMatchedSubstructure(mol, atomIdx, params)));
}

python::list getMolProperty(const DASH::DASHTree &self, const ROMol &mol,
                            const std::string &property,
                            const DASH::DASHParams &params) {
  return listFrom(self.getMolProperty(mol, property, params));
}

python::list getPartialCharges(const DASH::DASHTree &self, const ROMol &mol,
                               const DASH::NormalizationOptions &options) {
  return listFrom(self.getPartialCharges(mol, options));
}

python::list getNormalizedMolProperty(
    const DASH::DASHTree &self, const ROMol &mol, double target,
    const DASH::NormalizationOptions &options) {
  return listFrom(self.getNormalizedMolProperty(mol, target, options));
}

python::dict getPartialChargesDetails(
    const DASH::DASHTree &self, const ROMol &mol,
    const DASH::NormalizationOptions &options) {
  const auto [charges, rawValues, stds, matchDepths] =
      self.getPartialChargesDetails(mol, options);
  python::dict res;
  res["charges"] = listFrom(charges);
  res["raw"] = listFrom(rawValues);
  res["std"] = listFrom(stds);
  res["matchDepth"] = listFrom(matchDepths);
  return res;
}

//! borrows the molecules of \p mols, which has to outlive the result
std::vector<const ROMol *> molsFromPython(const python::list &mols) {
  std::vector<const ROMol *> res;
  pythonObjectToVect(mols, res);
  return res;
}

python::list nestedToPython(const std::vector<std::vector<double>> &values) {
  python::list res;
  for (const auto &row : values) {
    res.append(listFrom(row));
  }
  return res;
}

// A batch touches no python objects, so it runs without the GIL. A Ctrl-C
// meanwhile raises KeyboardInterrupt once python has the GIL again. The list
// keeps the molecules alive while it runs, whatever iterable they came in.

python::list getPartialChargesBatch(const DASH::DASHTree &self,
                                    const python::object &mols,
                                    const DASH::NormalizationOptions &options,
                                    int numThreads) {
  const python::list items(mols);
  const std::vector<const ROMol *> molPtrs = molsFromPython(items);
  std::vector<std::vector<double>> charges;
  {
    NOGIL gil;
    charges = self.getPartialChargesBatch(molPtrs, options, numThreads);
  }
  return nestedToPython(charges);
}

python::list getNormalizedMolPropertyBatch(
    const DASH::DASHTree &self, const python::object &mols,
    const python::object &targets, const DASH::NormalizationOptions &options,
    int numThreads) {
  const python::list items(mols);
  const std::vector<const ROMol *> molPtrs = molsFromPython(items);
  std::vector<double> targetValues;
  pythonObjectToVect(python::list(targets), targetValues);
  std::vector<std::vector<double>> values;
  {
    NOGIL gil;
    values = self.getNormalizedMolPropertyBatch(molPtrs, targetValues, options,
                                                numThreads);
  }
  return nestedToPython(values);
}

python::list getMolPropertyBatch(const DASH::DASHTree &self,
                                 const python::object &mols,
                                 const std::string &property,
                                 const DASH::DASHParams &params,
                                 int numThreads) {
  const python::list items(mols);
  const std::vector<const ROMol *> molPtrs = molsFromPython(items);
  std::vector<std::vector<double>> values;
  {
    NOGIL gil;
    values = self.getMolPropertyBatch(molPtrs, property, params, numThreads);
  }
  return nestedToPython(values);
}

python::tuple treeAtomFeature(const DASH::DASHTree &self, unsigned int index) {
  const DASH::AtomFeature &feature = self.getAtomFeature(index);
  return python::make_tuple(
      feature.atomicNum, feature.degree, feature.formalCharge,
      static_cast<bool>(feature.conjugated), feature.numHs);
}

//! \brief __getitem__ over the children
/*!
  Negative indices count from the end, as they do for a python sequence.
*/
DASH::DASHTreeNode nodeGetItem(const DASH::DASHTreeNode &self, int i) {
  const int n = static_cast<int>(self.getNumChildren());
  const int index = i < 0 ? i + n : i;
  if (index < 0 || index >= n) {
    throw_index_error(i);  // the index as the caller wrote it
  }
  return self.getChild(static_cast<unsigned int>(index));
}

const char *moduleDoc =
    "Module for DASH: per-atom properties and partial charges from a\n"
    "Dynamic Attention-based Substructure Hierarchy.\n\n"
    "DASHTree(path) opens a .dash container; rdkit.Chem.DASHTree has\n"
    "GetDASHTree() and the convert tool.\n\n"
    "NOTE: This functionality is experimental and the API and/or results may "
    "change in future releases.\n\n"
    "References:\n"
    "  M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6014\n"
    "  M. Lehner et al., J. Chem. Phys. 2024, 161, 074103\n";

}  // namespace

BOOST_PYTHON_MODULE(rdDASHTree) {
  python::scope().attr("__doc__") = moduleDoc;

  python::enum_<DASH::Normalization>("Normalization")
      .value("NONE", DASH::Normalization::NONE)
      .value("SYMMETRIC", DASH::Normalization::SYMMETRIC)
      .value("STD_WEIGHTED", DASH::Normalization::STD_WEIGHTED);

  python::class_<DASH::DASHParams>(
      "DASHParams", "Controls how a subgraph match descends into the tree.",
      python::init<>(python::args("self")))
      .def_readwrite("maxDepth", &DASH::DASHParams::maxDepth,
                     "maximum number of tree levels to descend (default 16)")
      .def_readwrite(
          "attentionThreshold", &DASH::DASHParams::attentionThreshold,
          "stop once the attention accumulated over the descent exceeds this "
          "(default 10.0)")
      .def_readwrite("attentionIncrementThreshold",
                     &DASH::DASHParams::attentionIncrementThreshold,
                     "stop once a single step contributes less attention than "
                     "this (default 0.0)")
      .def_readwrite(
          "legacyTieBreaking", &DASH::DASHParams::legacyTieBreaking,
          "where several atoms fit the next node, take the first in "
          "atom order, as the DASH-tree python package does. The result "
          "then depends on the atom order. By default, the descent "
          "follows each of them until the tree distinguishes them "
          "(default False)")
      .def("__setattr__", &safeSetattr);

  python::class_<DASH::NormalizationOptions>(
      "NormalizationOptions",
      "Controls a normalized assignment. Partial charges are the usual case.",
      python::init<>(python::args("self")))
      .def_readwrite("valueProperty",
                     &DASH::NormalizationOptions::valueProperty,
                     "property column that holds the values (default 'result')")
      .def_readwrite("stdProperty", &DASH::NormalizationOptions::stdProperty,
                     "property column that holds their standard deviations "
                     "(default 'std'). STD_WEIGHTED needs it. The others only "
                     "report it")
      .def_readwrite("normalization",
                     &DASH::NormalizationOptions::normalization,
                     "a Normalization (default STD_WEIGHTED)")
      .def_readwrite("defaultStdValue",
                     &DASH::NormalizationOptions::defaultStdValue,
                     "substituted for a stored deviation that is absent or "
                     "not positive (default 0.1)")
      .def_readwrite("params", &DASH::NormalizationOptions::params,
                     "the DASHParams that control the match")
      .def("__setattr__", &safeSetattr);

  python::class_<DASH::DASHTreeNode>(
      "DASHTreeNode",
      "One node of a DASH tree.\n\n"
      "A node describes itself relative to its parent.\n"
      "GetAtomFeatureIndex returns the atom-feature class of the atom it adds\n"
      "to the matched substructure. GetConAtom returns the position in that\n"
      "substructure of the atom it attaches to. GetConType returns the bond\n"
      "descriptor between them: 1, 2, 3 for the bond order, 4 for a\n"
      "conjugated bond. A branch root, and the heavy-atom child of a hydrogen\n"
      "root, attach to nothing. For them, GetConAtom and GetConType return\n"
      "-1.\n\n"
      "A node is a sequence of its children: len(node), node[i] and\n"
      "'for child in node' all work.\n\n"
      "Nodes read straight from the mapped file. The tree they came from\n"
      "stays alive for as long as one of its nodes exists.\n",
      python::no_init)
      .def("GetBranch", &DASH::DASHTreeNode::getBranch, python::args("self"),
           "Returns the branch the node belongs to.")
      .def("GetId", &DASH::DASHTreeNode::getId, python::args("self"),
           "Returns the node's id within its branch, as GetAtomNodePath "
           "reports it.")
      .def("GetAtomFeatureIndex", &DASH::DASHTreeNode::getAtomFeatureIndex,
           python::args("self"),
           "Returns the atom-feature class of the atom this node adds.")
      .def("GetConAtom", &DASH::DASHTreeNode::getConAtom, python::args("self"),
           "Returns the position of the substructure atom this node attaches "
           "to,\n-1 for none.")
      .def("GetConType", &DASH::DASHTreeNode::getConType, python::args("self"),
           "Returns the bond descriptor of that attachment, -1 for none.")
      .def("GetNumChildren", &DASH::DASHTreeNode::getNumChildren,
           python::args("self"), "Returns how many children the node has.")
      .def("GetChild", &DASH::DASHTreeNode::getChild,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "i"),
           "Returns the i-th child. Raises IndexError unless i is below\n"
           "GetNumChildren().")
      .def("__len__", &DASH::DASHTreeNode::getNumChildren, python::args("self"),
           "Returns the number of children.")
      .def("__getitem__", nodeGetItem,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "i"),
           "Returns the i-th child. A negative i counts from the end. Raises\n"
           "IndexError for an i past the end, which ends a for-loop over the\n"
           "node.")
      .def("GetAttention", &DASH::DASHTreeNode::getAttention,
           python::args("self"),
           "Returns the attention weight the tree assigned to this node.")
      .def("IsStopNode", &DASH::DASHTreeNode::isStopNode, python::args("self"),
           "Returns whether a descent at the default threshold stops here.")
      .def("GetValue", &DASH::DASHTreeNode::getValue,
           python::args("self", "property"),
           "Returns the value of a property at this node, NaN if it carries "
           "none.\nRaises ValueError if the property was not resolved.");

  python::class_<DASH::DASHTree, boost::noncopyable>(
      "DASHTree", "A memory-mapped DASH tree.", python::no_init)
      .def("__init__",
           python::make_constructor(
               makeDASHTree, python::default_call_policies(),
               (python::arg("filename"),
                python::arg("properties") = python::object())),
           "Maps a '.dash' container.\n\n"
           "  ARGUMENTS:\n"
           "    - filename: path to the container\n"
           "    - properties: the property columns to resolve. None or an "
           "empty\n"
           "      sequence means every column in the file\n")

      .def("GetNumBranches", &DASH::DASHTree::getNumBranches,
           python::args("self"),
           "Returns the number of atom-feature branches in the file.")
      .def("GetNumNodes", &DASH::DASHTree::getNumNodes, python::args("self"),
           "Returns the total number of tree nodes in the file.")
      .def("GetFilename", &DASH::DASHTree::getFilename,
           python::return_value_policy<python::copy_const_reference>(),
           python::args("self"), "Returns the path the tree was mapped from.")
      .def("GetPropertyNames", propertyNames, python::args("self"),
           "Returns the property columns that were resolved.")
      .def("GetAvailablePropertyNames", availablePropertyNames,
           python::args("self"),
           "Returns every property column the file holds, whether or not it "
           "was\nresolved at construction.")

      .def("GetAtomNodePath", getAtomNodePath,
           (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the nodes the match for one atom descended through: the "
           "branch\n"
           "index followed by the node ids of the descent.\n")
      .def("GetMatchedSubstructure", getMatchedSubstructure,
           (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the atom indices a match covers, as a tuple in the order "
           "the\n"
           "descent added them. With path = GetAtomNodePath(mol, atomIdx, "
           "params),\n"
           "entry i is the atom that GetNode(path[0], path[i + 1]) matched, "
           "and the\n"
           "position a node's GetConAtom refers to. A hydrogen is matched "
           "through\n"
           "its heavy neighbour. For a hydrogen, entry i is the atom of "
           "path[i + 2],\n"
           "and the hydrogen itself is not in the tuple.\n")

      .def("GetRoot", &DASH::DASHTree::getRoot,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "branch"), "Returns the root node of a branch.")
      .def("GetNode", &DASH::DASHTree::getNode,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "branch", "nodeId"),
           "Returns a node by branch and id, in the numbering GetAtomNodePath\n"
           "reports: GetNode(path[0], path[k]) is the k-th node of a match.\n")

      .def(
          "GetAtomProperty", &DASH::DASHTree::getAtomProperty,
          (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
           python::arg("property"), python::arg("params") = DASH::DASHParams()),
          "Returns the value of a property for one atom.\n\n"
          "  The deepest node of the match that carries a value gives\n"
          "  the result. Where the deepest nodes have none, a shallower,\n"
          "  more general substructure gives it. Returns NaN if no node\n"
          "  on the path carries one.\n")
      .def("GetMolProperty", getMolProperty,
           (python::arg("self"), python::arg("mol"), python::arg("property"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the value of a property for every atom of a molecule.\n")
      .def("GetNormalizedMolProperty", getNormalizedMolProperty,
           (python::arg("self"), python::arg("mol"), python::arg("target"),
            python::arg("options") = DASH::NormalizationOptions()),
           "Returns the value of options.valueProperty for every atom, "
           "adjusted\n"
           "so the values sum to target. GetPartialCharges is this with the\n"
           "molecule's formal charge as the target.\n")
      .def("GetNormalizedMolPropertyBatch", getNormalizedMolPropertyBatch,
           (python::arg("self"), python::arg("mols"), python::arg("targets"),
            python::arg("options") = DASH::NormalizationOptions(),
            python::arg("numThreads") = 1),
           "Runs GetNormalizedMolProperty for an iterable of molecules, one "
           "target\n"
           "each. See GetPartialChargesBatch for the threading semantics.\n")
      .def("GetPartialCharges", getPartialCharges,
           (python::arg("self"), python::arg("mol"),
            python::arg("options") = DASH::NormalizationOptions()),
           "Returns one partial charge per atom, normalized so they sum to "
           "the\n"
           "molecule's formal charge.\n")
      .def("GetPartialChargesDetails", getPartialChargesDetails,
           (python::arg("self"), python::arg("mol"),
            python::arg("options") = DASH::NormalizationOptions()),
           "Like GetPartialCharges, but returns a dict that also holds the raw "
           "tree\n"
           "values ('raw'), the deviations used ('std') and how deep each "
           "atom's\n"
           "match went ('matchDepth').\n")

      .def("GetPartialChargesBatch", getPartialChargesBatch,
           (python::arg("self"), python::arg("mols"),
            python::arg("options") = DASH::NormalizationOptions(),
            python::arg("numThreads") = 1),
           "Returns partial charges for an iterable of molecules.\n\n"
           "  ARGUMENTS:\n"
           "    - mols: an iterable of molecules. None entries give empty "
           "results\n"
           "    - numThreads: how many threads to use. 0 means every hardware\n"
           "      thread. A negative value leaves that many of them unused.\n\n"
           "  Each thread takes one molecule at a time, so a batch of unequal\n"
           "  molecules still balances across the threads.\n")
      .def("GetMolPropertyBatch", getMolPropertyBatch,
           (python::arg("self"), python::arg("mols"), python::arg("property"),
            python::arg("params") = DASH::DASHParams(),
            python::arg("numThreads") = 1),
           "Returns the value of a property for every atom of every molecule "
           "in\n"
           "an iterable. See GetPartialChargesBatch for the threading "
           "semantics.\n")

      .def("GetAtomFeature", treeAtomFeature,
           (python::arg("self"), python::arg("index")),
           "Returns the atom-feature class with the given branch index as\n"
           "(atomicNum, degree, formalCharge, conjugated, numHs).\n")
      .def("GetAtomFeatureIndex", &DASH::DASHTree::getAtomFeatureIndex,
           (python::arg("self"), python::arg("atom")),
           "Returns the branch index of an atom, or -1 if the atom is in none\n"
           "of the tree's atom-feature classes.\n");
}
