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

#include <memory>
#include <string>
#include <vector>

#include <RDBoost/Wrap.h>
#include <RDGeneral/ControlCHandler.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/DASHTree/AtomFeatures.h>
#include <GraphMol/DASHTree/DASHTree.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

//! turns a python sequence of strings into a vector
std::vector<std::string> stringsFromPython(const python::object &seq) {
  std::vector<std::string> res;
  if (seq.is_none()) {
    return res;
  }
  python::stl_input_iterator<std::string> begin(seq), end;
  res.insert(res.end(), begin, end);
  return res;
}

DASH::DASHTree *makeDASHTree(const std::string &filename,
                             const python::object &properties, bool prefetch) {
  return new DASH::DASHTree(filename, stringsFromPython(properties), prefetch);
}

python::list doublesToPython(const std::vector<double> &values) {
  python::list res;
  for (const auto value : values) {
    res.append(value);
  }
  return res;
}

python::list propertyNames(const DASH::DASHTree &self) {
  python::list res;
  for (const auto &name : self.propertyNames()) {
    res.append(name);
  }
  return res;
}

python::list availablePropertyNames(const DASH::DASHTree &self) {
  python::list res;
  for (const auto &name : self.availablePropertyNames()) {
    res.append(name);
  }
  return res;
}

python::list getAtomNodePath(const DASH::DASHTree &self, const ROMol &mol,
                             unsigned int atomIdx,
                             const DASH::DASHParams &params,
                             bool sourceNodeIds) {
  std::vector<std::uint32_t> path;
  self.getAtomNodePath(mol, atomIdx, path, params, sourceNodeIds);
  python::list res;
  for (const auto node : path) {
    res.append(node);
  }
  return res;
}

python::tuple getMatchedSubstructure(const DASH::DASHTree &self,
                                     const ROMol &mol, unsigned int atomIdx,
                                     const DASH::DASHParams &params) {
  std::vector<unsigned int> atoms;
  self.getMatchedSubstructure(mol, atomIdx, atoms, params);
  python::list res;
  for (const auto atom : atoms) {
    res.append(atom);
  }
  return python::tuple(res);
}

python::list getMolProperty(const DASH::DASHTree &self, const ROMol &mol,
                            const std::string &property,
                            const DASH::DASHParams &params) {
  std::vector<double> values;
  self.getMolProperty(mol, property, values, params);
  return doublesToPython(values);
}

python::list getPartialCharges(const DASH::DASHTree &self, const ROMol &mol,
                               const DASH::ChargeOptions &options) {
  std::vector<double> charges;
  self.getPartialCharges(mol, charges, options);
  return doublesToPython(charges);
}

python::dict getPartialChargesDetails(const DASH::DASHTree &self,
                                      const ROMol &mol,
                                      const DASH::ChargeOptions &options) {
  std::vector<double> charges, rawValues, stds;
  std::vector<unsigned int> matchDepths;
  self.getPartialCharges(mol, charges, rawValues, stds, matchDepths, options);
  python::list depths;
  for (const auto depth : matchDepths) {
    depths.append(depth);
  }
  python::dict res;
  res["charges"] = doublesToPython(charges);
  res["raw"] = doublesToPython(rawValues);
  res["std"] = doublesToPython(stds);
  res["match_depth"] = depths;
  return res;
}

//! borrows the molecules a python sequence holds, without copying them
std::vector<const ROMol *> molsFromPython(const python::object &seq) {
  std::vector<const ROMol *> mols;
  const unsigned int n = python::extract<unsigned int>(seq.attr("__len__")());
  mols.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    python::object item = seq[i];
    if (item.is_none()) {
      mols.push_back(nullptr);
      continue;
    }
    mols.push_back(python::extract<const ROMol *>(item));
  }
  return mols;
}

python::list nestedToPython(const std::vector<std::vector<double>> &values) {
  python::list res;
  for (const auto &row : values) {
    res.append(doublesToPython(row));
  }
  return res;
}

python::list getPartialChargesBatch(const DASH::DASHTree &self,
                                    const python::object &mols,
                                    const DASH::ChargeOptions &options,
                                    int numThreads) {
  const std::vector<const ROMol *> molPtrs = molsFromPython(mols);
  std::vector<std::vector<double>> charges;
  bool interrupted = false;
  {
    // A batch can run for a long time and touches no python objects, so the
    // GIL goes; Ctrl-C has to be raised after it comes back.
    NOGIL gil;
    ControlCHandler::reset();
    self.getPartialChargesBatch(molPtrs, charges, options, numThreads);
    interrupted = ControlCHandler::getGotSignal();
  }
  if (interrupted) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "DASH tree batch cancelled");
    python::throw_error_already_set();
  }
  return nestedToPython(charges);
}

python::list getMolPropertyBatch(const DASH::DASHTree &self,
                                 const python::object &mols,
                                 const std::string &property,
                                 const DASH::DASHParams &params,
                                 int numThreads) {
  const std::vector<const ROMol *> molPtrs = molsFromPython(mols);
  std::vector<std::vector<double>> values;
  bool interrupted = false;
  {
    NOGIL gil;
    ControlCHandler::reset();
    self.getMolPropertyBatch(molPtrs, property, values, params, numThreads);
    interrupted = ControlCHandler::getGotSignal();
  }
  if (interrupted) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "DASH tree batch cancelled");
    python::throw_error_already_set();
  }
  return nestedToPython(values);
}

python::tuple featureToPython(const DASH::AtomFeature &feature) {
  return python::make_tuple(feature.atomicNum, feature.degree,
                            feature.formalCharge, feature.conjugated,
                            feature.numHs);
}

python::tuple atomFeature(unsigned int index) {
  if (index >= DASH::numAtomFeatures) {
    throw_value_error("atom feature index out of range");
  }
  return featureToPython(DASH::getAtomFeatureTable()[index]);
}

python::tuple nodeFeature(const DASH::DASHTreeNode &self) {
  return featureToPython(self.getFeature());
}

//! \brief __getitem__ over the children
/*!
  Negative indices count from the end. Past the end is an IndexError rather
  than the ValueError GetChild raises, because that is what ends a for-loop
  over the node.
*/
DASH::DASHTreeNode nodeGetItem(const DASH::DASHTreeNode &self, int i) {
  const int n = static_cast<int>(self.getNumChildren());
  const int index = i < 0 ? i + n : i;
  if (index < 0 || index >= n) {
    throw_index_error(i);  // the index as the caller wrote it
  }
  return self.getChild(static_cast<unsigned int>(index));
}

int atomFeatureIndexForAtom(const Atom *atom) {
  return DASH::atomFeatureIndex(atom);
}

const char *moduleDoc =
    "DASH: per-atom properties and partial charges from a Dynamic "
    "Attention-based Substructure Hierarchy.\n\n"
    "The tree itself is a data file. Convert a DASH-tree distribution with\n"
    "Code/GraphMol/DASHTree/tools/dash_convert.py and hand the resulting\n"
    "'.dash' file to DASHTree:\n\n"
    "  >>> from rdkit import Chem\n"
    "  >>> from rdkit.Chem import rdDASHTree\n"
    "  >>> tree = rdDASHTree.DASHTree('default.dash', ['result', 'std'])\n"
    "  >>> mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))\n"
    "  >>> charges = tree.GetPartialCharges(mol)\n\n"
    "The file is memory mapped, so constructing a DASHTree is near-instant "
    "and\n"
    "only the parts of the tree a query visits are ever read from disk.\n\n"
    "References:\n"
    "  M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6296\n"
    "  M. Lehner et al., J. Chem. Phys. 2024, 161, 044113\n";

}  // namespace

BOOST_PYTHON_MODULE(rdDASHTree) {
  python::scope().attr("__doc__") = moduleDoc;

  RegisterVectorConverter<std::string>("_vectstring_dash");

  python::enum_<DASH::ChargeNormalization>("ChargeNormalization")
      .value("NONE", DASH::ChargeNormalization::NONE)
      .value("SYMMETRIC", DASH::ChargeNormalization::SYMMETRIC)
      .value("STD_WEIGHTED", DASH::ChargeNormalization::STD_WEIGHTED);

  python::class_<DASH::DASHParams>(
      "DASHParams",
      "Controls how far a subgraph match descends into the tree.",
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
                     "this (default 0.0)");

  python::class_<DASH::ChargeOptions>(
      "ChargeOptions", "Controls partial-charge assignment.",
      python::init<>(python::args("self")))
      .def_readwrite("valueProperty", &DASH::ChargeOptions::valueProperty,
                     "property column holding the charges (default 'result')")
      .def_readwrite(
          "stdProperty", &DASH::ChargeOptions::stdProperty,
          "property column holding their standard deviations (default 'std')")
      .def_readwrite("normalization", &DASH::ChargeOptions::normalization,
                     "a ChargeNormalization (default STD_WEIGHTED)")
      .def_readwrite("defaultStdValue", &DASH::ChargeOptions::defaultStdValue,
                     "substituted for a stored deviation that is not positive "
                     "(default 0.1)")
      .def_readwrite("params", &DASH::ChargeOptions::params,
                     "the DASHParams controlling the match");

  python::class_<DASH::DASHTreeNode>(
      "DASHTreeNode",
      "One node of a DASH tree.\n\n"
      "A node describes itself relative to its parent: the atom-feature "
      "class of\n"
      "the atom it adds to the matched substructure (GetAtomFeatureIndex, "
      "GetFeature),\n"
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
      "they came from is kept alive for as long as one exists.\n",
      python::no_init)
      .def("GetBranch", &DASH::DASHTreeNode::getBranch, python::args("self"),
           "Returns the branch the node belongs to.")
      .def("GetId", &DASH::DASHTreeNode::getId, python::args("self"),
           "Returns the node's id within its branch, as GetAtomNodePath "
           "reports it.")
      .def("GetKey", &DASH::DASHTreeNode::getKey, python::args("self"),
           "Returns the packed 16-bit match key.")
      .def("GetAtomFeatureIndex", &DASH::DASHTreeNode::getAtomFeatureIndex,
           python::args("self"),
           "Returns the atom-feature class of the atom this node adds.")
      .def("GetFeature", nodeFeature, python::args("self"),
           "Returns that class as (atomicNum, degree, formalCharge, "
           "conjugated, numHs).")
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
           "Returns the i-th child. Raises ValueError from GetNumChildren() "
           "on.")
      .def("__len__", &DASH::DASHTreeNode::getNumChildren,
           python::args("self"), "The number of children.")
      .def("__getitem__", nodeGetItem,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "i"),
           "The i-th child, counting from the end for negative i. IndexError "
           "past\nthe end, which is what ends a for-loop over the node.")
      .def("GetAttention", &DASH::DASHTreeNode::getAttention,
           python::args("self"),
           "Returns the attention weight the tree assigned to this node.")
      .def("Stops", &DASH::DASHTreeNode::stops, python::args("self"),
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
                python::arg("properties") = python::object(),
                python::arg("prefetch") = false)),
           "Maps a '.dash' container.\n\n"
           "  ARGUMENTS:\n"
           "    - filename: path to the container\n"
           "    - properties: the property columns to resolve; None or an "
           "empty\n"
           "      sequence means every column in the file\n"
           "    - prefetch: read through the mapped arrays once so later "
           "queries\n"
           "      do not pay page faults. Worth it before a large batch, "
           "pointless\n"
           "      for a handful of molecules.\n")

      .def("GetNumBranches", &DASH::DASHTree::numBranches, python::args("self"),
           "Returns the number of atom-feature branches in the file.")
      .def("GetNumNodes", &DASH::DASHTree::numNodes, python::args("self"),
           "Returns the total number of tree nodes in the file.")
      .def("GetMappedSize", &DASH::DASHTree::mappedSize, python::args("self"),
           "Returns how many bytes of the file were mapped.")
      .def("GetFileName", &DASH::DASHTree::filename,
           python::return_value_policy<python::copy_const_reference>(),
           python::args("self"),
           "Returns the path the tree was mapped from.")
      .def("GetPropertyNames", propertyNames, python::args("self"),
           "Returns the property columns that were resolved.")
      .def("HasProperty", &DASH::DASHTree::hasProperty,
           python::args("self", "name"),
           "Returns whether a property column was resolved.")

      .def("HasSourceNodeIds", &DASH::DASHTree::hasSourceNodeIds,
           python::args("self"),
           "Returns whether the file carries the map back to the source tree's "
           "node\nnumbering, which GetAtomNodePath can report instead of the "
           "container's own.")
      .def("GetAvailablePropertyNames", availablePropertyNames,
           python::args("self"),
           "Returns every property column the file holds, whether or not it "
           "was\nresolved at construction.")

      .def("GetAtomNodePath", getAtomNodePath,
           (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
            python::arg("params") = DASH::DASHParams(),
            python::arg("sourceNodeIds") = false),
           "Returns the nodes the match for one atom descended through: the "
           "branch\n"
           "index followed by the node ids of the descent.\n\n"
           "  Pass sourceNodeIds=True for the numbering the DASH-tree python\n"
           "  package uses rather than the container's own. That needs a file\n"
           "  carrying it -- see HasSourceNodeIds().\n")
      .def("GetMatchedSubstructure", getMatchedSubstructure,
           (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the atom indices a match covers, as a tuple in the order "
           "the\n"
           "descent added them: entry i is the atom that node i of the path\n"
           "matched, and the position a node's GetConAtom refers to. A "
           "hydrogen\n"
           "is matched through its heavy neighbour, so for one the tuple "
           "starts at\n"
           "that neighbour.\n")

      .def("GetRoot", &DASH::DASHTree::getRoot,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "branch"),
           "Returns the root node of a branch.")
      .def("GetNode", &DASH::DASHTree::getNode,
           python::with_custodian_and_ward_postcall<0, 1>(),
           python::args("self", "branch", "nodeId"),
           "Returns a node by branch and id, in the numbering GetAtomNodePath\n"
           "reports: GetNode(path[0], path[k]) is the k-th node of a match.\n")

      .def("GetAtomProperty", &DASH::DASHTree::getAtomProperty,
           (python::arg("self"), python::arg("mol"), python::arg("atomIdx"),
            python::arg("property"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the value of a property for one atom.\n\n"
           "  The deepest node of the match that carries a value wins, so a\n"
           "  sparsely populated column falls back to a more general\n"
           "  substructure. NaN if no node on the path carries one.\n")
      .def("GetMolProperty", getMolProperty,
           (python::arg("self"), python::arg("mol"), python::arg("property"),
            python::arg("params") = DASH::DASHParams()),
           "Returns the value of a property for every atom of a molecule.\n")
      .def("GetPartialCharges", getPartialCharges,
           (python::arg("self"), python::arg("mol"),
            python::arg("options") = DASH::ChargeOptions()),
           "Returns one partial charge per atom, normalised so they sum to "
           "the\n"
           "molecule's formal charge.\n")
      .def("GetPartialChargesDetails", getPartialChargesDetails,
           (python::arg("self"), python::arg("mol"),
            python::arg("options") = DASH::ChargeOptions()),
           "Like GetPartialCharges, but returns a dict also holding the raw "
           "tree\n"
           "values ('raw'), the deviations used ('std') and how deep each "
           "atom's\n"
           "match went ('match_depth').\n")

      .def("GetPartialChargesBatch", getPartialChargesBatch,
           (python::arg("self"), python::arg("mols"),
            python::arg("options") = DASH::ChargeOptions(),
            python::arg("numThreads") = 1),
           "Returns partial charges for a sequence of molecules.\n\n"
           "  ARGUMENTS:\n"
           "    - mols: a sequence of molecules; None entries give empty "
           "results\n"
           "    - numThreads: how many threads to use. Negative counts back "
           "from\n"
           "      the number of hardware threads, so -1 means all of them.\n\n"
           "  Molecules are handed out one at a time, so a batch of unequal\n"
           "  molecules still balances across the threads.\n")
      .def("GetMolPropertyBatch", getMolPropertyBatch,
           (python::arg("self"), python::arg("mols"), python::arg("property"),
            python::arg("params") = DASH::DASHParams(),
            python::arg("numThreads") = 1),
           "Returns the value of a property for every atom of every molecule "
           "in\n"
           "a sequence. See GetPartialChargesBatch for the threading "
           "semantics.\n");

  python::def("GetNumAtomFeatures",
              +[]() { return DASH::numAtomFeatures; },
              "Returns the number of atom-feature classes the DASH trees are "
              "built from.");
  python::def("GetAtomFeature", atomFeature, python::args("index"),
              "Returns the atom-feature class with the given branch index as\n"
              "(atomicNum, degree, formalCharge, conjugated, numHs).\n");
  python::def("GetAtomFeatureIndex", atomFeatureIndexForAtom,
              python::args("atom"),
              "Returns the DASH branch index of an atom, or -1 if its atom "
              "type\n"
              "is not one of the classes the DASH trees cover.\n");
}
