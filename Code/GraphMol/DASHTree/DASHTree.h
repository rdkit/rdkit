//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file DASHTree.h
/*!
  \brief DASH: per-atom properties and partial charges from a Dynamic
         Attention-based Substructure Hierarchy

  \b Note that this functionality is experimental and the API and/or results
  may change in future releases.

  For each atom, DASH walks a precomputed decision tree whose levels grow the
  matched substructure around the atom one neighbour at a time, following the
  attention weights of the graph neural network the tree was distilled from.
  Where the descent stops, a stored value is read off.

  References:
   - M. Lehner et al., <i>J. Chem. Inf. Model.</i> <b>2023</b>, 63, 6014
     (DOI 10.1021/acs.jcim.3c00800)
   - M. Lehner et al., <i>J. Chem. Phys.</i> <b>2024</b>, 161, 074103
     (DOI 10.1063/5.0218154)

  The tree itself is a data file, not part of the library: construct a DASHTree
  with the path to a \c .dash container, which \c python -m
  rdkit.Chem.DASHTree.convert writes from a legacy tree of the DASH-tree python
  package; the RDKit ships a pruned one as
  $RDBASE/Data/DASHTree/default_pruned.dash. The file is memory mapped, so
  construction costs microseconds and only the pages a query actually visits are
  ever read.

  Which properties a tree carries is a property of the file. Name the ones you
  want when you construct the object and only those are resolved:

  \code
    DASH::DASHTree tree("default.dash", {"result", "std"});
    std::vector<double> charges;
    tree.getPartialCharges(mol, charges);
  \endcode
*/

#include <RDGeneral/export.h>

#ifndef RD_DASHTREE_H
#define RD_DASHTREE_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
class Atom;
class ROMol;

namespace DASH {

class DASHTreeNode;

//! \brief one atom-feature class, laid out as the container stores it
/*!
  The 5-tuple DASH classifies atoms by. Which tuples exist is a property of
  the tree, not of the code: the container carries its table, one entry per
  branch, and an atom whose tuple is not in it cannot be assigned by that tree
  at all. The published trees have 122 classes.
*/
struct RDKIT_DASHTREE_EXPORT AtomFeature {
  std::uint8_t atomicNum;
  std::uint8_t degree;
  std::int8_t formalCharge;
  std::uint8_t conjugated;  //!< is any bond of the atom conjugated
  std::uint8_t numHs;       //!< total Hs, explicit neighbours included
  std::uint8_t reserved[3];
};
static_assert(sizeof(AtomFeature) == 8, "AtomFeature is an 8-byte file record");

//! \brief controls how a subgraph match descends
struct RDKIT_DASHTREE_EXPORT DASHParams {
  //! maximum number of tree levels to descend
  unsigned int maxDepth = 16;
  //! stop once the attention accumulated over the descent exceeds this
  double attentionThreshold = 10.0;
  //! stop once a single step contributes less attention than this
  double attentionIncrementThreshold = 0.0;
  //! where several atoms fit the next node, take the first in atom order, as
  //! the DASH-tree python package does, instead of following each until the
  //! tree tells them apart; the result then depends on the atom order
  bool legacyTieBreaking = false;
};

//! \brief how per-atom values are made to sum to a molecule-wide target
enum class Normalization {
  NONE,         //!< leave the raw tree values alone
  SYMMETRIC,    //!< spread the deficit equally over all atoms
  STD_WEIGHTED  //!< spread it proportionally to each atom's stored deviation
};

//! \brief controls a normalised assignment; partial charges are the usual case
struct RDKIT_DASHTREE_EXPORT NormalizationOptions {
  //! property column holding the values
  std::string valueProperty = "result";
  //! property column holding their standard deviations; needed by
  //! STD_WEIGHTED, the others only report it
  std::string stdProperty = "std";
  Normalization normalization = Normalization::STD_WEIGHTED;
  //! substituted for a stored deviation that is zero, negative or absent
  double defaultStdValue = 0.1;
  DASHParams params;
};

//! \brief a memory-mapped DASH tree
/*!
  Immutable once constructed and safe to query from many threads at once: every
  query reads the mapping and writes only into its own scratch space.

  Every query on a molecule first classifies all of its atoms and throws
  ValueErrorException if one is outside the tree's classes or the property
  cache is stale.
*/
class RDKIT_DASHTREE_EXPORT DASHTree {
 public:
  //! \brief maps a \c .dash container
  /*!
    \param filename path to the container
    \param properties the property columns to resolve; empty means every column
           in the file

    \throws BadFileException if the file cannot be opened or is not a container
    \throws ValueErrorException if it is corrupt or lacks a requested property
  */
  explicit DASHTree(
      const std::string &filename,
      const std::vector<std::string> &properties = std::vector<std::string>());
  ~DASHTree();

  DASHTree(const DASHTree &) = delete;
  DASHTree &operator=(const DASHTree &) = delete;
  //! nodes taken from the source stay valid; the moved-from tree may only be
  //! assigned to or destroyed
  DASHTree(DASHTree &&) noexcept;
  DASHTree &operator=(DASHTree &&) noexcept;

  //! number of atom-feature branches in the file
  unsigned int getNumBranches() const;
  //! total number of tree nodes in the file
  unsigned int getNumNodes() const;
  //! the property columns that were resolved
  std::vector<std::string> getPropertyNames() const;
  //! every property column the file holds, whether or not it was resolved
  std::vector<std::string> getAvailablePropertyNames() const;
  //! path the tree was mapped from
  const std::string &getFilename() const;

  // -------------------------------------------------------------------------
  //  the atom-feature classes, which the file carries
  // -------------------------------------------------------------------------
  //! \brief the class with the given branch index
  /*!
    \throws IndexErrorException if \p index is out of range
  */
  const AtomFeature &getAtomFeature(unsigned int index) const;
  //! branch index of an atom, or -1 if its tuple is none of the tree's classes
  int getAtomFeatureIndex(const Atom *atom) const;

  // -------------------------------------------------------------------------
  //  a single atom
  // -------------------------------------------------------------------------
  //! \brief matches one atom and reports the nodes it descended through
  /*!
    \param nodePath overwritten with the branch index followed by the node ids
                    of the descent, so it reads like the node paths the
                    DASH-tree python package reports

    \throws IndexErrorException if \p atomIdx is out of range
  */
  void getAtomNodePath(const ROMol &mol, unsigned int atomIdx,
                       std::vector<std::uint32_t> &nodePath,
                       const DASHParams &params = DASHParams()) const;

  //! \brief value of \p property for one atom
  /*!
    The deepest node of the match that carries a value wins; where the deepest
    nodes have none, the value of a shallower, more general substructure is
    used. NaN if no node on the path carries one. A float column marks a node
    without a value with NaN, an int32 column with its smallest value.
  */
  double getAtomProperty(const ROMol &mol, unsigned int atomIdx,
                         const std::string &property,
                         const DASHParams &params = DASHParams()) const;

  // -------------------------------------------------------------------------
  //  one molecule
  // -------------------------------------------------------------------------
  //! value of \p property for every atom of \p mol
  void getMolProperty(const ROMol &mol, const std::string &property,
                      std::vector<double> &res,
                      const DASHParams &params = DASHParams()) const;

  //! \brief values of \p options.valueProperty for every atom of \p mol,
  //! adjusted to sum to \p target
  /*!
    Partial charges are the published case: the target is the molecule's formal
    charge, which getPartialCharges() supplies. Any per-atom property with a
    molecule-wide total works the same way. \p options names the value column,
    the deviation column STD_WEIGHTED weights the adjustment by, and how to
    adjust; with NONE the target is ignored.

    \throws ValueErrorException if a column it needs was not resolved, if an
            atom's path carries no value, or if STD_WEIGHTED is asked for
            without a positive \c defaultStdValue
  */
  void getNormalizedMolProperty(
      const ROMol &mol, double target, std::vector<double> &res,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! \brief getNormalizedMolProperty(), also reporting the intermediates
  /*!
    \param res         the adjusted values
    \param rawValues   the values straight out of the tree, before normalisation
    \param stds        the deviations used, after \c defaultStdValue
                       substitution
    \param matchDepths how deep each atom's match went
  */
  void getNormalizedMolProperty(
      const ROMol &mol, double target, std::vector<double> &res,
      std::vector<double> &rawValues, std::vector<double> &stds,
      std::vector<unsigned int> &matchDepths,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! \brief partial charges for every atom of \p mol
  /*!
    getNormalizedMolProperty() with the molecule's formal charge as the target.
  */
  void getPartialCharges(
      const ROMol &mol, std::vector<double> &res,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! partial charges, also reporting the intermediates
  void getPartialCharges(
      const ROMol &mol, std::vector<double> &res,
      std::vector<double> &rawValues, std::vector<double> &stds,
      std::vector<unsigned int> &matchDepths,
      const NormalizationOptions &options = NormalizationOptions()) const;

  // -------------------------------------------------------------------------
  //  many molecules
  // -------------------------------------------------------------------------
  //! \brief value of \p property for every atom of every molecule
  /*!
    \param numThreads how many threads to use. If set to 0 this will use the
           maximum number of threads allowed on your system; a negative value
           leaves that many threads unused, so -1 keeps one core free. Ignored
           unless RDKit was built with RDK_BUILD_THREADSAFE_SSS.

    Molecules are handed out one at a time, so a batch of unequal molecules
    still balances across the threads. A null entry in \p mols yields an empty
    result vector rather than an error.
  */
  void getMolPropertyBatch(const std::vector<const ROMol *> &mols,
                           const std::string &property,
                           std::vector<std::vector<double>> &res,
                           const DASHParams &params = DASHParams(),
                           int numThreads = 1) const;

  //! \brief getNormalizedMolProperty() for every molecule, one target each
  /*!
    \throws ValueErrorException if \p targets and \p mols differ in length
  */
  void getNormalizedMolPropertyBatch(
      const std::vector<const ROMol *> &mols,
      const std::vector<double> &targets, std::vector<std::vector<double>> &res,
      const NormalizationOptions &options = NormalizationOptions(),
      int numThreads = 1) const;

  //! partial charges for every atom of every molecule
  void getPartialChargesBatch(
      const std::vector<const ROMol *> &mols,
      std::vector<std::vector<double>> &res,
      const NormalizationOptions &options = NormalizationOptions(),
      int numThreads = 1) const;

  // -------------------------------------------------------------------------
  //  the tree itself
  // -------------------------------------------------------------------------
  //! the root node of a branch
  DASHTreeNode getRoot(unsigned int branch) const;

  //! \brief a node by branch and id
  /*!
    Ids are the container's own, counted from the branch root -- the numbering
    getAtomNodePath() reports -- so \c getNode(path[0], path[k]) is the k-th
    node a match descended through.

    \throws IndexErrorException if the branch or the id is out of range
  */
  DASHTreeNode getNode(unsigned int branch, std::uint32_t nodeId) const;

  //! \brief the molecule atoms a match covers
  /*!
    \param atoms overwritten with the atom indices in the order the descent
           added them: position i is the atom that node i of the path matched,
           and the position a node key's conAtom refers to. A hydrogen is
           matched through its heavy neighbour, so for one the list starts at
           that neighbour, which node 1 matched, and the hydrogen itself is not
           in it.
  */
  void getMatchedSubstructure(const ROMol &mol, unsigned int atomIdx,
                              std::vector<unsigned int> &atoms,
                              const DASHParams &params = DASHParams()) const;

  //! defined in DASHTreeImpl.h, which is not installed
  class Impl;

 private:
  std::unique_ptr<Impl> d_impl;
};

//! \brief one node of a mapped tree
/*!
  A handle -- the tree and a node id -- that reads straight out of the mapped
  file, so it costs nothing to copy and one load to query. It borrows the tree
  and must not outlive it.

  A node describes itself relative to its parent: the atom-feature class of the
  atom it adds to the substructure, the position in the substructure of the
  atom that one attaches to, and the bond descriptor between them (1, 2, 3 for
  the bond order, 4 for a conjugated bond). A branch root, and the heavy-atom
  child of a hydrogen root, attach to nothing and report -1 for both.
*/
class RDKIT_DASHTREE_EXPORT DASHTreeNode {
 public:
  //! the branch the node belongs to
  unsigned int getBranch() const;
  //! the node's id within its branch, as getAtomNodePath() reports it
  std::uint32_t getId() const;
  //! atom-feature class of the atom this node adds
  int getAtomFeatureIndex() const;
  //! position of the substructure atom it attaches to, -1 for none
  int getConAtom() const;
  //! bond descriptor of that attachment, -1 for none
  int getConType() const;
  unsigned int getNumChildren() const;
  //! \throws IndexErrorException if \p i is not below getNumChildren()
  DASHTreeNode getChild(unsigned int i) const;
  //! the attention weight the tree assigned to this node
  float getAttention() const;
  //! does a descent at the default threshold stop here
  bool isStopNode() const;
  //! \brief the value of a property at this node, NaN if it carries none
  /*!
    \throws ValueErrorException if \p property was not resolved
  */
  double getValue(const std::string &property) const;

 private:
  DASHTreeNode(const DASHTree::Impl &impl, std::uint32_t absoluteId);
  std::uint16_t getKey() const;
  const DASHTree::Impl *dp_impl;
  std::uint32_t d_abs;  //!< the container-wide id
  friend class DASHTree;
};

}  // namespace DASH
}  // namespace RDKit

#endif
