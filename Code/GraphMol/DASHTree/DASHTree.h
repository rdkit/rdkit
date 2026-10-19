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
  matched substructure around the atom one neighbour at a time. The growth
  follows the attention weights of the graph neural network the tree was
  distilled from. Where the descent stops, DASH reads a stored value.

  References:
   - M. Lehner et al., <i>J. Chem. Inf. Model.</i> <b>2023</b>, 63, 6014
     (DOI 10.1021/acs.jcim.3c00800)
   - M. Lehner et al., <i>J. Chem. Phys.</i> <b>2024</b>, 161, 074103
     (DOI 10.1063/5.0218154)

  The tree itself is a data file, not part of the library: construct a DASHTree
  with the path to a \c .dash container, which \c python -m
  rdkit.Chem.DASHTree.convert writes from a legacy tree of the DASH-tree python
  package. The file is memory mapped, so construction costs microseconds and
  only the pages a query actually visits are ever read.

  Which properties a tree carries depends on the file. When you construct the
  object, name the ones you want. Only those are resolved:

  \code
    DASH::DASHTree tree("default.dash", {"result", "std"});
    std::vector<double> charges = tree.getPartialCharges(mol);
  \endcode
*/

#include <RDGeneral/export.h>

#ifndef RD_DASHTREE_H
#define RD_DASHTREE_H

#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace RDKit {
class Atom;
class ROMol;

namespace DASH {

class DASHTreeNode;

//! \brief one atom-feature class, in the record layout of the container
/*!
  DASH classifies atoms by this 5-tuple. Which tuples exist depends on the
  tree, not on the code. The container carries its table, one entry per
  branch. The tree cannot assign values to an atom whose tuple is not in the
  table. The published trees have 122 classes.
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
  //! the DASH-tree python package does. The result then depends on the atom
  //! order. By default, the descent follows each of them until the tree
  //! distinguishes them.
  bool legacyTieBreaking = false;
};

//! \brief how per-atom values are made to sum to a molecule-wide target
enum class Normalization {
  NONE,         //!< keep the raw tree values
  SYMMETRIC,    //!< spread the deficit equally over all atoms
  STD_WEIGHTED  //!< spread it proportionally to each atom's stored deviation
};

//! \brief controls a normalized assignment, usually of partial charges
struct RDKIT_DASHTREE_EXPORT NormalizationOptions {
  //! property column that holds the values
  std::string valueProperty = "result";
  //! property column that holds their standard deviations. STD_WEIGHTED needs
  //! it. The other modes only report it.
  std::string stdProperty = "std";
  Normalization normalization = Normalization::STD_WEIGHTED;
  //! substituted for a stored deviation that is zero, negative or absent
  double defaultStdValue = 0.1;
  DASHParams params;
};

//! \brief a memory-mapped DASH tree
/*!
  The tree is immutable once constructed and safe to query from many threads at
  once. Every query reads the mapping and writes only into its own scratch
  space.

  Every query on a molecule first classifies all of its atoms. It throws
  ValueErrorException if an atom is outside the tree's classes or carries
  hydrogens that are not explicit atoms, or if the property cache is stale. A
  query also throws IndexErrorException for an index out of range and
  ValueErrorException for a property that was not resolved.
*/
class RDKIT_DASHTREE_EXPORT DASHTree {
 public:
  //! \brief maps a \c .dash container
  /*!
    \param filename path to the container
    \param properties the property columns to resolve, or empty for every
           column in the file

    \throws BadFileException if the file cannot be opened or is not a container
    \throws ValueErrorException if it is corrupt or lacks a requested property
  */
  explicit DASHTree(
      const std::string &filename,
      const std::vector<std::string> &properties = std::vector<std::string>());
  ~DASHTree();

  DASHTree(const DASHTree &) = delete;
  DASHTree &operator=(const DASHTree &) = delete;
  //! nodes taken from the source stay valid. The moved-from tree may only be
  //! assigned to or destroyed.
  DASHTree(DASHTree &&) noexcept;
  DASHTree &operator=(DASHTree &&) noexcept;

  //! number of atom-feature branches in the file
  unsigned int getNumBranches() const;
  //! total number of tree nodes in the file
  unsigned int getNumNodes() const;
  //! the property columns that were resolved
  const std::vector<std::string> &getPropertyNames() const;
  //! every property column the file holds, whether or not it was resolved
  const std::vector<std::string> &getAvailablePropertyNames() const;
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
    \return the branch index followed by the node ids of the descent, so it
            reads like the node paths the DASH-tree python package reports

    \throws IndexErrorException if \p atomIdx is out of range
  */
  std::vector<std::uint32_t> getAtomNodePath(
      const ROMol &mol, unsigned int atomIdx,
      const DASHParams &params = DASHParams()) const;

  //! \brief value of \p property for one atom
  /*!
    The deepest node of the match that carries a value gives the result. Where
    the deepest nodes have none, a shallower, more general substructure gives
    it. The result is NaN if no node on the path carries a value. A float
    column marks a node without a value with NaN. An int32 column marks it with
    INT32_MIN.
  */
  double getAtomProperty(const ROMol &mol, unsigned int atomIdx,
                         const std::string &property,
                         const DASHParams &params = DASHParams()) const;

  // -------------------------------------------------------------------------
  //  one molecule
  // -------------------------------------------------------------------------
  //! value of \p property for every atom of \p mol
  std::vector<double> getMolProperty(
      const ROMol &mol, const std::string &property,
      const DASHParams &params = DASHParams()) const;

  //! \brief values of \p options.valueProperty for every atom of \p mol,
  //! adjusted to sum to \p target
  /*!
    Partial charges are the published case. Their target is the molecule's
    formal charge, which getPartialCharges() supplies. Any per-atom property
    with a molecule-wide total works the same way. \p options names the value
    column and the deviation column, and selects the normalization mode.
    STD_WEIGHTED weights the adjustment by the deviations. NONE ignores the
    target.

    \throws ValueErrorException if a column it needs was not resolved, if an
            atom's path carries no value, or if \p options selects STD_WEIGHTED
            without a positive \c defaultStdValue
  */
  std::vector<double> getNormalizedMolProperty(
      const ROMol &mol, double target,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! \brief getNormalizedMolProperty() and its intermediates
  /*!
    \return the adjusted values, the raw tree values before normalization, the
            deviations used after \c defaultStdValue substitution, and the
            match depth of each atom
  */
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
             std::vector<unsigned int>>
  getNormalizedMolPropertyDetails(
      const ROMol &mol, double target,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! \brief partial charges for every atom of \p mol
  /*!
    getNormalizedMolProperty() with the molecule's formal charge as the target.
  */
  std::vector<double> getPartialCharges(
      const ROMol &mol,
      const NormalizationOptions &options = NormalizationOptions()) const;

  //! partial charges and their intermediates, as
  //! getNormalizedMolPropertyDetails() reports them
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
             std::vector<unsigned int>>
  getPartialChargesDetails(
      const ROMol &mol,
      const NormalizationOptions &options = NormalizationOptions()) const;

  // -------------------------------------------------------------------------
  //  many molecules
  // -------------------------------------------------------------------------
  //! \brief value of \p property for every atom of every molecule
  /*!
    \param numThreads how many threads to use. If set to 0, this will use the
           maximum number of threads allowed on your system. A negative value
           leaves that many threads unused, so -1 keeps one core free. Ignored
           unless RDKit was built with RDK_BUILD_THREADSAFE_SSS.

    The threads take the molecules one at a time, so a batch of unequal
    molecules still balances across the threads. A null entry in \p mols
    yields an empty result vector rather than an error.
  */
  std::vector<std::vector<double>> getMolPropertyBatch(
      const std::vector<const ROMol *> &mols, const std::string &property,
      const DASHParams &params = DASHParams(), int numThreads = 1) const;

  //! \brief getNormalizedMolProperty() for every molecule, one target each
  /*!
    \throws ValueErrorException if \p targets and \p mols differ in length
  */
  std::vector<std::vector<double>> getNormalizedMolPropertyBatch(
      const std::vector<const ROMol *> &mols,
      const std::vector<double> &targets,
      const NormalizationOptions &options = NormalizationOptions(),
      int numThreads = 1) const;

  //! partial charges for every atom of every molecule
  std::vector<std::vector<double>> getPartialChargesBatch(
      const std::vector<const ROMol *> &mols,
      const NormalizationOptions &options = NormalizationOptions(),
      int numThreads = 1) const;

  // -------------------------------------------------------------------------
  //  the tree itself
  // -------------------------------------------------------------------------
  //! the root node of a branch
  DASHTreeNode getRoot(unsigned int branch) const;

  //! \brief a node by branch and id
  /*!
    Ids are the container's own, counted from the branch root. This is the
    numbering getAtomNodePath() reports, so \c getNode(path[0], path[k]) is the
    k-th node a match descended through.

    \throws IndexErrorException if the branch or the id is out of range
  */
  DASHTreeNode getNode(unsigned int branch, std::uint32_t nodeId) const;

  //! \brief the molecule atoms a match covers
  /*!
    \return the atom indices in the order the descent added them. With
            path = getAtomNodePath(mol, atomIdx, params), position i is the
            atom that getNode(path[0], path[i + 1]) matched, and the position
            a node key's conAtom refers to. A hydrogen is matched through its
            heavy neighbour. For a hydrogen, position i is the atom of
            path[i + 2], and the hydrogen itself is not in the list.
  */
  std::vector<unsigned int> getMatchedSubstructure(
      const ROMol &mol, unsigned int atomIdx,
      const DASHParams &params = DASHParams()) const;

  //! defined in DASHTreeImpl.h, which is not installed
  class Impl;

 private:
  std::unique_ptr<Impl> d_impl;
};

//! \brief one node of a mapped tree
/*!
  A node is a handle that holds a pointer to the tree and a node id. It reads
  straight from the mapped file, so it costs nothing to copy. It borrows the
  tree and must not outlive it.

  A node describes itself relative to its parent:
   - the atom-feature class of the atom it adds to the substructure
   - the position in the substructure of the atom that the added atom
     attaches to
   - the bond descriptor between them (1, 2, 3 for the bond order, 4 for a
     conjugated bond)

  A branch root, and the heavy-atom child of a hydrogen root, attach to nothing.
  For them, getConAtom() and getConType() return -1.
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
  double getAttention() const;
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
