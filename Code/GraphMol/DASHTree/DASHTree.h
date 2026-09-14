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

  For each atom, DASH walks a precomputed decision tree whose levels grow the
  matched substructure around the atom one neighbour at a time, following the
  attention weights of the graph neural network the tree was distilled from.
  Where the descent stops, a stored value is read off.

  References:
   - M. Lehner et al., <i>J. Chem. Inf. Model.</i> <b>2023</b>, 63, 6296
     (DOI 10.1021/acs.jcim.3c00800)
   - M. Lehner et al., <i>J. Chem. Phys.</i> <b>2024</b>, 161, 044113
     (DOI 10.1063/5.0218154)

  The tree itself is a data file, not part of the library: construct a DASHTree
  with the path to a \c .dash container produced by \c tools/dash_convert.py
  from a DASH-tree distribution. The file is memory mapped, so construction
  costs microseconds and only the pages a query actually visits are ever read.

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
class ROMol;

namespace DASH {

//! \brief controls how far a subgraph match descends
struct RDKIT_DASHTREE_EXPORT DASHParams {
  //! maximum number of tree levels to descend
  unsigned int maxDepth = 16;
  //! stop once the attention accumulated over the descent exceeds this
  double attentionThreshold = 10.0;
  //! stop once a single step contributes less attention than this
  double attentionIncrementThreshold = 0.0;
};

//! \brief how partial charges are made to sum to the molecule's formal charge
enum class ChargeNormalization {
  NONE,         //!< leave the raw tree values alone
  SYMMETRIC,    //!< spread the deficit equally over all atoms
  STD_WEIGHTED  //!< spread it proportionally to each atom's stored deviation
};

//! \brief controls partial-charge assignment
struct RDKIT_DASHTREE_EXPORT ChargeOptions {
  //! property column holding the charges
  std::string valueProperty = "result";
  //! property column holding their standard deviations
  std::string stdProperty = "std";
  ChargeNormalization normalization = ChargeNormalization::STD_WEIGHTED;
  //! substituted for a stored deviation that is zero or negative
  double defaultStdValue = 0.1;
  DASHParams params;
};

//! \brief a memory-mapped DASH tree
/*!
  Immutable once constructed and safe to query from many threads at once: every
  query reads the mapping and writes only into its own scratch space.
*/
class RDKIT_DASHTREE_EXPORT DASHTree {
 public:
  //! \brief maps a \c .dash container
  /*!
    \param filename path to the container
    \param properties the property columns to resolve; empty means every column
           in the file
    \param prefetch read through the mapped arrays once, so later queries do not
           pay page faults. Useful before a benchmark or a large batch; pointless
           for a handful of molecules, which is the case mapping is there for.

    \throws BadFileException if the file cannot be opened or is not a container
    \throws ValueErrorException if it is corrupt or lacks a requested property
  */
  DASHTree(const std::string &filename,
           const std::vector<std::string> &properties =
               std::vector<std::string>(),
           bool prefetch = false);
  ~DASHTree();

  DASHTree(const DASHTree &) = delete;
  DASHTree &operator=(const DASHTree &) = delete;

  //! number of atom-feature branches in the file
  unsigned int numBranches() const;
  //! total number of tree nodes in the file
  std::uint64_t numNodes() const;
  //! bytes of the file that were mapped
  std::uint64_t mappedSize() const;
  //! the property columns that were resolved
  std::vector<std::string> propertyNames() const;
  //! every property column the file holds, whether or not it was resolved
  std::vector<std::string> availablePropertyNames() const;
  //! was \p name resolved
  bool hasProperty(const std::string &name) const;
  //! path the tree was mapped from
  const std::string &filename() const;
  //! \brief does the file carry the map back to the source tree's node numbering
  /*!
    The container renumbers nodes breadth-first, which is what lets a node's
    children be a contiguous range. The original numbering is kept as an extra
    array -- never read unless asked for -- so node paths can be compared
    against the python implementation.
  */
  bool hasSourceNodeIds() const;

  // -------------------------------------------------------------------------
  //  a single atom
  // -------------------------------------------------------------------------
  //! \brief matches one atom and reports the nodes it descended through
  /*!
    \param nodePath      overwritten with the branch index followed by the node
                         ids of the descent, so it reads like the node paths the
                         DASH-tree python package reports
    \param sourceNodeIds report the source tree's node numbering instead of the
                         container's own; needs a file that carries it, see
                         hasSourceNodeIds()

    \throws ValueErrorException if the atom's feature tuple is not one of the
            DASH classes
  */
  void getAtomNodePath(const ROMol &mol, unsigned int atomIdx,
                       std::vector<std::uint32_t> &nodePath,
                       const DASHParams &params = DASHParams(),
                       bool sourceNodeIds = false) const;

  //! \brief value of \p property for one atom
  /*!
    The deepest node of the match that carries a value wins; where the deepest
    nodes have none, the value of a shallower, more general substructure is
    used. NaN if no node on the path carries one.
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

  //! \brief partial charges for every atom of \p mol
  /*!
    \param res overwritten with one charge per atom, normalised as
           \c options.normalization asks
  */
  void getPartialCharges(const ROMol &mol, std::vector<double> &res,
                         const ChargeOptions &options = ChargeOptions()) const;

  //! \brief partial charges, also reporting the intermediates
  /*!
    \param res         normalised charges
    \param rawValues   the values straight out of the tree, before normalisation
    \param stds        the deviations used, after \c defaultStdValue substitution
    \param matchDepths how deep each atom's match went
  */
  void getPartialCharges(const ROMol &mol, std::vector<double> &res,
                         std::vector<double> &rawValues,
                         std::vector<double> &stds,
                         std::vector<unsigned int> &matchDepths,
                         const ChargeOptions &options = ChargeOptions()) const;

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

  //! partial charges for every atom of every molecule
  void getPartialChargesBatch(const std::vector<const ROMol *> &mols,
                              std::vector<std::vector<double>> &res,
                              const ChargeOptions &options = ChargeOptions(),
                              int numThreads = 1) const;

  //! \brief the mapped container
  /*!
    Named here so the library's own translation units can refer to it; it is
    defined in DASHTreeImpl.h, which is not installed.
  */
  class Impl;

 private:
  std::unique_ptr<Impl> d_impl;
};

//! \brief the current container format version
/*!
  A file written by a newer converter than the library understands is rejected
  rather than misread.
*/
const std::uint32_t dashFormatVersion = 2;

}  // namespace DASH
}  // namespace RDKit

#endif
