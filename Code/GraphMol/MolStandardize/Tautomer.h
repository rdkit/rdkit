//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_TAUTOMER_H
#define RD_TAUTOMER_H

#include <boost/function.hpp>
#include <string>
#include <iterator>
#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogEntry.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogParams.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
class ROMol;
class RWMol;

namespace MolStandardize {

typedef RDCatalog::HierarchCatalog<TautomerCatalogEntry, TautomerCatalogParams,
                                   int>
    TautomerCatalog;

namespace TautomerScoringFunctions {
const std::string tautomerScoringVersion = "1.0.0";

RDKIT_MOLSTANDARDIZE_EXPORT int scoreRings(const ROMol &mol);
RDKIT_MOLSTANDARDIZE_EXPORT int scoreSubstructs(const ROMol &mol);
RDKIT_MOLSTANDARDIZE_EXPORT int scoreHeteroHs(const ROMol &mol);

inline int scoreTautomer(const ROMol &mol) {
  return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
}
}  // namespace TautomerScoringFunctions

enum class TautomerEnumeratorStatus {
  Completed = 0,
  MaxTautomersReached,
  MaxTransformsReached,
  Canceled
};

class Tautomer {
  friend class TautomerEnumerator;

 public:
  Tautomer() : d_numModifiedAtoms(0), d_numModifiedBonds(0) {}
  Tautomer(const ROMOL_SPTR &t, const ROMOL_SPTR &k, size_t a = 0, size_t b = 0)
      : tautomer(t),
        kekulized(k),
        d_numModifiedAtoms(a),
        d_numModifiedBonds(b) {}
  ROMOL_SPTR tautomer;
  ROMOL_SPTR kekulized;

 private:
  size_t d_numModifiedAtoms;
  size_t d_numModifiedBonds;
};

typedef std::map<std::string, Tautomer> SmilesTautomerMap;
typedef std::pair<std::string, Tautomer> SmilesTautomerPair;

//! Contains results of tautomer enumeration
class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumeratorResult {
  friend class TautomerEnumerator;

 public:
  class const_iterator {
   public:
    typedef ROMOL_SPTR value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const ROMol *pointer;
    typedef const ROMOL_SPTR &reference;
    typedef std::bidirectional_iterator_tag iterator_category;

    explicit const_iterator(const SmilesTautomerMap::const_iterator &it)
        : d_it(it) {}
    reference operator*() const { return d_it->second.tautomer; }
    pointer operator->() const { return d_it->second.tautomer.get(); }
    bool operator==(const const_iterator &other) const {
      return (d_it == other.d_it);
    }
    bool operator!=(const const_iterator &other) const {
      return !(*this == other);
    }
    const_iterator operator++(int) {
      const_iterator copy(d_it);
      operator++();
      return copy;
    }
    const_iterator &operator++() {
      ++d_it;
      return *this;
    }
    const_iterator operator--(int) {
      const_iterator copy(d_it);
      operator--();
      return copy;
    }
    const_iterator &operator--() {
      --d_it;
      return *this;
    }

   private:
    SmilesTautomerMap::const_iterator d_it;
  };
  TautomerEnumeratorResult() : d_status(TautomerEnumeratorStatus::Completed) {}
  TautomerEnumeratorResult(const TautomerEnumeratorResult &other)
      : d_tautomers(other.d_tautomers),
        d_status(other.d_status),
        d_modifiedAtoms(other.d_modifiedAtoms),
        d_modifiedBonds(other.d_modifiedBonds) {
    fillTautomersItVec();
  }
  const const_iterator begin() const {
    return const_iterator(d_tautomers.begin());
  }
  const const_iterator end() const { return const_iterator(d_tautomers.end()); }
  size_t size() const { return d_tautomers.size(); }
  bool empty() const { return d_tautomers.empty(); }
  const ROMOL_SPTR &at(size_t pos) const {
    PRECONDITION(pos < d_tautomers.size(), "index out of bounds");
    return d_tautomersItVec.at(pos)->second.tautomer;
  }
  const ROMOL_SPTR &operator[](size_t pos) const { return at(pos); }
  const boost::dynamic_bitset<> &modifiedAtoms() const {
    return d_modifiedAtoms;
  }
  const boost::dynamic_bitset<> &modifiedBonds() const {
    return d_modifiedBonds;
  }
  TautomerEnumeratorStatus status() const { return d_status; }
  std::vector<ROMOL_SPTR> tautomers() const {
    std::vector<ROMOL_SPTR> tautomerVec;
    tautomerVec.reserve(d_tautomers.size());
    std::transform(
        d_tautomers.begin(), d_tautomers.end(), std::back_inserter(tautomerVec),
        [](const SmilesTautomerPair &t) { return t.second.tautomer; });
    return tautomerVec;
  }
  std::vector<ROMOL_SPTR> operator()() const { return tautomers(); }
  std::vector<std::string> smiles() const {
    std::vector<std::string> smilesVec;
    smilesVec.reserve(d_tautomers.size());
    std::transform(d_tautomers.begin(), d_tautomers.end(),
                   std::back_inserter(smilesVec),
                   [](const SmilesTautomerPair &t) { return t.first; });
    return smilesVec;
  }
  const SmilesTautomerMap &smilesTautomerMap() const { return d_tautomers; }

 private:
  void fillTautomersItVec() {
    for (auto it = d_tautomers.begin(); it != d_tautomers.end(); ++it) {
      d_tautomersItVec.push_back(it);
    }
  }
  // the enumerated tautomers
  SmilesTautomerMap d_tautomers;
  // internal; vector of iterators into map items to enable random
  // access to map items by index
  std::vector<SmilesTautomerMap::const_iterator> d_tautomersItVec;
  // status of the enumeration: did it complete? did it hit a limit?
  // was it canceled?
  TautomerEnumeratorStatus d_status;
  // bit vector: flags atoms modified by the transforms
  boost::dynamic_bitset<> d_modifiedAtoms;
  // bit vector: flags bonds modified by the transforms
  boost::dynamic_bitset<> d_modifiedBonds;
};

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumeratorCallback {
 public:
  TautomerEnumeratorCallback() {}
  virtual ~TautomerEnumeratorCallback() {}
  virtual bool operator()(const ROMol &, const TautomerEnumeratorResult &) = 0;
};

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumerator {
 public:
  TautomerEnumerator(TautomerCatalog *tautCat)
      : dp_catalog(tautCat),
        d_maxTautomers(1000),
        d_maxTransforms(1000),
        d_removeSp3Stereo(true),
        d_removeBondStereo(true),
        d_removeIsotopicHs(true),
        d_reassignStereo(true) {}
  TautomerEnumerator(const CleanupParameters &params = CleanupParameters());
  TautomerEnumerator(const TautomerEnumerator &other)
      : dp_catalog(other.dp_catalog),
        d_callback(other.d_callback.get()),
        d_maxTautomers(other.d_maxTautomers),
        d_maxTransforms(other.d_maxTransforms),
        d_removeSp3Stereo(other.d_removeSp3Stereo),
        d_removeBondStereo(other.d_removeBondStereo),
        d_removeIsotopicHs(other.d_removeIsotopicHs),
        d_reassignStereo(other.d_reassignStereo) {}
  TautomerEnumerator &operator=(const TautomerEnumerator &other) {
    if (this == &other) return *this;
    dp_catalog = other.dp_catalog;
    d_callback.reset(other.d_callback.get());
    d_maxTautomers = other.d_maxTautomers;
    d_maxTransforms = other.d_maxTransforms;
    d_removeSp3Stereo = other.d_removeSp3Stereo;
    d_removeBondStereo = other.d_removeBondStereo;
    d_removeIsotopicHs = other.d_removeIsotopicHs;
    d_reassignStereo = other.d_reassignStereo;
    return *this;
  }
  //! \param maxTautomers maximum number of tautomers to be generated
  void setMaxTautomers(unsigned int maxTautomers) {
    d_maxTautomers = maxTautomers;
  }
  //! \return maximum number of tautomers to be generated
  unsigned int getMaxTautomers() { return d_maxTautomers; }
  /*! \param maxTransforms maximum number of transformations to be applied
      this limit is usually hit earlier than the maxTautomers limit
      and leads to a more linear scaling of CPU time with increasing
      number of tautomeric centers (see Sitzmann et al.)
   */
  void setMaxTransforms(unsigned int maxTransforms) {
    d_maxTransforms = maxTransforms;
  }
  //! \return maximum number of transformations to be applied
  unsigned int getMaxTransforms() { return d_maxTransforms; }
  /*! \param removeSp3Stereo; if set to true, stereochemistry information
      will be removed from sp3 atoms involved in tautomerism.
      This means that S-aminoacids will lose their stereochemistry after going
      through tautomer enumeration because of the amido-imidol tautomerism.
      This defaults to true in RDKit, false in the workflow described
      by Sitzmann et al.
   */
  void setRemoveSp3Stereo(bool removeSp3Stereo) {
    d_removeSp3Stereo = removeSp3Stereo;
  }
  /*! \return whether stereochemistry information will be removed from
      sp3 atoms involved in tautomerism
   */
  bool getRemoveSp3Stereo() { return d_removeSp3Stereo; }
  /*! \param removeBondStereo; if set to true, stereochemistry information
      will be removed from double bonds involved in tautomerism.
      This means that enols will lose their E/Z stereochemistry after going
      through tautomer enumeration because of the keto-enolic tautomerism.
      This defaults to true in RDKit and also in the workflow described
      by Sitzmann et al.
   */
  void setRemoveBondStereo(bool removeBondStereo) {
    d_removeBondStereo = removeBondStereo;
  }
  /*! \return whether stereochemistry information will be removed from
      double bonds involved in tautomerism
   */
  bool getRemoveBondStereo() { return d_removeBondStereo; }
  /*! \param removeIsotopicHs; if set to true, isotopic Hs
      will be removed from centers involved in tautomerism.
   */
  void setRemoveIsotopicHs(bool removeIsotopicHs) {
    d_removeIsotopicHs = removeIsotopicHs;
  }
  /*! \return whether isotpoic Hs will be removed from
      centers involved in tautomerism
   */
  bool getRemoveIsotopicHs() { return d_removeIsotopicHs; }
  /*! \param reassignStereo; if set to true, assignStereochemistry
      will be called on each tautomer generated by the enumerate() method.
      This defaults to true.
   */
  void setReassignStereo(bool reassignStereo) {
    d_reassignStereo = reassignStereo;
  }
  /*! \return whether assignStereochemistry will be called on each
      tautomer generated by the enumerate() method
   */
  bool getReassignStereo() { return d_reassignStereo; }
  /*! set this to an instance of a class derived from
      TautomerEnumeratorCallback where operator() is overridden.
      DO NOT delete the instance as ownership of the pointer is transferred
      to the TautomerEnumerator
   */
  void setCallback(TautomerEnumeratorCallback *callback) {
    d_callback.reset(callback);
  }
  /*! \return pointer to an instance of a class derived from
      TautomerEnumeratorCallback.
      DO NOT delete the instance as ownership of the pointer is transferred
      to the TautomerEnumerator
   */
  TautomerEnumeratorCallback *getCallback() const { return d_callback.get(); }

  //! returns a \c TautomerEnumeratorResult structure for the input molecule
  /*!
    The enumeration rules are inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

    \param mol: the molecule to be enumerated

    Note: the definitions used here are that the atoms modified during
    tautomerization are the atoms at the beginning and end of each tautomer
    transform (the H "donor" and H "acceptor" in the transform) and the bonds
    modified during transformation are any bonds whose order is changed during
    the tautomer transform (these are the bonds between the "donor" and the
    "acceptor")

  */
  TautomerEnumeratorResult enumerate(const ROMol &mol) const;

  //! Deprecated, please use the form returning a \c TautomerEnumeratorResult
  //! instead
  [
      [deprecated("please use the form returning a TautomerEnumeratorResult "
                  "instead")]] std::vector<ROMOL_SPTR>
  enumerate(const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms,
            boost::dynamic_bitset<> *modifiedBonds = nullptr) const;

  //! returns the canonical tautomer from a \c TautomerEnumeratorResult
  ROMol *pickCanonical(const TautomerEnumeratorResult &tautRes,
                       boost::function<int(const ROMol &mol)> scoreFunc =
                           TautomerScoringFunctions::scoreTautomer) const;

  //! returns the canonical tautomer from an iterable of possible tautomers
  // When Iterable is TautomerEnumeratorResult we use the other non-templated
  // overload for efficiency (TautomerEnumeratorResult already has SMILES so no
  // need to recompute them)
  template <class Iterable,
            typename std::enable_if<
                !std::is_same<Iterable, TautomerEnumeratorResult>::value,
                int>::type = 0>
  ROMol *pickCanonical(const Iterable &tautomers,
                       boost::function<int(const ROMol &mol)> scoreFunc =
                           TautomerScoringFunctions::scoreTautomer) const {
    ROMOL_SPTR bestMol;
    if (tautomers.size() == 1) {
      bestMol = *tautomers.begin();
    } else {
      // Calculate score for each tautomer
      int bestScore = std::numeric_limits<int>::min();
      std::string bestSmiles = "";
      for (const auto &t : tautomers) {
        auto score = scoreFunc(*t);
#ifdef VERBOSE_ENUMERATION
        std::cerr << "  " << MolToSmiles(*t) << " " << score << std::endl;
#endif
        if (score > bestScore) {
          bestScore = score;
          bestSmiles = MolToSmiles(*t);
          bestMol = t;
        } else if (score == bestScore) {
          auto smiles = MolToSmiles(*t);
          if (smiles < bestSmiles) {
            bestSmiles = smiles;
            bestMol = t;
          }
        }
      }
    }
    ROMol *res = new ROMol(*bestMol);
    static const bool cleanIt = true;
    static const bool force = true;
    MolOps::assignStereochemistry(*res, cleanIt, force);

    return res;
  }

  //! returns the canonical tautomer for a molecule
  /*!
    Note that the canonical tautomer is very likely not the most stable tautomer
    for any given conditions. The default scoring rules are designed to produce
    "reasonable" tautomers, but the primary concern is that the results are
    canonical: you always get the same canonical tautomer for a molecule
    regardless of what the input tautomer or atom ordering were.

    The default scoring scheme is inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

  */
  ROMol *canonicalize(const ROMol &mol,
                      boost::function<int(const ROMol &mol)> scoreFunc =
                          TautomerScoringFunctions::scoreTautomer) const;

 private:
  bool setTautomerStereoAndIsoHs(const ROMol &mol, ROMol &taut,
                                 const TautomerEnumeratorResult &res) const;
  std::shared_ptr<TautomerCatalog> dp_catalog;
  std::unique_ptr<TautomerEnumeratorCallback> d_callback;
  unsigned int d_maxTautomers;
  unsigned int d_maxTransforms;
  bool d_removeSp3Stereo;
  bool d_removeBondStereo;
  bool d_removeIsotopicHs;
  bool d_reassignStereo;
};  // TautomerEnumerator class

}  // namespace MolStandardize
}  // namespace RDKit

#endif
