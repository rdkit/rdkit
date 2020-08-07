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
#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogEntry.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogParams.h>
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

enum TautomerEnumeratorStatus {
  Completed = 0,
  MaxTautomersReached,
  MaxTransformsReached,
  Canceled
};

//! Contains results of tautomer enumeration
struct RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumeratorResult {
  // the enumerated tautomers
  std::vector<ROMOL_SPTR> tautomers;
  // status of the enumeration: did it complete? did it hit a limit?
  // was it canceled?
  TautomerEnumeratorStatus status{Completed};
  // bit vector: flags atoms modified by the transforms
  boost::dynamic_bitset<> modifiedAtoms;
  // bit vector: flags bonds modified by the transforms
  boost::dynamic_bitset<> modifiedBonds;
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
        d_removeBondStereo(true) {}
  TautomerEnumerator(const CleanupParameters &params = CleanupParameters());
  TautomerEnumerator(const TautomerEnumerator &other)
      : dp_catalog(other.dp_catalog),
        d_callback(other.d_callback.get()),
        d_maxTautomers(other.d_maxTautomers),
        d_maxTransforms(other.d_maxTransforms),
        d_removeSp3Stereo(other.d_removeSp3Stereo),
        d_removeBondStereo(other.d_removeBondStereo) {}
  TautomerEnumerator &operator=(const TautomerEnumerator &other) {
    if (this == &other) return *this;
    dp_catalog = other.dp_catalog;
    d_callback.reset(other.d_callback.get());
    d_maxTautomers = other.d_maxTautomers;
    d_maxTransforms = other.d_maxTransforms;
    d_removeSp3Stereo = other.d_removeSp3Stereo;
    d_removeBondStereo = other.d_removeBondStereo;
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

  //! returns a TautomerEnumeratorResult structure for the input molecule
  /*!
    The enumeration rules are inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

    \param mol: the molecule to be enumerated
    \param reassignStereo: whether AssignStereochemistry should be called
    on each generated tautomer (not required for canonicalization)
    \param modifiedBonds: if provided this is used to return which bonds are
    modified during the tautomerization

    Note: the definitions used here are that the atoms modified during
    tautomerization are the atoms at the beginning and end of each tautomer
    transform (the H "donor" and H "acceptor" in the transform) and the bonds
    modified during transformation are any bonds whose order is changed during
    the tautomer transform (these are the bonds between the "donor" and the
    "acceptor")

  */
  TautomerEnumeratorResult enumerate(const ROMol &mol,
                                     bool reassignStereo = true) const;

  //! returns the canonical tautomer from a set of possible tautomers
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
  ROMol *pickCanonical(std::vector<ROMOL_SPTR> &tautomers,
                       boost::function<int(const ROMol &mol)> scoreFunc =
                           TautomerScoringFunctions::scoreTautomer) const;

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
                          TautomerScoringFunctions::scoreTautomer) const {
    auto res = enumerate(mol, false);
    if (res.tautomers.empty()) {
      BOOST_LOG(rdWarningLog)
          << "no tautomers found, returning input molecule" << std::endl;
      return new ROMol(mol);
    }
    return pickCanonical(res.tautomers, scoreFunc);
  };

 private:
  bool setTautomerStereo(const ROMol &mol, ROMol &taut,
                         const TautomerEnumeratorResult &res,
                         bool reassignStereo) const;
  std::shared_ptr<TautomerCatalog> dp_catalog;
  std::unique_ptr<TautomerEnumeratorCallback> d_callback;
  unsigned int d_maxTautomers;
  unsigned int d_maxTransforms;
  bool d_removeSp3Stereo;
  bool d_removeBondStereo;
};  // TautomerEnumerator class

}  // namespace MolStandardize
}  // namespace RDKit

#endif
