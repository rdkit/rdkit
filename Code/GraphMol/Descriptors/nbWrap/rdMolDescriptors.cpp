//
//  Copyright (C) 2007-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/shared_ptr.h>
#include <RDBoost/Wrap_nb.h>
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>

#include <RDGeneral/RDLog.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomFeat.h>
#include <GraphMol/Descriptors/OxidationNumbers.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Descriptors/DCLV.h>
#include <DataStructs/BitVects.h>

#include <GraphMol/Descriptors/USRDescriptor.h>

#ifdef RDK_HAS_EIGEN3
#include <GraphMol/Descriptors/BCUT.h>
#endif

#ifdef RDK_BUILD_DESCRIPTORS3D
#include <GraphMol/Descriptors/MolDescriptors3D.h>
#endif

#include <RDGeneral/Exceptions.h>

#include <boost/dynamic_bitset.hpp>
#include <optional>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

struct AtomPairsParameters {};

namespace {

std::vector<unsigned int> atomPairTypes(
    RDKit::AtomPairs::atomNumberTypes,
    RDKit::AtomPairs::atomNumberTypes +
        sizeof(RDKit::AtomPairs::atomNumberTypes) / sizeof(unsigned int));

std::pair<std::vector<double>, double> computeASAContribs(
    const RDKit::ROMol &mol, bool includeHs = true, bool force = false) {
  std::vector<double> contribs(mol.getNumAtoms());
  double hContrib = 0.0;
  RDKit::Descriptors::getLabuteAtomContribs(mol, contribs, hContrib, includeHs,
                                            force);
  return std::make_pair(contribs, hContrib);
}

std::vector<double> computeTPSAContribs(const RDKit::ROMol &mol, bool force,
                                        bool includeSandP) {
  std::vector<double> contribs(mol.getNumAtoms());
  RDKit::Descriptors::getTPSAAtomContribs(mol, contribs, force, includeSandP);
  return contribs;
}

std::vector<std::pair<double, double>> computeCrippenContribs(
    const RDKit::ROMol &mol, bool force = false,
    nb::object atomTypes = nb::none(),
    nb::object atomTypeLabels = nb::none()) {
  std::optional<std::vector<unsigned int>> tAtomTypes;
  std::optional<std::vector<std::string>> tAtomTypeLabels;

  if (!atomTypes.is_none()) {
    auto atList = nb::cast<nb::list>(atomTypes);
    auto n = nb::len(atList);
    if (n != 0) {
      if (n != mol.getNumAtoms()) {
        throw ValueErrorException(
            "if atomTypes vector is provided, it must be as long as the number "
            "of atoms");
      }
      tAtomTypes.emplace(mol.getNumAtoms(), 0u);
    }
  }
  if (!atomTypeLabels.is_none()) {
    auto atLabelList = nb::cast<nb::list>(atomTypeLabels);
    auto n = nb::len(atLabelList);
    if (n != 0) {
      if (n != mol.getNumAtoms()) {
        throw ValueErrorException(
            "if atomTypeLabels vector is provided, it must be as long as the "
            "number of atoms");
      }
      tAtomTypeLabels.emplace(mol.getNumAtoms(), "");
    }
  }

  std::vector<double> logpContribs(mol.getNumAtoms());
  std::vector<double> mrContribs(mol.getNumAtoms());

  RDKit::Descriptors::getCrippenAtomContribs(
      mol, logpContribs, mrContribs, force,
      tAtomTypes ? &(*tAtomTypes) : nullptr,
      tAtomTypeLabels ? &(*tAtomTypeLabels) : nullptr);

  if (tAtomTypes) {
    auto atList = nb::cast<nb::list>(atomTypes);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      atList[i] = nb::int_((*tAtomTypes)[i]);
    }
  }
  if (tAtomTypeLabels) {
    auto atLabelList = nb::cast<nb::list>(atomTypeLabels);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      atLabelList[i] = nb::str((*tAtomTypeLabels)[i].c_str());
    }
  }

  std::vector<std::pair<double, double>> pycontribs;
  pycontribs.reserve(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    pycontribs.emplace_back(logpContribs[i], mrContribs[i]);
  }
  return pycontribs;
}

std::pair<double, double> calcCrippenDescriptors(const RDKit::ROMol &mol,
                                                 bool includeHs = true,
                                                 bool force = false) {
  double logp, mr;
  RDKit::Descriptors::calcCrippenDescriptors(mol, logp, mr, includeHs, force);
  return std::make_pair(logp, mr);
}

#ifdef RDK_BUILD_DESCRIPTORS3D

std::vector<std::vector<double>> calcCoulombMat(const RDKit::ROMol &mol,
                                                int confId) {
  std::vector<std::vector<double>> results;
  RDKit::Descriptors::CoulombMat(mol, results, confId);
  return results;
}

std::vector<double> calcEEMcharges(RDKit::ROMol &mol, int confId) {
  std::vector<double> res;
  RDKit::Descriptors::EEM(mol, res, confId);
  return res;
}

std::vector<double> calcWHIMs(const RDKit::ROMol &mol, int confId,
                               double thresh,
                               const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::WHIM(mol, res, confId, thresh, CustomAtomProperty);
  return res;
}

std::vector<double> calcGETAWAYs(const RDKit::ROMol &mol, int confId,
                                  double precision,
                                  const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::GETAWAY(mol, res, confId, precision, CustomAtomProperty);
  return res;
}

std::vector<double> calcRDFs(const RDKit::ROMol &mol, int confId,
                              const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::RDF(mol, res, confId, CustomAtomProperty);
  return res;
}

std::vector<double> calcMORSEs(const RDKit::ROMol &mol, int confId,
                                const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::MORSE(mol, res, confId, CustomAtomProperty);
  return res;
}

std::vector<double> calcAUTOCORR3Ds(const RDKit::ROMol &mol, int confId,
                                     const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::AUTOCORR3D(mol, res, confId, CustomAtomProperty);
  return res;
}

std::vector<double> calcAUTOCORR2Ds(const RDKit::ROMol &mol,
                                     const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::AUTOCORR2D(mol, res, CustomAtomProperty);
  return res;
}

#endif

RDKit::SparseIntVect<std::int32_t> *GetAtomPairFingerprint(
    const RDKit::ROMol &mol, unsigned int minLength, unsigned int maxLength,
    nb::object fromAtoms, nb::object ignoreAtoms, nb::object atomInvariants,
    bool includeChirality, bool use2D, int confId) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  RDKit::SparseIntVect<std::int32_t> *res;
  res = RDKit::AtomPairs::getAtomPairFingerprint(
      mol, minLength, maxLength, fvect.get(), ivect.get(), invvect.get(),
      includeChirality, use2D, confId);
  return res;
}

RDKit::SparseIntVect<std::int32_t> *GetHashedAtomPairFingerprint(
    const RDKit::ROMol &mol, unsigned int nBits, unsigned int minLength,
    unsigned int maxLength, nb::object fromAtoms, nb::object ignoreAtoms,
    nb::object atomInvariants, bool includeChirality, bool use2D, int confId) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  RDKit::SparseIntVect<std::int32_t> *res;
  res = RDKit::AtomPairs::getHashedAtomPairFingerprint(
      mol, nBits, minLength, maxLength, fvect.get(), ivect.get(), invvect.get(),
      includeChirality, use2D, confId);
  return res;
}

RDKit::SparseIntVect<boost::int64_t> *GetTopologicalTorsionFingerprint(
    const RDKit::ROMol &mol, unsigned int targetSize, nb::object fromAtoms,
    nb::object ignoreAtoms, nb::object atomInvariants, bool includeChirality) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  if (targetSize * RDKit::AtomPairs::codeSize > 64) {
    std::ostringstream errout;
    errout << "Maximum supported topological torsion path length is "
           << 64 / RDKit::AtomPairs::codeSize << std::endl;
    throw ValueErrorException(errout.str());
  }

  RDKit::SparseIntVect<boost::int64_t> *res;
  res = RDKit::AtomPairs::getTopologicalTorsionFingerprint(
      mol, targetSize, fvect.get(), ivect.get(), invvect.get(),
      includeChirality);
  return res;
}

RDKit::SparseIntVect<boost::int64_t> *GetHashedTopologicalTorsionFingerprint(
    const RDKit::ROMol &mol, unsigned int nBits, unsigned int targetSize,
    nb::object fromAtoms, nb::object ignoreAtoms, nb::object atomInvariants,
    bool includeChirality) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  RDKit::SparseIntVect<boost::int64_t> *res;
  res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(
      mol, nBits, targetSize, fvect.get(), ivect.get(), invvect.get(),
      includeChirality);
  return res;
}

ExplicitBitVect *GetHashedTopologicalTorsionFingerprintAsBitVect(
    const RDKit::ROMol &mol, unsigned int nBits, unsigned int targetSize,
    nb::object fromAtoms, nb::object ignoreAtoms, nb::object atomInvariants,
    unsigned int nBitsPerEntry, bool includeChirality) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  ExplicitBitVect *res;
  res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(
      mol, nBits, targetSize, fvect.get(), ivect.get(), invvect.get(),
      nBitsPerEntry, includeChirality);
  return res;
}

ExplicitBitVect *GetHashedAtomPairFingerprintAsBitVect(
    const RDKit::ROMol &mol, unsigned int nBits, unsigned int minLength,
    unsigned int maxLength, nb::object fromAtoms, nb::object ignoreAtoms,
    nb::object atomInvariants, unsigned int nBitsPerEntry, bool includeChirality,
    bool use2D, int confId) {
  std::unique_ptr<std::vector<std::uint32_t>> fvect =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> ivect =
      pythonObjectToVect(ignoreAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<std::uint32_t>> invvect = pythonObjectToVect(
      atomInvariants,
      static_cast<unsigned int>(1 << RDKit::AtomPairs::codeSize));
  ExplicitBitVect *res;
  res = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(
      mol, nBits, minLength, maxLength, fvect.get(), ivect.get(), invvect.get(),
      nBitsPerEntry, includeChirality, use2D, confId);
  return res;
}

double kappaHelper(double (*fn)(const RDKit::ROMol &, std::vector<double> *),
                   const RDKit::ROMol &mol, nb::object atomContribs) {
  std::optional<std::vector<double>> lContribs;
  if (!atomContribs.is_none()) {
    auto acl = nb::cast<nb::list>(atomContribs);
    if (nb::len(acl) != mol.getNumAtoms()) {
      throw ValueErrorException(
          "length of atomContribs list != number of atoms");
    }
    lContribs.emplace(mol.getNumAtoms());
  }
  double res = fn(mol, lContribs ? &(*lContribs) : nullptr);
  if (lContribs) {
    auto acl = nb::cast<nb::list>(atomContribs);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      acl[i] = nb::float_((*lContribs)[i]);
    }
  }
  return res;
}

double hkAlphaHelper(const RDKit::ROMol &mol, nb::object atomContribs) {
  return kappaHelper(RDKit::Descriptors::calcHallKierAlpha, mol, atomContribs);
}

[[deprecated(
    "please use MorganGenerator")]] RDKit::SparseIntVect<std::uint32_t> *
MorganFingerprintHelper(const RDKit::ROMol &mol, unsigned int radius, int nBits,
                        nb::object invariants, nb::object fromAtoms,
                        bool useChirality, bool useBondTypes, bool useFeatures,
                        bool useCounts, nb::object bitInfo,
                        bool includeRedundantEnvironments) {
  RDLog::deprecationWarning("please use MorganGenerator");
  std::vector<boost::uint32_t> *invars = nullptr;
  bool haveInvars = false;
  if (!invariants.is_none()) {
    unsigned int nInvar = nb::len(invariants);
    if (nInvar) {
      haveInvars = true;
      if (nInvar != mol.getNumAtoms()) {
        throw ValueErrorException(
            "length of invariant vector != number of atoms");
      }
      invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        (*invars)[i] = nb::cast<std::uint32_t>(invariants[i]);
      }
    }
  }
  if (!haveInvars && useFeatures) {
    invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol, *invars);
  }
  std::vector<std::uint32_t> *froms = nullptr;
  if (!fromAtoms.is_none()) {
    unsigned int nFrom = nb::len(fromAtoms);
    if (nFrom) {
      froms = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < nFrom; ++i) {
        froms->push_back(nb::cast<std::uint32_t>(fromAtoms[i]));
      }
    }
  }
  RDKit::MorganFingerprints::BitInfoMap *bitInfoMap = nullptr;
  if (!bitInfo.is_none()) {
    bitInfoMap = new RDKit::MorganFingerprints::BitInfoMap();
  }
  RDKit::SparseIntVect<std::uint32_t> *res;
  if (nBits < 0) {
    res = RDKit::MorganFingerprints::getFingerprint(
        mol, static_cast<unsigned int>(radius), invars, froms, useChirality,
        useBondTypes, useCounts, false, bitInfoMap,
        includeRedundantEnvironments);
  } else {
    res = RDKit::MorganFingerprints::getHashedFingerprint(
        mol, static_cast<unsigned int>(radius),
        static_cast<unsigned int>(nBits), invars, froms, useChirality,
        useBondTypes, false, bitInfoMap, includeRedundantEnvironments);
  }
  if (bitInfoMap) {
    bitInfo.attr("clear")();
    nb::dict biDict = nb::cast<nb::dict>(bitInfo);
    for (RDKit::MorganFingerprints::BitInfoMap::const_iterator iter =
             bitInfoMap->begin();
         iter != bitInfoMap->end(); ++iter) {
      const std::vector<std::pair<std::uint32_t, std::uint32_t>> &v =
          iter->second;
      nb::list localL;
      for (const auto &vIt : v) {
        localL.append(nb::make_tuple(vIt.first, vIt.second));
      }
      biDict[nb::int_(iter->first)] = nb::tuple(localL);
    }
    delete bitInfoMap;
  }
  delete invars;
  delete froms;
  return res;
}

#ifdef RDK_HAS_EIGEN3
std::vector<double> BCUT(const RDKit::ROMol &mol) {
  return RDKit::Descriptors::BCUT2D(mol);
}

std::pair<double, double> BCUT2D_list(const RDKit::ROMol &m,
                                      nb::list atomprops) {
  std::vector<double> dvec;
  for (size_t i = 0; i < nb::len(atomprops); ++i) {
    dvec.push_back(nb::cast<double>(atomprops[i]));
  }
  return RDKit::Descriptors::BCUT2D(m, dvec);
}

std::pair<double, double> BCUT2D_tuple(const RDKit::ROMol &m,
                                       nb::tuple atomprops) {
  std::vector<double> dvec;
  for (size_t i = 0; i < nb::len(atomprops); ++i) {
    dvec.push_back(nb::cast<double>(atomprops[i]));
  }
  return RDKit::Descriptors::BCUT2D(m, dvec);
}
#endif

RDKit::SparseIntVect<std::uint32_t> *GetMorganFingerprint(
    const RDKit::ROMol &mol, unsigned int radius, nb::object invariants,
    nb::object fromAtoms, bool useChirality, bool useBondTypes,
    bool useFeatures, bool useCounts, nb::object bitInfo,
    bool includeRedundantEnvironments) {
  return MorganFingerprintHelper(
      mol, radius, -1, invariants, fromAtoms, useChirality, useBondTypes,
      useFeatures, useCounts, bitInfo, includeRedundantEnvironments);
}

RDKit::SparseIntVect<std::uint32_t> *GetHashedMorganFingerprint(
    const RDKit::ROMol &mol, unsigned int radius, unsigned int nBits,
    nb::object invariants, nb::object fromAtoms, bool useChirality,
    bool useBondTypes, bool useFeatures, nb::object bitInfo,
    bool includeRedundantEnvironments) {
  return MorganFingerprintHelper(mol, radius, nBits, invariants, fromAtoms,
                                 useChirality, useBondTypes, useFeatures, true,
                                 bitInfo, includeRedundantEnvironments);
}

[[deprecated("please use MorganGenerator")]] ExplicitBitVect *
GetMorganFingerprintBV(const RDKit::ROMol &mol, unsigned int radius,
                       unsigned int nBits, nb::object invariants,
                       nb::object fromAtoms, bool useChirality,
                       bool useBondTypes, bool useFeatures, nb::object bitInfo,
                       bool includeRedundantEnvironments) {
  RDLog::deprecationWarning("please use MorganGenerator");
  std::vector<boost::uint32_t> *invars = nullptr;
  bool haveInvars = false;
  if (!invariants.is_none()) {
    unsigned int nInvar = nb::len(invariants);
    if (nInvar) {
      haveInvars = true;
      if (nInvar != mol.getNumAtoms()) {
        throw ValueErrorException(
            "length of invariant vector != number of atoms");
      }
      invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        (*invars)[i] = nb::cast<std::uint32_t>(invariants[i]);
      }
    }
  }
  if (!haveInvars && useFeatures) {
    invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol, *invars);
  }

  std::unique_ptr<std::vector<std::uint32_t>> froms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  RDKit::MorganFingerprints::BitInfoMap *bitInfoMap = nullptr;
  if (!bitInfo.is_none()) {
    bitInfoMap = new RDKit::MorganFingerprints::BitInfoMap();
  }
  ExplicitBitVect *res;
  res = RDKit::MorganFingerprints::getFingerprintAsBitVect(
      mol, radius, nBits, invars, froms.get(), useChirality, useBondTypes,
      false, bitInfoMap, includeRedundantEnvironments);
  if (bitInfoMap) {
    bitInfo.attr("clear")();
    nb::dict biDict = nb::cast<nb::dict>(bitInfo);
    for (RDKit::MorganFingerprints::BitInfoMap::const_iterator iter =
             bitInfoMap->begin();
         iter != bitInfoMap->end(); ++iter) {
      const std::vector<std::pair<std::uint32_t, std::uint32_t>> &v =
          iter->second;
      nb::list localL;
      for (const auto &vIt : v) {
        localL.append(nb::make_tuple(vIt.first, vIt.second));
      }
      biDict[nb::int_(iter->first)] = nb::tuple(localL);
    }
    delete bitInfoMap;
  }
  delete invars;
  return res;
}

std::vector<double> GetAtomFeatures(const RDKit::ROMol &mol, int atomid,
                                    bool addchiral) {
  std::vector<double> res;
  RDKit::Descriptors::AtomFeatVect(mol, res, atomid, addchiral);
  return res;
}

std::vector<unsigned int> GetConnectivityInvariants(
    const RDKit::ROMol &mol, bool includeRingMembership) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  RDKit::MorganFingerprints::getConnectivityInvariants(mol, invars,
                                                       includeRingMembership);
  return invars;
}

std::vector<unsigned int> GetFeatureInvariants(const RDKit::ROMol &mol) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  RDKit::MorganFingerprints::getFeatureInvariants(mol, invars);
  return invars;
}

std::vector<double> GetUSR(const RDKit::ROMol &mol, int confId) {
  if (mol.getNumConformers() == 0) {
    throw ValueErrorException("no conformers");
  }
  if (mol.getNumAtoms() < 3) {
    throw ValueErrorException("too few atoms (minimum three)");
  }
  std::vector<double> descriptor(12);
  RDKit::Descriptors::USR(mol, descriptor, confId);
  return descriptor;
}

nb::list GetUSRDistributions(nb::object coords, nb::object points) {
  unsigned int numCoords = nb::len(coords);
  if (numCoords == 0) {
    throw ValueErrorException("no coordinates");
  }
  RDGeom::Point3DConstPtrVect c(numCoords);
  for (unsigned int i = 0; i < numCoords; ++i) {
    auto *pt = new RDGeom::Point3D;
    *pt = nb::cast<RDGeom::Point3D>(coords[i]);
    c[i] = pt;
  }
  std::vector<RDGeom::Point3D> pts(4);
  std::vector<std::vector<double>> distances(4);
  RDKit::Descriptors::calcUSRDistributions(c, distances, pts);
  if (!points.is_none()) {
    nb::list tmpPts = nb::cast<nb::list>(points);
    for (const auto &p : pts) {
      tmpPts.append(p);
    }
  }
  nb::list pyDist;
  for (const auto &dist : distances) {
    nb::list pytmp;
    for (const auto d : dist) {
      pytmp.append(d);
    }
    pyDist.append(pytmp);
  }
  for (const auto *pt : c) {
    delete pt;
  }
  return pyDist;
}

nb::list GetUSRDistributionsFromPoints(nb::object coords, nb::object points) {
  unsigned int numCoords = nb::len(coords);
  unsigned int numPts = nb::len(points);
  if (numCoords == 0) {
    throw ValueErrorException("no coordinates");
  }
  RDGeom::Point3DConstPtrVect c(numCoords);
  for (unsigned int i = 0; i < numCoords; ++i) {
    c[i] = nb::cast<RDGeom::Point3D *>(coords[i]);
  }
  std::vector<RDGeom::Point3D> p(numPts);
  if (numPts == 0) {
    throw ValueErrorException("no points");
  }
  for (unsigned int i = 0; i < numPts; ++i) {
    p[i] = nb::cast<RDGeom::Point3D>(points[i]);
  }
  std::vector<std::vector<double>> distances(numPts);
  RDKit::Descriptors::calcUSRDistributionsFromPoints(c, p, distances);
  nb::list pyDist;
  for (const auto &dist : distances) {
    nb::list pytmp;
    for (const auto d : dist) {
      pytmp.append(d);
    }
    pyDist.append(pytmp);
  }
  return pyDist;
}

std::vector<double> GetUSRFromDistributions(nb::object distances) {
  unsigned int numDist = nb::len(distances);
  if (numDist == 0) {
    throw ValueErrorException("no distances");
  }
  std::vector<std::vector<double>> dist(numDist);
  for (unsigned int i = 0; i < numDist; ++i) {
    nb::object inner = distances[i];
    unsigned int numPts = nb::len(inner);
    if (numPts == 0) {
      throw ValueErrorException("distances missing");
    }
    std::vector<double> tmpDist(numPts);
    for (unsigned int j = 0; j < numPts; ++j) {
      tmpDist[j] = nb::cast<double>(inner[j]);
    }
    dist[i] = tmpDist;
  }
  std::vector<double> descriptor(12);
  RDKit::Descriptors::calcUSRFromDistributions(dist, descriptor);
  return descriptor;
}

double GetUSRScore(nb::object descriptor1, nb::object descriptor2,
                   nb::object weights) {
  unsigned int numElements = nb::len(descriptor1);
  if (numElements != nb::len(descriptor2)) {
    throw ValueErrorException("descriptors must have the same length");
  }
  unsigned int numWeights = numElements / 12;
  unsigned int numPyWeights = nb::len(weights);
  std::vector<double> w(numWeights, 1.0);
  if ((numPyWeights > 0) && (numPyWeights != numWeights)) {
    throw ValueErrorException("number of weights is not correct");
  } else if (numPyWeights == numWeights) {
    for (unsigned int i = 0; i < numWeights; ++i) {
      w[i] = nb::cast<double>(weights[i]);
    }
  }
  std::vector<double> d1(numElements);
  std::vector<double> d2(numElements);
  for (unsigned int i = 0; i < numElements; ++i) {
    d1[i] = nb::cast<double>(descriptor1[i]);
    d2[i] = nb::cast<double>(descriptor2[i]);
  }
  return RDKit::Descriptors::calcUSRScore(d1, d2, w);
}

std::vector<double> GetUSRCAT(const RDKit::ROMol &mol,
                               nb::object atomSelections, int confId) {
  if (mol.getNumConformers() == 0) {
    throw ValueErrorException("no conformers");
  }
  if (mol.getNumAtoms() < 3) {
    throw ValueErrorException("too few atoms (minimum three)");
  }

  std::vector<std::vector<unsigned int>> atomIds;
  unsigned int sizeDescriptor = 60;
  if (!atomSelections.is_none()) {
    unsigned int numSel = nb::len(atomSelections);
    if (numSel == 0) {
      throw ValueErrorException("empty atom selections");
    }
    atomIds.resize(numSel);
    for (unsigned int i = 0; i < numSel; ++i) {
      nb::object inner = atomSelections[i];
      unsigned int numPts = nb::len(inner);
      std::vector<unsigned int> tmpIds(numPts);
      for (unsigned int j = 0; j < numPts; ++j) {
        tmpIds[j] = nb::cast<unsigned int>(inner[j]) - 1;
      }
      atomIds[i] = tmpIds;
    }
    sizeDescriptor = 12 * (numSel + 1);
  }
  std::vector<double> descriptor(sizeDescriptor);
  RDKit::Descriptors::USRCAT(mol, descriptor, atomIds, confId);
  return descriptor;
}

std::vector<double> CalcSlogPVSA(const RDKit::ROMol &mol, nb::object bins,
                                  bool force) {
  std::vector<double> *lbins = nullptr;
  if (!bins.is_none()) {
    unsigned int nBins = nb::len(bins);
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = nb::cast<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcSlogP_VSA(mol, lbins, force);
  delete lbins;
  return res;
}

std::vector<double> CalcSMRVSA(const RDKit::ROMol &mol, nb::object bins,
                                bool force) {
  std::vector<double> *lbins = nullptr;
  if (!bins.is_none()) {
    unsigned int nBins = nb::len(bins);
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = nb::cast<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcSMR_VSA(mol, lbins, force);
  delete lbins;
  return res;
}

std::vector<double> CalcPEOEVSA(const RDKit::ROMol &mol, nb::object bins,
                                  bool force) {
  std::vector<double> *lbins = nullptr;
  if (!bins.is_none()) {
    unsigned int nBins = nb::len(bins);
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = nb::cast<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcPEOE_VSA(mol, lbins, force);
  delete lbins;
  return res;
}

std::vector<double> CalcCustomPropVSA(const RDKit::ROMol &mol,
                                       const std::string customPropName,
                                       nb::list bins, bool force) {
  unsigned int nBins = nb::len(bins);
  std::vector<double> lbins(nBins, 0.0);
  for (unsigned int i = 0; i < nBins; ++i) {
    lbins[i] = nb::cast<double>(bins[i]);
  }
  return RDKit::Descriptors::calcCustomProp_VSA(mol, customPropName, lbins,
                                                force);
}

std::vector<unsigned int> CalcMQNs(const RDKit::ROMol &mol, bool force) {
  return RDKit::Descriptors::calcMQNs(mol, force);
}

unsigned int numSpiroAtoms(const RDKit::ROMol &mol, nb::object pyatoms) {
  std::vector<unsigned int> ats;
  unsigned int res = RDKit::Descriptors::calcNumSpiroAtoms(
      mol, !pyatoms.is_none() ? &ats : nullptr);
  if (!pyatoms.is_none()) {
    nb::list pyres = nb::cast<nb::list>(pyatoms);
    for (const auto d : ats) {
      pyres.append(d);
    }
  }
  return res;
}

unsigned int numBridgeheadAtoms(const RDKit::ROMol &mol, nb::object pyatoms) {
  std::vector<unsigned int> ats;
  unsigned int res = RDKit::Descriptors::calcNumBridgeheadAtoms(
      mol, !pyatoms.is_none() ? &ats : nullptr);
  if (!pyatoms.is_none()) {
    nb::list pyres = nb::cast<nb::list>(pyatoms);
    for (const auto d : ats) {
      pyres.append(d);
    }
  }
  return res;
}

// Python-callable property functor that delegates __call__ to a Python callback
struct PythonPropertyFunctor : public RDKit::Descriptors::PropertyFunctor {
  nb::object callback;

  PythonPropertyFunctor(nb::object cb, const std::string &name,
                        const std::string &version)
      : PropertyFunctor(name, version), callback(cb) {}

  double operator()(const RDKit::ROMol &mol) const override {
    return nb::cast<double>(callback(mol));
  }
};

int registerPropertyHelper(std::shared_ptr<PythonPropertyFunctor> ptr) {
  return RDKit::Descriptors::Properties::registerProperty(
      new PythonPropertyFunctor(*ptr));
}

// Convert nb::object of atom indices to boost::dynamic_bitset
boost::dynamic_bitset<> objectToDynBitset(
    const nb::object &obj,
    boost::dynamic_bitset<>::size_type maxV) {
  boost::dynamic_bitset<> res(maxV);
  if (!obj.is_none()) {
    for (auto item : obj) {
      auto idx = nb::cast<boost::dynamic_bitset<>::size_type>(item);
      if (idx < maxV) {
        res.set(idx);
      }
    }
  }
  return res;
}

double getPartialSurfaceAreaHelper(
    RDKit::Descriptors::DoubleCubicLatticeVolume &self,
    const nb::object &atomIdxs) {
  unsigned int numAtoms = self.mol.getNumAtoms();
  auto atoms = objectToDynBitset(atomIdxs, numAtoms);

  if (atoms.empty()) {
    throw ValueErrorException(
        "No atom indices supplied for Partial Surface Area");
  }

  return self.getPartialSurfaceArea(atoms);
}

double getPartialVolumeHelper(
    RDKit::Descriptors::DoubleCubicLatticeVolume &self,
    const nb::object &atomIdxs) {
  unsigned int numAtoms = self.mol.getNumAtoms();
  auto atoms = objectToDynBitset(atomIdxs, numAtoms);
  if (atoms.empty()) {
    throw ValueErrorException(
        "No atom indices supplied for Partial Surface Area");
  }

  return self.getPartialVolume(atoms);
}

nb::dict getSurfacePointsHelper(
    RDKit::Descriptors::DoubleCubicLatticeVolume &self) {
  const std::map<unsigned int, std::vector<RDGeom::Point3D>> &points =
      self.getSurfacePoints();
  nb::dict surfacePoints;

  for (const auto &it : points) {
    nb::list points3D;
    for (const auto &p : it.second) {
      points3D.append(p);
    }
    surfacePoints[nb::int_(it.first)] = points3D;
  }
  return surfacePoints;
}

}  // namespace

NB_MODULE(rdMolDescriptors, m) {
  m.doc() = "Module containing functions to compute molecular descriptors";

  nb::class_<AtomPairsParameters>(m, "AtomPairsParameters")
      .def_prop_ro_static(
          "version",
          [](nb::object) {
            return std::string(RDKit::AtomPairs::atomPairsVersion);
          })
      .def_prop_ro_static(
          "numTypeBits",
          [](nb::object) { return RDKit::AtomPairs::numTypeBits; })
      .def_prop_ro_static(
          "numPiBits",
          [](nb::object) { return RDKit::AtomPairs::numPiBits; })
      .def_prop_ro_static(
          "numBranchBits",
          [](nb::object) { return RDKit::AtomPairs::numBranchBits; })
      .def_prop_ro_static(
          "numChiralBits",
          [](nb::object) { return RDKit::AtomPairs::numChiralBits; })
      .def_prop_ro_static(
          "codeSize",
          [](nb::object) { return RDKit::AtomPairs::codeSize; })
      .def_prop_ro_static("atomTypes",
                          [](nb::object) { return atomPairTypes; })
      .def_prop_ro_static(
          "numPathBits",
          [](nb::object) { return RDKit::AtomPairs::numPathBits; })
      .def_prop_ro_static(
          "numAtomPairFingerprintBits",
          [](nb::object) {
            return RDKit::AtomPairs::numAtomPairFingerprintBits;
          })
      .def("__setattr__", &safeSetattr);

  m.def("GetAtomPairAtomCode", RDKit::AtomPairs::getAtomCode,
        "atom"_a, "branchSubtract"_a = 0, "includeChirality"_a = false,
        "Returns the atom code (hash) for an atom");

  m.def("GetAtomPairCode", RDKit::AtomPairs::getAtomPairCode,
        "atom1Code"_a, "atom2Code"_a, "distance"_a,
        "includeChirality"_a = false,
        R"DOC(Returns the atom-pair code (hash) for a pair of atoms separated by a
certain number of bonds)DOC");

  m.def("GetAtomPairFingerprint", GetAtomPairFingerprint,
        "mol"_a, "minLength"_a = 1,
        "maxLength"_a = (unsigned int)(RDKit::AtomPairs::maxPathLen - 1),
        "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
        "atomInvariants"_a = nb::none(), "includeChirality"_a = false,
        "use2D"_a = true, "confId"_a = -1,
        "Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect",
        nb::rv_policy::take_ownership);

  m.def("GetHashedAtomPairFingerprint", GetHashedAtomPairFingerprint,
        "mol"_a, "nBits"_a = 2048, "minLength"_a = 1,
        "maxLength"_a = (unsigned int)(RDKit::AtomPairs::maxPathLen - 1),
        "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
        "atomInvariants"_a = nb::none(), "includeChirality"_a = false,
        "use2D"_a = true, "confId"_a = -1,
        "Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect",
        nb::rv_policy::take_ownership);

  m.def("GetHashedAtomPairFingerprintAsBitVect",
        GetHashedAtomPairFingerprintAsBitVect,
        "mol"_a, "nBits"_a = 2048, "minLength"_a = 1,
        "maxLength"_a = (unsigned int)(RDKit::AtomPairs::maxPathLen - 1),
        "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
        "atomInvariants"_a = nb::none(), "nBitsPerEntry"_a = 4,
        "includeChirality"_a = false, "use2D"_a = true, "confId"_a = -1,
        "Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect",
        nb::rv_policy::take_ownership);

  m.def("GetTopologicalTorsionFingerprint", GetTopologicalTorsionFingerprint,
        "mol"_a, "targetSize"_a = 4, "fromAtoms"_a = nb::none(),
        "ignoreAtoms"_a = nb::none(), "atomInvariants"_a = nb::none(),
        "includeChirality"_a = false,
        "Returns the topological-torsion fingerprint for a molecule as a "
        "LongIntSparseIntVect",
        nb::rv_policy::take_ownership);

  m.def("GetHashedTopologicalTorsionFingerprint",
        GetHashedTopologicalTorsionFingerprint,
        "mol"_a, "nBits"_a = 2048, "targetSize"_a = 4,
        "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
        "atomInvariants"_a = nb::none(), "includeChirality"_a = false,
        "Returns the hashed topological-torsion fingerprint for a molecule as "
        "a LongIntSparseIntVect",
        nb::rv_policy::take_ownership);

  m.def("GetHashedTopologicalTorsionFingerprintAsBitVect",
        GetHashedTopologicalTorsionFingerprintAsBitVect,
        "mol"_a, "nBits"_a = 2048, "targetSize"_a = 4,
        "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
        "atomInvariants"_a = nb::none(), "nBitsPerEntry"_a = 4,
        "includeChirality"_a = false,
        "Returns the topological-torsion fingerprint for a molecule as an "
        "ExplicitBitVect",
        nb::rv_policy::take_ownership);

  m.def("GetMorganFingerprint", GetMorganFingerprint,
        "mol"_a, "radius"_a, "invariants"_a = nb::list(),
        "fromAtoms"_a = nb::list(), "useChirality"_a = false,
        "useBondTypes"_a = true, "useFeatures"_a = false, "useCounts"_a = true,
        "bitInfo"_a = nb::none(), "includeRedundantEnvironments"_a = false,
        "Returns a Morgan fingerprint for a molecule",
        nb::rv_policy::take_ownership);

  m.def("GetHashedMorganFingerprint", GetHashedMorganFingerprint,
        "mol"_a, "radius"_a, "nBits"_a = 2048, "invariants"_a = nb::list(),
        "fromAtoms"_a = nb::list(), "useChirality"_a = false,
        "useBondTypes"_a = true, "useFeatures"_a = false,
        "bitInfo"_a = nb::none(), "includeRedundantEnvironments"_a = false,
        "Returns a hashed Morgan fingerprint for a molecule",
        nb::rv_policy::take_ownership);

  m.def("GetMorganFingerprintAsBitVect", GetMorganFingerprintBV,
        "mol"_a, "radius"_a, "nBits"_a = 2048, "invariants"_a = nb::list(),
        "fromAtoms"_a = nb::list(), "useChirality"_a = false,
        "useBondTypes"_a = true, "useFeatures"_a = false,
        "bitInfo"_a = nb::none(), "includeRedundantEnvironments"_a = false,
        "Returns a Morgan fingerprint for a molecule as a bit vector",
        nb::rv_policy::take_ownership);

  m.attr("_MorganFingerprint_version") =
      RDKit::MorganFingerprints::morganFingerprintVersion;

  m.def("GetConnectivityInvariants", GetConnectivityInvariants,
        "mol"_a, "includeRingMembership"_a = true,
        "Returns connectivity invariants (ECFP-like) for a molecule.");
  m.attr("_ConnectivityInvariants_version") =
      RDKit::MorganFingerprints::morganConnectivityInvariantVersion;

  m.def("GetFeatureInvariants", GetFeatureInvariants,
        "mol"_a, "Returns feature invariants (FCFP-like) for a molecule.");
  m.attr("_FeatureInvariants_version") =
      RDKit::MorganFingerprints::morganFeatureInvariantVersion;

  // USR descriptor
  m.def("GetUSR", GetUSR, "mol"_a, "confId"_a = -1,
        "Returns a USR descriptor for one conformer of a molecule");

  m.def("GetUSRDistributions", GetUSRDistributions,
        "coords"_a, "points"_a = nb::none(),
        "Returns the four USR distance distributions for a set of coordinates");

  m.def("GetUSRDistributionsFromPoints", GetUSRDistributionsFromPoints,
        "coords"_a, "points"_a,
        "Returns the USR distance distributions for a set of coordinates and points");

  m.def("GetUSRFromDistributions", GetUSRFromDistributions,
        "distances"_a,
        "Returns the USR descriptor from a set of distance distributions");

  m.def("GetUSRScore", GetUSRScore,
        "descriptor1"_a, "descriptor2"_a, "weights"_a = nb::list(),
        "Returns the USR score for two USR or USRCAT descriptors");

  m.def("GetUSRCAT", GetUSRCAT,
        "mol"_a, "atomSelections"_a = nb::none(), "confId"_a = -1,
        "Returns a USRCAT descriptor for one conformer of a molecule");

  m.def("_CalcCrippenContribs", computeCrippenContribs,
        "mol"_a, "force"_a = false, "atomTypes"_a = nb::none(),
        "atomTypeLabels"_a = nb::none(),
        R"DOC(returns (as a list of 2-tuples) the contributions of each atom to
the Wildman-Cripppen logp and mr value)DOC");

  m.def("CalcCrippenDescriptors", calcCrippenDescriptors,
        "mol"_a, "includeHs"_a = true, "force"_a = false,
        "returns a 2-tuple with the Wildman-Crippen logp,mr values");
  m.attr("_CalcCrippenDescriptors_version") =
      RDKit::Descriptors::crippenVersion;

  m.def("CalcLabuteASA", RDKit::Descriptors::calcLabuteASA,
        "mol"_a, "includeHs"_a = true, "force"_a = false,
        "returns the Labute ASA value for a molecule");
  m.attr("_CalcLabuteASA_version") = RDKit::Descriptors::labuteASAVersion;

  m.def("_CalcLabuteASAContribs", computeASAContribs,
        "mol"_a, "includeHs"_a = true, "force"_a = false,
        "returns a list of atomic contributions to the Labute ASA");

  m.def("CalcTPSA", RDKit::Descriptors::calcTPSA,
        "mol"_a, "force"_a = false, "includeSandP"_a = false,
        "returns the TPSA value for a molecule");
  m.attr("_CalcTPSA_version") = RDKit::Descriptors::tpsaVersion;

  m.def("_CalcTPSAContribs", computeTPSAContribs,
        "mol"_a, "force"_a = false, "includeSandP"_a = false,
        "returns a list of atomic contributions to the TPSA");

  m.def("_CalcMolWt", RDKit::Descriptors::calcAMW,
        "mol"_a, "onlyHeavy"_a = false,
        "returns the molecule's molecular weight");
  m.attr("_CalcMolWt_version") = "1.0.0";

  m.def("CalcExactMolWt", RDKit::Descriptors::calcExactMW,
        "mol"_a, "onlyHeavy"_a = false,
        "returns the molecule's exact molecular weight");
  m.attr("_CalcExactMolWt_version") = "1.0.0";

  m.def("CalcMolFormula", RDKit::Descriptors::calcMolFormula,
        "mol"_a, "separateIsotopes"_a = false, "abbreviateHIsotopes"_a = true,
        "returns the molecule's formula");
  m.attr("_CalcMolFormula_version") = "1.3.0";

  m.def("CalcNumLipinskiHBD", RDKit::Descriptors::calcLipinskiHBD,
        "mol"_a, "returns the number of Lipinski H-bond donors for a molecule");
  m.attr("_CalcNumLipinskiHBD_version") =
      RDKit::Descriptors::lipinskiHBDVersion;

  m.def("CalcNumLipinskiHBA", RDKit::Descriptors::calcLipinskiHBA,
        "mol"_a,
        "returns the number of Lipinski H-bond acceptors for a molecule");
  m.attr("_CalcNumLipinskiHBA_version") =
      RDKit::Descriptors::lipinskiHBAVersion;

  m.def("CalcNumHBD", RDKit::Descriptors::calcNumHBD,
        "mol"_a, "returns the number of H-bond donors for a molecule");
  m.attr("_CalcNumHBD_version") = RDKit::Descriptors::NumHBDVersion;

  m.def("CalcNumHBA", RDKit::Descriptors::calcNumHBA,
        "mol"_a, "returns the number of H-bond acceptors for a molecule");
  m.attr("_CalcNumHBA_version") = RDKit::Descriptors::NumHBAVersion;

  // exposes NumRotatableBondsOptions enum
  nb::enum_<RDKit::Descriptors::NumRotatableBondsOptions>(
      m, "NumRotatableBondsOptions",
      R"DOC(Options for generating rotatable bonds
NonStrict - standard loose definitions
Strict - stricter definition excluding amides, esters, etc
StrictLinkages - adds rotors between rotatable bonds
Default - Current RDKit default)DOC")
      .value("NonStrict", RDKit::Descriptors::NonStrict)
      .value("Strict", RDKit::Descriptors::Strict)
      .value("StrictLinkages", RDKit::Descriptors::StrictLinkages)
      .value("Default", RDKit::Descriptors::Default);

#ifdef RDK_USE_STRICT_ROTOR_DEFINITION
  std::string calcNumRotatableBondsDoc =
      R"DOC(returns the number of rotatable bonds for a molecule.
strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
strict = NumRotatableBondsOptions.Strict - (default) does not count things like
         amide or ester bonds
strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
   systems.
   - Single bonds between aliphatic ring Cs are always rotatable. This
     means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now
     considered rotatable; it was not before
   - Heteroatoms in the linked rings no longer affect whether or not
     the linking bond is rotatable
   - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
      considered non-rotatable)DOC";
#else
  std::string calcNumRotatableBondsDoc =
      R"DOC(returns the number of rotatable bonds for a molecule.
strict = NumRotatableBondsOptions.NonStrict - (default) Simple rotatable bond definition.
strict = NumRotatableBondsOptions.Strict - does not count things like
         amide or ester bonds
strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
   systems.
   - Single bonds between aliphatic ring Cs are always rotatable. This
     means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now
     considered rotatable; it was not before
   - Heteroatoms in the linked rings no longer affect whether or not
     the linking bond is rotatable
   - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
      considered non-rotatable)DOC";
#endif

  m.def("CalcNumRotatableBonds",
        (unsigned int (*)(const RDKit::ROMol &,
                          bool))RDKit::Descriptors::calcNumRotatableBonds,
        "mol"_a, "strict"_a, calcNumRotatableBondsDoc.c_str());

  m.def("CalcNumRotatableBonds",
        (unsigned int (*)(const RDKit::ROMol &,
                          RDKit::Descriptors::NumRotatableBondsOptions))
            RDKit::Descriptors::calcNumRotatableBonds,
        "mol"_a, "strict"_a = RDKit::Descriptors::Default,
        calcNumRotatableBondsDoc.c_str());
  m.attr("_CalcNumRotatableBonds_version") =
      RDKit::Descriptors::NumRotatableBondsVersion;

  m.def("CalcNumRings", RDKit::Descriptors::calcNumRings,
        "mol"_a, "returns the number of rings for a molecule");
  m.attr("_CalcNumRings_version") = RDKit::Descriptors::NumRingsVersion;

  m.def("CalcNumAromaticRings", RDKit::Descriptors::calcNumAromaticRings,
        "mol"_a, "returns the number of aromatic rings for a molecule");
  m.attr("_CalcNumAromaticRings_version") =
      RDKit::Descriptors::NumAromaticRingsVersion;

  m.def("CalcNumSaturatedRings", RDKit::Descriptors::calcNumSaturatedRings,
        "mol"_a, "returns the number of saturated rings for a molecule");
  m.attr("_CalcNumSaturatedRings_version") =
      RDKit::Descriptors::NumSaturatedRingsVersion;

  m.def("CalcNumHeterocycles", RDKit::Descriptors::calcNumHeterocycles,
        "mol"_a, "returns the number of heterocycles for a molecule");
  m.attr("_CalcNumHeterocycles_version") =
      RDKit::Descriptors::NumHeterocyclesVersion;

  m.def("CalcNumAromaticHeterocycles",
        RDKit::Descriptors::calcNumAromaticHeterocycles,
        "mol"_a, "returns the number of aromatic heterocycles for a molecule");
  m.attr("_CalcNumAromaticHeterocycles_version") =
      RDKit::Descriptors::NumAromaticHeterocyclesVersion;

  m.def("CalcNumAromaticCarbocycles",
        RDKit::Descriptors::calcNumAromaticCarbocycles,
        "mol"_a, "returns the number of aromatic carbocycles for a molecule");
  m.attr("_CalcNumAromaticCarbocycles_version") =
      RDKit::Descriptors::NumAromaticCarbocyclesVersion;

  m.def("CalcNumSaturatedHeterocycles",
        RDKit::Descriptors::calcNumSaturatedHeterocycles,
        "mol"_a, "returns the number of saturated heterocycles for a molecule");
  m.attr("_CalcNumSaturatedHeterocycles_version") =
      RDKit::Descriptors::NumSaturatedHeterocyclesVersion;

  m.def("CalcNumSaturatedCarbocycles",
        RDKit::Descriptors::calcNumSaturatedCarbocycles,
        "mol"_a, "returns the number of saturated carbocycles for a molecule");
  m.attr("_CalcNumSaturatedCarbocycles_version") =
      RDKit::Descriptors::NumSaturatedCarbocyclesVersion;

  m.def("CalcNumAliphaticRings", RDKit::Descriptors::calcNumAliphaticRings,
        "mol"_a,
        "returns the number of aliphatic (containing at least one non-aromatic "
        "bond) rings for a molecule");
  m.attr("_CalcNumAliphaticRings_version") =
      RDKit::Descriptors::NumAliphaticRingsVersion;

  m.def("CalcNumAliphaticHeterocycles",
        RDKit::Descriptors::calcNumAliphaticHeterocycles,
        "mol"_a,
        "returns the number of aliphatic (containing at least one non-aromatic "
        "bond) heterocycles for a molecule");
  m.attr("_CalcNumAliphaticHeterocycles_version") =
      RDKit::Descriptors::NumAliphaticHeterocyclesVersion;

  m.def("CalcNumAliphaticCarbocycles",
        RDKit::Descriptors::calcNumAliphaticCarbocycles,
        "mol"_a,
        "returns the number of aliphatic (containing at least one non-aromatic "
        "bond) carbocycles for a molecule");
  m.attr("_CalcNumAliphaticCarbocycles_version") =
      RDKit::Descriptors::NumAliphaticCarbocyclesVersion;

  m.def("CalcNumHeavyAtoms", RDKit::Descriptors::calcNumHeavyAtoms,
        "mol"_a, "returns the number of heavy atoms for a molecule");
  m.attr("_CalcNumHeavyAtoms_version") =
      RDKit::Descriptors::NumHeavyAtomsVersion;

  m.def("CalcNumAtoms", RDKit::Descriptors::calcNumAtoms,
        "mol"_a, "returns the total number of atoms for a molecule");
  m.attr("_CalcNumAtoms_version") = RDKit::Descriptors::NumAtomsVersion;

  m.def("CalcNumHeteroatoms", RDKit::Descriptors::calcNumHeteroatoms,
        "mol"_a, "returns the number of heteroatoms for a molecule");
  m.attr("_CalcNumHeteroatoms_version") =
      RDKit::Descriptors::NumHeteroatomsVersion;

  m.def("CalcNumAmideBonds", RDKit::Descriptors::calcNumAmideBonds,
        "mol"_a, "returns the number of amide bonds in a molecule");
  m.attr("_CalcNumAmideBonds_version") =
      RDKit::Descriptors::NumAmideBondsVersion;

  m.def("CalcFractionCSP3", RDKit::Descriptors::calcFractionCSP3,
        "mol"_a, "returns the fraction of C atoms that are SP3 hybridized");
  m.attr("_CalcFractionCSP3_version") =
      RDKit::Descriptors::FractionCSP3Version;

  m.def("SlogP_VSA_", CalcSlogPVSA,
        "mol"_a, "bins"_a = nb::none(), "force"_a = false,
        "returns the SlogP VSA contributions for a molecule");

  m.def("SMR_VSA_", CalcSMRVSA,
        "mol"_a, "bins"_a = nb::none(), "force"_a = false,
        "returns the SMR VSA contributions for a molecule");

  m.def("PEOE_VSA_", CalcPEOEVSA,
        "mol"_a, "bins"_a = nb::none(), "force"_a = false,
        "returns the PEOE VSA contributions for a molecule");

  m.def("CustomProp_VSA_", CalcCustomPropVSA,
        "mol"_a, "customPropName"_a, "bins"_a, "force"_a = false,
        "returns the VSA contributions based on a custom property for a molecule");

  m.def("MQNs_", CalcMQNs,
        "mol"_a, "force"_a = false,
        "returns the MQN descriptors for a molecule");

  m.def("CalcChiNv", RDKit::Descriptors::calcChiNv,
        "mol"_a, "n"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChiNv_version") = RDKit::Descriptors::chiNvVersion;

  m.def("CalcChi0v", RDKit::Descriptors::calcChi0v,
        "mol"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChi0v_version") = RDKit::Descriptors::chi0vVersion;

  m.def("CalcChi1v", RDKit::Descriptors::calcChi1v,
        "mol"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChi1v_version") = RDKit::Descriptors::chi1vVersion;

  m.def("CalcChi2v", RDKit::Descriptors::calcChi2v,
        "mol"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChi2v_version") = RDKit::Descriptors::chi2vVersion;

  m.def("CalcChi3v", RDKit::Descriptors::calcChi3v,
        "mol"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChi3v_version") = RDKit::Descriptors::chi3vVersion;

  m.def("CalcChi4v", RDKit::Descriptors::calcChi4v,
        "mol"_a, "force"_a = false,
        "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
        "(1991)");
  m.attr("_CalcChi4v_version") = RDKit::Descriptors::chi4vVersion;

  m.def("CalcChiNn", RDKit::Descriptors::calcChiNn,
        "mol"_a, "n"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChiNn_version") = RDKit::Descriptors::chiNnVersion;

  m.def("CalcChi0n", RDKit::Descriptors::calcChi0n,
        "mol"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChi0n_version") = RDKit::Descriptors::chi0nVersion;

  m.def("CalcChi1n", RDKit::Descriptors::calcChi1n,
        "mol"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChi1n_version") = RDKit::Descriptors::chi1nVersion;

  m.def("CalcChi2n", RDKit::Descriptors::calcChi2n,
        "mol"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChi2n_version") = RDKit::Descriptors::chi2nVersion;

  m.def("CalcChi3n", RDKit::Descriptors::calcChi3n,
        "mol"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChi3n_version") = RDKit::Descriptors::chi3nVersion;

  m.def("CalcChi4n", RDKit::Descriptors::calcChi4n,
        "mol"_a, "force"_a = false,
        "Similar to ChiXv, but uses uses nVal instead of valence. This makes "
        "a big difference after we get out of the first row.");
  m.attr("_CalcChi4n_version") = RDKit::Descriptors::chi4nVersion;

  m.def("CalcHallKierAlpha", hkAlphaHelper,
        "mol"_a, "atomContribs"_a = nb::none(),
        R"DOC(From equation (58) of Rev. Comp. Chem. vol 2, 367-422, (1991).
NOTE: Because hybridization is used to calculate this, results may
differ from other implementations which have different conventions for
assigning hybridization)DOC");
  m.attr("_CalcHallKierAlpha_version") =
      RDKit::Descriptors::hallKierAlphaVersion;

  m.def("CalcKappa1", RDKit::Descriptors::calcKappa1,
        "mol"_a,
        R"DOC(From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991)
NOTE: Because hybridization is used to calculate this, results may
differ from other implementations which have different conventions for
assigning hybridization)DOC");
  m.attr("_CalcKappa1_version") = RDKit::Descriptors::kappa1Version;

  m.def("CalcKappa2", RDKit::Descriptors::calcKappa2,
        "mol"_a,
        R"DOC(From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991)
NOTE: Because hybridization is used to calculate this, results may
differ from other implementations which have different conventions for
assigning hybridization)DOC");
  m.attr("_CalcKappa2_version") = RDKit::Descriptors::kappa2Version;

  m.def("CalcKappa3", RDKit::Descriptors::calcKappa3,
        "mol"_a,
        R"DOC(From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, (1991)
NOTE: Because hybridization is used to calculate this, results may
differ from other implementations which have different conventions for
assigning hybridization)DOC");
  m.attr("_CalcKappa3_version") = RDKit::Descriptors::kappa3Version;

  m.def("CalcPhi", RDKit::Descriptors::calcPhi,
        "mol"_a,
        "From Quantitative Structure-Activity Relationships 8, 221-224 (1989). "
        "NOTE: Because hybridization is used to calculate this, results may "
        "differ from other implementations which have different conventions for "
        "assigning hybridization");
  m.attr("_CalcPhi_version") = RDKit::Descriptors::PhiVersion;

  m.def("GetMACCSKeysFingerprint",
        RDKit::MACCSFingerprints::getFingerprintAsBitVect,
        "mol"_a,
        "Returns the MACCS keys for a molecule as an ExplicitBitVect",
        nb::rv_policy::take_ownership);

  m.attr("_GetAtomFeatures_version") = RDKit::Descriptors::AtomFeatVersion;
  m.def("GetAtomFeatures", GetAtomFeatures,
        "mol"_a, "atomid"_a, "addchiral"_a = false,
        "Returns the Atom Features vector");

  m.attr("_CalcNumSpiroAtoms_version") =
      RDKit::Descriptors::NumSpiroAtomsVersion;
  m.def("CalcNumSpiroAtoms", numSpiroAtoms,
        "mol"_a, "atoms"_a = nb::none(),
        "Returns the number of spiro atoms (atoms shared between rings that "
        "share exactly one atom)");

  m.attr("_CalcNumBridgeheadAtoms_version") =
      RDKit::Descriptors::NumBridgeheadAtomsVersion;
  m.def("CalcNumBridgeheadAtoms", numBridgeheadAtoms,
        "mol"_a, "atoms"_a = nb::none(),
        "Returns the number of bridgehead atoms (atoms shared between rings "
        "that share at least two bonds)");

  m.attr("_CalcNumAtomStereoCenters_version") =
      RDKit::Descriptors::NumAtomStereoCentersVersion;
  m.def("CalcNumAtomStereoCenters",
        RDKit::Descriptors::numAtomStereoCenters,
        "mol"_a,
        "Returns the total number of atomic stereocenters (specified and "
        "unspecified)");

  m.attr("_CalcNumUnspecifiedAtomStereoCenters_version") =
      RDKit::Descriptors::NumUnspecifiedAtomStereoCentersVersion;
  m.def("CalcNumUnspecifiedAtomStereoCenters",
        RDKit::Descriptors::numUnspecifiedAtomStereoCenters,
        "mol"_a, "Returns the number of unspecified atomic stereocenters");

  nb::class_<RDKit::Descriptors::PropertyFunctor>(
      m, "PropertyFunctor",
      R"DOC(Property computation class stored in the property registry.
See rdkit.Chem.rdMolDescriptor.Properties.GetProperty and
rdkit.Chem.Descriptor.Properties.PropertyFunctor for creating new ones)DOC")
      .def("__call__", &RDKit::Descriptors::PropertyFunctor::operator(),
           "mol"_a, "Compute the property for the specified molecule")
      .def("GetName", &RDKit::Descriptors::PropertyFunctor::getName,
           "Return the name of the property to calculate")
      .def("GetVersion", &RDKit::Descriptors::PropertyFunctor::getVersion,
           "Return the version of the calculated property");

  nb::class_<RDKit::Descriptors::Properties>(
      m, "Properties",
      R"DOC(Property computation and registry system.  To compute all registered properties:
mol = Chem.MolFromSmiles('c1ccccc1')
properties = rdMolDescriptors.Properties()
for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
  print(name, value)

To compute a subset
properties = rdMolDescriptors.Properties(['exactmw', 'lipinskiHBA'])
for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
  print(name, value))DOC")
      .def(nb::init<>())
      .def(nb::init<const std::vector<std::string> &>(), "propNames"_a)
      .def("GetPropertyNames",
           &RDKit::Descriptors::Properties::getPropertyNames,
           "Return the property names computed by this instance")
      .def("ComputeProperties",
           &RDKit::Descriptors::Properties::computeProperties,
           "mol"_a, "annotateMol"_a = false,
           "Return a list of computed properties, if annotateMol==True, "
           "annotate the molecule with the computed properties.")
      .def("AnnotateProperties",
           &RDKit::Descriptors::Properties::annotateProperties,
           "mol"_a,
           "Annotate the molecule with the computed properties.  These "
           "properties will be available as SDData or from mol.GetProp(prop)")
      .def_static("GetAvailableProperties",
                  &RDKit::Descriptors::Properties::getAvailableProperties,
                  "Return all available property names that can be computed")
      .def_static("GetProperty",
                  [](const std::string &name)
                      -> RDKit::Descriptors::PropertyFunctor * {
                    return RDKit::Descriptors::Properties::getProperty(name)
                        .get();
                  },
                  "propName"_a, nb::rv_policy::reference,
                  "Return the named property if it exists")
      .def_static("RegisterProperty", &registerPropertyHelper,
                  "propertyFunctor"_a,
                  "Register a new property object (not thread safe)");

  nb::class_<PythonPropertyFunctor,
             RDKit::Descriptors::PropertyFunctor>(
      m, "PythonPropertyFunctor", "")
      .def(nb::init<nb::object, const std::string &, const std::string &>(),
           "callback"_a, "name"_a, "version"_a)
      .def("__call__", &PythonPropertyFunctor::operator(),
           "mol"_a, "Compute the property for the specified molecule");

  nb::class_<Queries::RangeQuery<double, RDKit::ROMol const &, true>>(
      m, "PropertyRangeQuery",
      "Property Range Query for a molecule.  Match(mol) -> true if in range")
      .def("Match",
           &Queries::RangeQuery<double, RDKit::ROMol const &, true>::Match,
           "what"_a);

  m.def("MakePropertyRangeQuery",
        RDKit::Descriptors::makePropertyRangeQuery,
        "name"_a, "min"_a, "max"_a,
        R"DOC(Generates a Range property for the specified property, between min and max
query = MakePropertyRangeQuery('exactmw', 0, 500)
query.Match( mol ))DOC",
        nb::rv_policy::take_ownership);

  nb::class_<RDKit::Descriptors::DoubleCubicLatticeVolume>(
      m, "DoubleCubicLatticeVolume",
      "Class for the Double Cubic Lattice Volume method")
      .def(nb::init<const RDKit::ROMol &, bool, bool, double, int>(),
           "mol"_a, "isProtein"_a = false, "includeLigand"_a = true,
           "probeRadius"_a = 1.4, "confId"_a = -1)
      .def("__init__",
           [](RDKit::Descriptors::DoubleCubicLatticeVolume *self,
              const RDKit::ROMol &mol, const nb::list &radii, bool isProtein,
              bool includeLigand, double probeRadius, int confId) {
             std::vector<double> radiiAsVector;
             radiiAsVector.reserve(mol.getNumAtoms());
             pythonObjectToVect<double>(nb::cast<nb::object>(radii),
                                       radiiAsVector);
             new (self) RDKit::Descriptors::DoubleCubicLatticeVolume(
                 mol, std::move(radiiAsVector), isProtein, includeLigand,
                 probeRadius, confId);
           },
           "mol"_a, "radii"_a, "isProtein"_a = false,
           "includeLigand"_a = true, "probeRadius"_a = 1.4, "confId"_a = -1,
           R"DOC(ARGUMENTS:
   - mol: molecule or protein under consideration
   - radii: radii for atoms of input mol (get using GetPeriodicTable or provide custom list)
   - isProtein: flag to indicate if the input is a protein (default=False, free ligand).
   - includeLigand: flag to include or exclude a bound ligand when input is a protein (default=True)
   - probeRadius: radius of the solvent probe (default=1.2)
   - confId: conformer ID to consider (default=-1))DOC")
      .def("GetSurfaceArea",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getSurfaceArea,
           "Get the Surface Area of the Molecule or Protein")
      .def("GetAtomSurfaceArea",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getAtomSurfaceArea,
           "atom_idx"_a, "Get the surface area of atom with atom_idx")
      .def("GetPolarSurfaceArea",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getPolarSurfaceArea,
           "includeSandP"_a = false, "includeHs"_a = false,
           "Get the Polar Surface Area of the Molecule or Protein")
      .def("GetPartialSurfaceArea", &getPartialSurfaceAreaHelper,
           "atomIndices"_a,
           "Get the Partial Surface Area of the Molecule or Protein for "
           "specified subset of atoms")
      .def("GetSurfacePoints", &getSurfacePointsHelper,
           "Get the set of points representing the surface")
      .def("GetVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getVolume,
           "Get the Total Volume of the Molecule or Protein")
      .def("GetVDWVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getVDWVolume,
           "Get the van der Waals Volume of the Molecule or Protein")
      .def("GetAtomVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getAtomVolume,
           "atomIdx"_a, "solventRadius"_a,
           "Get the volume atom of atom_idx with volume for specified Probe Radius")
      .def("GetPolarVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getPolarVolume,
           "includeSandP"_a = false, "includeHs"_a = false,
           "Get the Polar Volume of the Molecule or Protein")
      .def("GetPartialVolume", &getPartialVolumeHelper,
           "atomIdx"_a,
           "Get the Partial Volume of the Molecule or Protein for specified "
           "subset of atoms")
      .def("GetCompactness",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getCompactness,
           "Get the Compactness of the Protein")
      .def("GetPackingDensity",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getPackingDensity,
           "Get the PackingDensity of the Protein");

#ifdef RDK_BUILD_DESCRIPTORS3D
  m.attr("_CalcCoulombMat_version") = RDKit::Descriptors::CoulombMatVersion;
  m.def("CalcCoulombMat", calcCoulombMat,
        "mol"_a, "confId"_a = -1,
        "Returns severals Coulomb randomized matrices");

  m.attr("_CalcEMMcharges_version") = RDKit::Descriptors::EEMVersion;
  m.def("CalcEEMcharges", calcEEMcharges,
        "mol"_a, "confId"_a = -1, "Returns EEM atomic partial charges");

  m.attr("_CalcWHIM_version") = RDKit::Descriptors::WHIMVersion;
  m.def("CalcWHIM", calcWHIMs,
        "mol"_a, "confId"_a = -1, "thresh"_a = 0.001,
        "CustomAtomProperty"_a = "",
        "Returns the WHIM descriptors vector");

  m.attr("_CalcGETAWAY_version") = RDKit::Descriptors::GETAWAYVersion;
  m.def("CalcGETAWAY", calcGETAWAYs,
        "mol"_a, "confId"_a = -1, "precision"_a = 2,
        "CustomAtomProperty"_a = "",
        "Returns the GETAWAY descriptors vector");

  m.attr("_CalcRDF_version") = RDKit::Descriptors::RDFVersion;
  m.def("CalcRDF", calcRDFs,
        "mol"_a, "confId"_a = -1, "CustomAtomProperty"_a = "",
        "Returns radial distribution fonction descriptors (RDF)");

  m.attr("_CalcMORSE_version") = RDKit::Descriptors::MORSEVersion;
  m.def("CalcMORSE", calcMORSEs,
        "mol"_a, "confId"_a = -1, "CustomAtomProperty"_a = "",
        "Returns Molecule Representation of Structures based on Electron "
        "diffraction descriptors");

  m.attr("_CalcAUTOCORR3D_version") = RDKit::Descriptors::AUTOCORR3DVersion;
  m.def("CalcAUTOCORR3D", calcAUTOCORR3Ds,
        "mol"_a, "confId"_a = -1, "CustomAtomProperty"_a = "",
        "Returns 3D Autocorrelation descriptors vector");

  m.attr("_CalcPBF_version") = RDKit::Descriptors::PBFVersion;
  m.def("CalcPBF", RDKit::Descriptors::PBF,
        "mol"_a, "confId"_a = -1,
        "Returns the PBF (plane of best fit) descriptor "
        "(https://doi.org/10.1021/ci300293f)");

  m.attr("_CalcNPR1_version") = RDKit::Descriptors::NPR1Version;
  m.def("CalcNPR1", RDKit::Descriptors::NPR1,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcNPR2_version") = RDKit::Descriptors::NPR2Version;
  m.def("CalcNPR2", RDKit::Descriptors::NPR2,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcPMI1_version") = RDKit::Descriptors::PMI1Version;
  m.def("CalcPMI1", RDKit::Descriptors::PMI1,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcPMI2_version") = RDKit::Descriptors::PMI2Version;
  m.def("CalcPMI2", RDKit::Descriptors::PMI2,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcPMI3_version") = RDKit::Descriptors::PMI3Version;
  m.def("CalcPMI3", RDKit::Descriptors::PMI3,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcRadiusOfGyration_version") =
      RDKit::Descriptors::radiusOfGyrationVersion;
  m.def("CalcRadiusOfGyration", RDKit::Descriptors::radiusOfGyration,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcInertialShapeFactor_version") =
      RDKit::Descriptors::inertialShapeFactorVersion;
  m.def("CalcInertialShapeFactor", RDKit::Descriptors::inertialShapeFactor,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcEccentricity_version") = RDKit::Descriptors::eccentricityVersion;
  m.def("CalcEccentricity", RDKit::Descriptors::eccentricity,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcAsphericity_version") = RDKit::Descriptors::asphericityVersion;
  m.def("CalcAsphericity", RDKit::Descriptors::asphericity,
        "mol"_a, "confId"_a = -1, "useAtomicMasses"_a = true,
        "force"_a = true, "");

  m.attr("_CalcSpherocityIndex_version") =
      RDKit::Descriptors::spherocityIndexVersion;
  m.def("CalcSpherocityIndex", RDKit::Descriptors::spherocityIndex,
        "mol"_a, "confId"_a = -1, "force"_a = true, "");

  m.attr("_CalcAUTOCORR2D_version") = RDKit::Descriptors::AUTOCORR2DVersion;
  m.def("CalcAUTOCORR2D", calcAUTOCORR2Ds,
        "mol"_a, "CustomAtomProperty"_a = "",
        "Returns 2D Autocorrelation descriptors vector");
#endif

#ifdef RDK_HAS_EIGEN3
  m.attr("_BCUT2D_version") = RDKit::Descriptors::BCUT2DVersion;

  std::pair<double, double> (*BCUT_atomprops)(
      const RDKit::ROMol &, const std::string &) = &RDKit::Descriptors::BCUT2D;

  m.def("BCUT2D", BCUT, "mol"_a,
        R"DOC(Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
Diagonal elements are (currently) atomic mass, gasteiger charge,
crippen logP and crippen MRReturns the 2D BCUT2D descriptors vector as described in
returns [mass eigen value high, mass eigen value low,
         gasteiger charge eigenvalue high, gasteiger charge low,
         crippen lowgp  eigenvalue high, crippen lowgp  low,
         crippen mr eigenvalue high, crippen mr low])DOC");

  m.def("BCUT2D", BCUT2D_list, "mol"_a, "atom_props"_a,
        R"DOC(Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule
and the specified atom props
 atom_props must be a list or tuple of floats equal in
size to the number of atoms in mol)DOC");

  m.def("BCUT2D", BCUT2D_tuple, "mol"_a, "atom_props"_a,
        R"DOC(Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule
and the specified atom props
 atom_props must be a list or tuple of floats equal in
size to the number of atoms in mol)DOC");

  m.def("BCUT2D", BCUT_atomprops, "mol"_a, "atom_propname"_a,
        R"DOC(Returns a 2D BCUT (eigen value high, eigen value low) given the
molecule and the specified atom prop name
atom_propname must exist on each atom and be convertible to a float)DOC");

  m.def("CalcOxidationNumbers", RDKit::Descriptors::calcOxidationNumbers,
        "mol"_a,
        "Adds the oxidation number/state to the atoms of a molecule as"
        " property OxidationNumber on each atom.  Use Pauling"
        " electronegativities."
        "  This is experimental code, still under development.");
#endif
}
