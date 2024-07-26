//
//  Copyright (C) 2007-2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <RDBoost/Wrap.h>
#include <GraphMol/Atom.h>
#include <GraphMol/GraphMol.h>

#include <RDGeneral/BoostStartInclude.h>
#include <RDGeneral/BoostEndInclude.h>
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

#include <vector>

namespace python = boost::python;

namespace {
std::vector<unsigned int> atomPairTypes(
    RDKit::AtomPairs::atomNumberTypes,
    RDKit::AtomPairs::atomNumberTypes +
        sizeof(RDKit::AtomPairs::atomNumberTypes) / sizeof(unsigned int));
python::tuple computeASAContribs(const RDKit::ROMol &mol, bool includeHs = true,
                                 bool force = false) {
  std::vector<double> contribs(mol.getNumAtoms());
  double hContrib = 0.0;
  RDKit::Descriptors::getLabuteAtomContribs(mol, contribs, hContrib, includeHs,
                                            force);
  python::tuple pycontribs(contribs);
  return python::make_tuple(contribs, hContrib);
}
python::tuple computeTPSAContribs(const RDKit::ROMol &mol, bool force,
                                  bool includeSandP) {
  std::vector<double> contribs(mol.getNumAtoms());
  RDKit::Descriptors::getTPSAAtomContribs(mol, contribs, force, includeSandP);
  python::tuple pycontribs(contribs);
  return pycontribs;
}

python::list computeCrippenContribs(
    const RDKit::ROMol &mol, bool force = false,
    python::list atomTypes = python::list(),
    python::list atomTypeLabels = python::list()) {
  std::vector<unsigned int> *tAtomTypes = nullptr;
  std::vector<std::string> *tAtomTypeLabels = nullptr;
  if (python::extract<unsigned int>(atomTypes.attr("__len__")()) != 0) {
    if (python::extract<unsigned int>(atomTypes.attr("__len__")()) !=
        mol.getNumAtoms()) {
      throw_value_error(
          "if atomTypes vector is provided, it must be as long as the number "
          "of atoms");
    } else {
      tAtomTypes = new std::vector<unsigned int>(mol.getNumAtoms(), 0);
    }
  }
  if (python::extract<unsigned int>(atomTypeLabels.attr("__len__")()) != 0) {
    if (python::extract<unsigned int>(atomTypeLabels.attr("__len__")()) !=
        mol.getNumAtoms()) {
      throw_value_error(
          "if atomTypeLabels vector is provided, it must be as long as the "
          "number of atoms");
    } else {
      tAtomTypeLabels = new std::vector<std::string>(mol.getNumAtoms(), "");
    }
  }

  std::vector<double> logpContribs(mol.getNumAtoms());
  std::vector<double> mrContribs(mol.getNumAtoms());

  RDKit::Descriptors::getCrippenAtomContribs(
      mol, logpContribs, mrContribs, force, tAtomTypes, tAtomTypeLabels);
  python::list pycontribs;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    pycontribs.append(python::make_tuple(logpContribs[i], mrContribs[i]));
  }
  if (tAtomTypes) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      atomTypes[i] = (*tAtomTypes)[i];
    }
    delete tAtomTypes;
  }
  if (tAtomTypeLabels) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      atomTypeLabels[i] = (*tAtomTypeLabels)[i];
    }
    delete tAtomTypeLabels;
  }
  return pycontribs;
}
python::tuple calcCrippenDescriptors(const RDKit::ROMol &mol,
                                     bool includeHs = true,
                                     bool force = false) {
  double logp, mr;
  RDKit::Descriptors::calcCrippenDescriptors(mol, logp, mr, includeHs, force);
  return python::make_tuple(logp, mr);
}

#ifdef RDK_BUILD_DESCRIPTORS3D

python::tuple calcCoulombMat(const RDKit::ROMol &mol, int confId) {
  std::vector<std::vector<double>> results;
  RDKit::Descriptors::CoulombMat(mol, results, confId);
  python::list result;
  for (auto &res : results) {
    result.append(res);
  }
  return python::tuple(result);
}

python::list calcEEMcharges(RDKit::ROMol &mol, int confId) {
  std::vector<double> res;
  RDKit::Descriptors::EEM(mol, res, confId);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcWHIMs(const RDKit::ROMol &mol, int confId, double thresh,
                       const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::WHIM(mol, res, confId, thresh, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcGETAWAYs(const RDKit::ROMol &mol, int confId, double precision,
                          const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::GETAWAY(mol, res, confId, precision, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcRDFs(const RDKit::ROMol &mol, int confId,
                      const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::RDF(mol, res, confId, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcMORSEs(const RDKit::ROMol &mol, int confId,
                        const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::MORSE(mol, res, confId, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcAUTOCORR3Ds(const RDKit::ROMol &mol, int confId,
                             const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::AUTOCORR3D(mol, res, confId, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list calcAUTOCORR2Ds(const RDKit::ROMol &mol,
                             const std::string CustomAtomProperty) {
  std::vector<double> res;
  RDKit::Descriptors::AUTOCORR2D(mol, res, CustomAtomProperty);
  python::list pyres;
  for (const auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

#endif

RDKit::SparseIntVect<std::int32_t> *GetAtomPairFingerprint(
    const RDKit::ROMol &mol, unsigned int minLength, unsigned int maxLength,
    python::object fromAtoms, python::object ignoreAtoms,
    python::object atomInvariants, bool includeChirality, bool use2D,
    int confId) {
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
    unsigned int maxLength, python::object fromAtoms,
    python::object ignoreAtoms, python::object atomInvariants,
    bool includeChirality, bool use2D, int confId) {
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
    const RDKit::ROMol &mol, unsigned int targetSize, python::object fromAtoms,
    python::object ignoreAtoms, python::object atomInvariants,
    bool includeChirality) {
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
    throw_value_error(errout.str());
  }

  RDKit::SparseIntVect<boost::int64_t> *res;
  res = RDKit::AtomPairs::getTopologicalTorsionFingerprint(
      mol, targetSize, fvect.get(), ivect.get(), invvect.get(),
      includeChirality);
  return res;
}

RDKit::SparseIntVect<boost::int64_t> *GetHashedTopologicalTorsionFingerprint(
    const RDKit::ROMol &mol, unsigned int nBits, unsigned int targetSize,
    python::object fromAtoms, python::object ignoreAtoms,
    python::object atomInvariants, bool includeChirality) {
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
    python::object fromAtoms, python::object ignoreAtoms,
    python::object atomInvariants, unsigned int nBitsPerEntry,
    bool includeChirality) {
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
    unsigned int maxLength, python::object fromAtoms,
    python::object ignoreAtoms, python::object atomInvariants,
    unsigned int nBitsPerEntry, bool includeChirality, bool use2D, int confId) {
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

namespace {
double kappaHelper(double (*fn)(const RDKit::ROMol &, std::vector<double> *),
                   const RDKit::ROMol &mol, python::object atomContribs) {
  std::vector<double> *lContribs = nullptr;
  if (atomContribs != python::object()) {
    // make sure the optional argument actually was a list
    python::list typecheck = python::extract<python::list>(atomContribs);

    if (python::extract<unsigned int>(typecheck.attr("__len__")()) !=
        mol.getNumAtoms()) {
      throw_value_error("length of atomContribs list != number of atoms");
    }

    lContribs = new std::vector<double>(mol.getNumAtoms());
  }
  double res = fn(mol, lContribs);
  if (lContribs) {
    python::list acl = python::extract<python::list>(atomContribs);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      acl[i] = (*lContribs)[i];
    }
    delete lContribs;
  }
  return res;
}
double hkAlphaHelper(const RDKit::ROMol &mol, python::object atomContribs) {
  return kappaHelper(RDKit::Descriptors::calcHallKierAlpha, mol, atomContribs);
}

[[deprecated("please use MorganGenerator")]] RDKit::SparseIntVect<std::uint32_t>
    *MorganFingerprintHelper(const RDKit::ROMol &mol, unsigned int radius,
                             int nBits, python::object invariants,
                             python::object fromAtoms, bool useChirality,
                             bool useBondTypes, bool useFeatures,
                             bool useCounts, python::object bitInfo,
                             bool includeRedundantEnvironments) {
  RDLog::deprecationWarning("please use MorganGenerator");
  std::vector<boost::uint32_t> *invars = nullptr;
  if (invariants) {
    unsigned int nInvar =
        python::extract<unsigned int>(invariants.attr("__len__")());
    if (nInvar) {
      if (nInvar != mol.getNumAtoms()) {
        throw_value_error("length of invariant vector != number of atoms");
      }
      invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        (*invars)[i] = python::extract<std::uint32_t>(invariants[i]);
      }
    }
  } else if (useFeatures) {
    invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol, *invars);
  }
  std::vector<std::uint32_t> *froms = nullptr;
  if (fromAtoms) {
    unsigned int nFrom =
        python::extract<unsigned int>(fromAtoms.attr("__len__")());
    if (nFrom) {
      froms = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < nFrom; ++i) {
        froms->push_back(python::extract<std::uint32_t>(fromAtoms[i]));
      }
    }
  }
  RDKit::MorganFingerprints::BitInfoMap *bitInfoMap = nullptr;
  if (bitInfo != python::object()) {
    // make sure the optional argument actually was a dictionary
    python::dict typecheck = python::extract<python::dict>(bitInfo);
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
    for (RDKit::MorganFingerprints::BitInfoMap::const_iterator iter =
             bitInfoMap->begin();
         iter != bitInfoMap->end(); ++iter) {
      const std::vector<std::pair<std::uint32_t, std::uint32_t>> &v =
          iter->second;
      python::list localL;
      for (const auto &vIt : v) {
        localL.append(python::make_tuple(vIt.first, vIt.second));
      }
      bitInfo[iter->first] = python::tuple(localL);
    }
    delete bitInfoMap;
  }
  if (invars) {
    delete invars;
  }
  if (froms) {
    delete froms;
  }
  return res;
}

#ifdef RDK_HAS_EIGEN3
python::list BCUT(const RDKit::ROMol &mol) {
  return python::list(RDKit::Descriptors::BCUT2D(mol));
}

std::pair<double, double> BCUT2D_list(const RDKit::ROMol &m,
                                      python::list atomprops) {
  std::vector<double> dvec;
  for (int i = 0; i < len(atomprops); ++i) {
    dvec.push_back(boost::python::extract<double>(atomprops[i]));
  }
  return RDKit::Descriptors::BCUT2D(m, dvec);
}

std::pair<double, double> BCUT2D_tuple(const RDKit::ROMol &m,
                                       python::tuple atomprops) {
  std::vector<double> dvec;
  for (int i = 0; i < len(atomprops); ++i) {
    dvec.push_back(boost::python::extract<double>(atomprops[i]));
  }
  return RDKit::Descriptors::BCUT2D(m, dvec);
}

// From boost::python examples
// Converts a std::pair instance to a Python tuple.
template <typename T1, typename T2>
struct std_pair_to_tuple {
  static PyObject *convert(std::pair<T1, T2> const &p) {
    return boost::python::incref(
        boost::python::make_tuple(p.first, p.second).ptr());
  }
  static PyTypeObject const *get_pytype() { return &PyTuple_Type; }
};

// Helper for convenience.
template <typename T1, typename T2>
struct std_pair_to_python_converter {
  std_pair_to_python_converter() {
    boost::python::to_python_converter<std::pair<T1, T2>,
                                       std_pair_to_tuple<T1, T2>,
                                       true  // std_pair_to_tuple has get_pytype
                                       >();
  }
};
#endif
}  // namespace
RDKit::SparseIntVect<std::uint32_t> *GetMorganFingerprint(
    const RDKit::ROMol &mol, unsigned int radius, python::object invariants,
    python::object fromAtoms, bool useChirality, bool useBondTypes,
    bool useFeatures, bool useCounts, python::object bitInfo,
    bool includeRedundantEnvironments) {
  return MorganFingerprintHelper(
      mol, radius, -1, invariants, fromAtoms, useChirality, useBondTypes,
      useFeatures, useCounts, bitInfo, includeRedundantEnvironments);
}
RDKit::SparseIntVect<std::uint32_t> *GetHashedMorganFingerprint(
    const RDKit::ROMol &mol, unsigned int radius, unsigned int nBits,
    python::object invariants, python::object fromAtoms, bool useChirality,
    bool useBondTypes, bool useFeatures, python::object bitInfo,
    bool includeRedundantEnvironments) {
  return MorganFingerprintHelper(mol, radius, nBits, invariants, fromAtoms,
                                 useChirality, useBondTypes, useFeatures, true,
                                 bitInfo, includeRedundantEnvironments);
}

[[deprecated("please use MorganGenerator")]] ExplicitBitVect *
GetMorganFingerprintBV(const RDKit::ROMol &mol, unsigned int radius,
                       unsigned int nBits, python::object invariants,
                       python::object fromAtoms, bool useChirality,
                       bool useBondTypes, bool useFeatures,
                       python::object bitInfo,
                       bool includeRedundantEnvironments) {
  RDLog::deprecationWarning("please use MorganGenerator");
  std::vector<boost::uint32_t> *invars = nullptr;
  if (invariants) {
    unsigned int nInvar =
        python::extract<unsigned int>(invariants.attr("__len__")());
    if (nInvar) {
      if (nInvar != mol.getNumAtoms()) {
        throw_value_error("length of invariant vector != number of atoms");
      }
      invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        (*invars)[i] = python::extract<std::uint32_t>(invariants[i]);
      }
    }
  } else if (useFeatures) {
    invars = new std::vector<std::uint32_t>(mol.getNumAtoms());
    RDKit::MorganFingerprints::getFeatureInvariants(mol, *invars);
  }

  std::unique_ptr<std::vector<std::uint32_t>> froms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  RDKit::MorganFingerprints::BitInfoMap *bitInfoMap = nullptr;
  if (bitInfo != python::object()) {
    // make sure the optional argument actually was a dictionary
    python::dict typecheck = python::extract<python::dict>(bitInfo);
    bitInfoMap = new RDKit::MorganFingerprints::BitInfoMap();
  }
  ExplicitBitVect *res;
  res = RDKit::MorganFingerprints::getFingerprintAsBitVect(
      mol, radius, nBits, invars, froms.get(), useChirality, useBondTypes,
      false, bitInfoMap, includeRedundantEnvironments);
  if (bitInfoMap) {
    bitInfo.attr("clear")();
    for (RDKit::MorganFingerprints::BitInfoMap::const_iterator iter =
             bitInfoMap->begin();
         iter != bitInfoMap->end(); ++iter) {
      const std::vector<std::pair<std::uint32_t, std::uint32_t>> &v =
          iter->second;
      python::list localL;
      for (const auto &vIt : v) {
        localL.append(python::make_tuple(vIt.first, vIt.second));
      }
      bitInfo[iter->first] = python::tuple(localL);
    }
    delete bitInfoMap;
  }
  delete invars;
  return res;
}

python::list GetAtomFeatures(const RDKit::ROMol &mol, int atomid,
                             bool addchiral) {
  std::vector<double> res;
  RDKit::Descriptors::AtomFeatVect(mol, res, atomid, addchiral);
  python::list pyres;
  for (auto iv : res) {
    pyres.append(iv);
  }
  return pyres;
}

python::list GetConnectivityInvariants(const RDKit::ROMol &mol,
                                       bool includeRingMembership) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  RDKit::MorganFingerprints::getConnectivityInvariants(mol, invars,
                                                       includeRingMembership);
  python::list res;
  for (const auto iv : invars) {
    res.append(python::long_(iv));
  }
  return res;
}
python::list GetFeatureInvariants(const RDKit::ROMol &mol) {
  std::vector<std::uint32_t> invars(mol.getNumAtoms());
  RDKit::MorganFingerprints::getFeatureInvariants(mol, invars);
  python::list res;
  for (const auto iv : invars) {
    res.append(python::long_(iv));
  }
  return res;
}

python::list GetUSR(const RDKit::ROMol &mol, int confId) {
  if (mol.getNumConformers() == 0) {
    throw_value_error("no conformers");
  }
  if (mol.getNumAtoms() < 3) {
    throw_value_error("too few atoms (minimum three)");
  }
  std::vector<double> descriptor(12);
  RDKit::Descriptors::USR(mol, descriptor, confId);
  python::list pyDescr;
  for (const auto d : descriptor) {
    pyDescr.append(d);
  }
  return pyDescr;
}

python::list GetUSRDistributions(python::object coords, python::object points) {
  unsigned int numCoords =
      python::extract<unsigned int>(coords.attr("__len__")());
  if (numCoords == 0) {
    throw_value_error("no coordinates");
  }
  RDGeom::Point3DConstPtrVect c(numCoords);
  for (unsigned int i = 0; i < numCoords; ++i) {
    auto *pt = new RDGeom::Point3D;
    *pt = python::extract<RDGeom::Point3D>(coords[i]);
    c[i] = pt;
  }
  std::vector<RDGeom::Point3D> pts(4);
  std::vector<std::vector<double>> distances(4);
  RDKit::Descriptors::calcUSRDistributions(c, distances, pts);
  if (points != python::object()) {
    // make sure the optional argument actually was a list
    python::list tmpPts = python::extract<python::list>(points);
    for (const auto &p : pts) {
      tmpPts.append(p);
    }
    points = tmpPts;
  }
  python::list pyDist;
  for (const auto &dist : distances) {
    python::list pytmp;
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

python::list GetUSRDistributionsFromPoints(python::object coords,
                                           python::object points) {
  unsigned int numCoords =
      python::extract<unsigned int>(coords.attr("__len__")());
  unsigned int numPts = python::extract<unsigned int>(points.attr("__len__")());
  if (numCoords == 0) {
    throw_value_error("no coordinates");
  }
  RDGeom::Point3DConstPtrVect c(numCoords);
  for (unsigned int i = 0; i < numCoords; ++i) {
    auto *pt = new RDGeom::Point3D;
    *pt = python::extract<RDGeom::Point3D>(coords[i]);
    c[i] = pt;
  }
  std::vector<RDGeom::Point3D> p(numPts);
  if (numPts == 0) {
    throw_value_error("no points");
  }
  for (unsigned int i = 0; i < numPts; ++i) {
    p[i] = python::extract<RDGeom::Point3D>(points[i]);
  }
  std::vector<std::vector<double>> distances(numPts);
  RDKit::Descriptors::calcUSRDistributionsFromPoints(c, p, distances);
  python::list pyDist;
  for (const auto &dist : distances) {
    python::list pytmp;
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

python::list GetUSRFromDistributions(python::object distances) {
  unsigned int numDist =
      python::extract<unsigned int>(distances.attr("__len__")());
  if (numDist == 0) {
    throw_value_error("no distances");
  }
  std::vector<std::vector<double>> dist(numDist);
  for (unsigned int i = 0; i < numDist; ++i) {
    unsigned int numPts =
        python::extract<unsigned int>(distances[i].attr("__len__")());
    if (numPts == 0) {
      throw_value_error("distances missing");
    }
    std::vector<double> tmpDist(numPts);
    for (unsigned int j = 0; j < numPts; ++j) {
      tmpDist[j] = python::extract<double>(distances[i][j]);
    }
    dist[i] = tmpDist;
  }
  std::vector<double> descriptor(12);
  RDKit::Descriptors::calcUSRFromDistributions(dist, descriptor);
  python::list pyDescr;
  for (const auto d : descriptor) {
    pyDescr.append(d);
  }
  return pyDescr;
}

double GetUSRScore(python::object descriptor1, python::object descriptor2,
                   python::object weights) {
  unsigned int numElements =
      python::extract<unsigned int>(descriptor1.attr("__len__")());
  if (numElements !=
      python::extract<unsigned int>(descriptor2.attr("__len__")())) {
    throw_value_error("descriptors must have the same length");
  }
  unsigned int numWeights = numElements / 12;
  unsigned int numPyWeights =
      python::extract<unsigned int>(weights.attr("__len__")());
  std::vector<double> w(numWeights, 1.0);  // default weights: all to 1.0
  if ((numPyWeights > 0) && (numPyWeights != numWeights)) {
    throw_value_error("number of weights is not correct");
  } else if (numPyWeights == numWeights) {
    for (unsigned int i = 0; i < numWeights; ++i) {
      w[i] = python::extract<double>(weights[i]);
    }
  }
  std::vector<double> d1(numElements);
  std::vector<double> d2(numElements);
  for (unsigned int i = 0; i < numElements; ++i) {
    d1[i] = python::extract<double>(descriptor1[i]);
    d2[i] = python::extract<double>(descriptor2[i]);
  }
  double res = RDKit::Descriptors::calcUSRScore(d1, d2, w);
  return res;
}

python::list GetUSRCAT(const RDKit::ROMol &mol, python::object atomSelections,
                       int confId) {
  if (mol.getNumConformers() == 0) {
    throw_value_error("no conformers");
  }
  if (mol.getNumAtoms() < 3) {
    throw_value_error("too few atoms (minimum three)");
  }

  // check if there is an atom selection provided
  std::vector<std::vector<unsigned int>> atomIds;
  unsigned int sizeDescriptor = 60;
  if (atomSelections != python::object()) {
    // make sure the optional argument actually was a list
    python::list typecheck = python::extract<python::list>(atomSelections);
    unsigned int numSel =
        python::extract<unsigned int>(atomSelections.attr("__len__")());
    if (numSel == 0) {
      throw_value_error("empty atom selections");
    }
    atomIds.resize(numSel);
    for (unsigned int i = 0; i < numSel; ++i) {
      unsigned int numPts =
          python::extract<unsigned int>(atomSelections[i].attr("__len__")());
      std::vector<unsigned int> tmpIds(numPts);
      for (unsigned int j = 0; j < numPts; ++j) {
        tmpIds[j] = python::extract<unsigned int>(atomSelections[i][j]) - 1;
      }
      atomIds[i] = tmpIds;
    }
    sizeDescriptor = 12 * (numSel + 1);
  }
  std::vector<double> descriptor(sizeDescriptor);
  RDKit::Descriptors::USRCAT(mol, descriptor, atomIds, confId);
  python::list pyDescr;
  for (const auto d : descriptor) {
    pyDescr.append(d);
  }
  return pyDescr;
}

python::list CalcSlogPVSA(const RDKit::ROMol &mol, python::object bins,
                          bool force) {
  std::vector<double> *lbins = nullptr;
  if (bins) {
    unsigned int nBins = python::extract<unsigned int>(bins.attr("__len__")());
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = python::extract<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcSlogP_VSA(mol, lbins, force);

  python::list pyres;
  for (const auto d : res) {
    pyres.append(d);
  }
  return pyres;
}
python::list CalcSMRVSA(const RDKit::ROMol &mol, python::object bins,
                        bool force) {
  std::vector<double> *lbins = nullptr;
  if (bins) {
    unsigned int nBins = python::extract<unsigned int>(bins.attr("__len__")());
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = python::extract<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcSMR_VSA(mol, lbins, force);

  python::list pyres;
  for (const auto d : res) {
    pyres.append(d);
  }
  return pyres;
}
python::list CalcPEOEVSA(const RDKit::ROMol &mol, python::object bins,
                         bool force) {
  std::vector<double> *lbins = nullptr;
  if (bins) {
    unsigned int nBins = python::extract<unsigned int>(bins.attr("__len__")());
    if (nBins) {
      lbins = new std::vector<double>(nBins, 0.0);
      for (unsigned int i = 0; i < nBins; ++i) {
        (*lbins)[i] = python::extract<double>(bins[i]);
      }
    }
  }
  std::vector<double> res;
  res = RDKit::Descriptors::calcPEOE_VSA(mol, lbins, force);

  python::list pyres;
  for (const auto d : res) {
    pyres.append(d);
  }
  return pyres;
}
python::list CalcCustomPropVSA(const RDKit::ROMol &mol,
                               const std::string customPropName,
                               python::object bins, bool force) {
  unsigned int nBins = python::extract<unsigned int>(bins.attr("__len__")());
  std::vector<double> lbins = std::vector<double>(nBins, 0.0);
  for (unsigned int i = 0; i < nBins; ++i) {
    lbins[i] = python::extract<double>(bins[i]);
  }
  std::vector<double> res;
  res =
      RDKit::Descriptors::calcCustomProp_VSA(mol, customPropName, lbins, force);

  python::list pyres;
  for (const auto d : res) {
    pyres.append(d);
  }
  return pyres;
}
python::list CalcMQNs(const RDKit::ROMol &mol, bool force) {
  std::vector<unsigned int> res;
  res = RDKit::Descriptors::calcMQNs(mol, force);

  python::list pyres;
  for (const auto d : res) {
    pyres.append(d);
  }
  return pyres;
}

unsigned int numSpiroAtoms(const RDKit::ROMol &mol, python::object pyatoms) {
  std::vector<unsigned int> ats;
  unsigned int res = RDKit::Descriptors::calcNumSpiroAtoms(
      mol, pyatoms != python::object() ? &ats : nullptr);
  if (pyatoms != python::object()) {
    python::list pyres = python::extract<python::list>(pyatoms);
    for (const auto d : ats) {
      pyres.append(d);
    }
  }
  return res;
}
unsigned int numBridgeheadAtoms(const RDKit::ROMol &mol,
                                python::object pyatoms) {
  std::vector<unsigned int> ats;
  unsigned int res = RDKit::Descriptors::calcNumBridgeheadAtoms(
      mol, pyatoms != python::object() ? &ats : nullptr);
  if (pyatoms != python::object()) {
    python::list pyres = python::extract<python::list>(pyatoms);
    for (const auto d : ats) {
      pyres.append(d);
    }
  }
  return res;
}

struct PythonPropertyFunctor : public RDKit::Descriptors::PropertyFunctor {
  PyObject *self;

  // n.b. until we switch the query d_dataFunc over to boost::function
  //  we can't use python props in functions.
  PythonPropertyFunctor(PyObject *self, const std::string &name,
                        const std::string &version)
      : PropertyFunctor(name, version), self(self) {
    python::incref(self);
  }

  ~PythonPropertyFunctor() override { python::decref(self); }

  double operator()(const RDKit::ROMol &mol) const override {
    return python::call_method<double>(self, "__call__", boost::ref(mol));
  }
};
}  // namespace

BOOST_PYTHON_MODULE(rdMolDescriptors) {
  python::scope().attr("__doc__") =
      "Module containing functions to compute molecular descriptors";

  std::string docString = "";

  python::class_<python::object>("AtomPairsParameters")
      .setattr("version", RDKit::AtomPairs::atomPairsVersion)
      .setattr("numTypeBits", RDKit::AtomPairs::numTypeBits)
      .setattr("numPiBits", RDKit::AtomPairs::numPiBits)
      .setattr("numBranchBits", RDKit::AtomPairs::numBranchBits)
      .setattr("numChiralBits", RDKit::AtomPairs::numChiralBits)
      .setattr("codeSize", RDKit::AtomPairs::codeSize)
      .setattr("atomTypes", atomPairTypes)
      .setattr("numPathBits", RDKit::AtomPairs::numPathBits)
      .setattr("numAtomPairFingerprintBits",
               RDKit::AtomPairs::numAtomPairFingerprintBits);
  docString = "Returns the atom code (hash) for an atom";
  python::def("GetAtomPairAtomCode", RDKit::AtomPairs::getAtomCode,
              (python::arg("atom"), python::arg("branchSubtract") = 0,
               python::arg("includeChirality") = false),
              docString.c_str());
  docString =
      "Returns the atom-pair code (hash) for a pair of atoms separated by a "
      "certain number of bonds";
  python::def(
      "GetAtomPairCode", RDKit::AtomPairs::getAtomPairCode,
      (python::arg("atom1Code"), python::arg("atom2Code"),
       python::arg("distance"), python::arg("includeChirality") = false),
      docString.c_str());
  docString =
      "Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect";
  python::def("GetAtomPairFingerprint", GetAtomPairFingerprint,
              (python::arg("mol"), python::arg("minLength") = 1,
               python::arg("maxLength") = RDKit::AtomPairs::maxPathLen - 1,
               python::arg("fromAtoms") = 0, python::arg("ignoreAtoms") = 0,
               python::arg("atomInvariants") = 0,
               python::arg("includeChirality") = false,
               python::arg("use2D") = true, python::arg("confId") = -1),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the hashed atom-pair fingerprint for a molecule as an "
      "IntSparseIntVect";
  python::def("GetHashedAtomPairFingerprint", GetHashedAtomPairFingerprint,
              (python::arg("mol"), python::arg("nBits") = 2048,
               python::arg("minLength") = 1,
               python::arg("maxLength") = RDKit::AtomPairs::maxPathLen - 1,
               python::arg("fromAtoms") = 0, python::arg("ignoreAtoms") = 0,
               python::arg("atomInvariants") = 0,
               python::arg("includeChirality") = false,
               python::arg("use2D") = true, python::arg("confId") = -1),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect";
  python::def(
      "GetHashedAtomPairFingerprintAsBitVect",
      GetHashedAtomPairFingerprintAsBitVect,
      (python::arg("mol"), python::arg("nBits") = 2048,
       python::arg("minLength") = 1,
       python::arg("maxLength") = RDKit::AtomPairs::maxPathLen - 1,
       python::arg("fromAtoms") = 0, python::arg("ignoreAtoms") = 0,
       python::arg("atomInvariants") = 0, python::arg("nBitsPerEntry") = 4,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("confId") = -1),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the topological-torsion fingerprint for a molecule as a "
      "LongIntSparseIntVect";
  python::def("GetTopologicalTorsionFingerprint",
              GetTopologicalTorsionFingerprint,
              (python::arg("mol"), python::arg("targetSize") = 4,
               python::arg("fromAtoms") = 0, python::arg("ignoreAtoms") = 0,
               python::arg("atomInvariants") = 0,
               python::arg("includeChirality") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Returns the hashed topological-torsion fingerprint for a molecule as a "
      "LongIntSparseIntVect";
  python::def(
      "GetHashedTopologicalTorsionFingerprint",
      GetHashedTopologicalTorsionFingerprint,
      (python::arg("mol"), python::arg("nBits") = 2048,
       python::arg("targetSize") = 4, python::arg("fromAtoms") = 0,
       python::arg("ignoreAtoms") = 0, python::arg("atomInvariants") = 0,
       python::arg("includeChirality") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  docString =
      "Returns the topological-torsion fingerprint for a molecule as an "
      "ExplicitBitVect";
  python::def(
      "GetHashedTopologicalTorsionFingerprintAsBitVect",
      GetHashedTopologicalTorsionFingerprintAsBitVect,
      (python::arg("mol"), python::arg("nBits") = 2048,
       python::arg("targetSize") = 4, python::arg("fromAtoms") = 0,
       python::arg("ignoreAtoms") = 0, python::arg("atomInvariants") = 0,
       python::arg("nBitsPerEntry") = 4,
       python::arg("includeChirality") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  docString = "Returns a Morgan fingerprint for a molecule";
  python::def(
      "GetMorganFingerprint", GetMorganFingerprint,
      (python::arg("mol"), python::arg("radius"),
       python::arg("invariants") = python::list(),
       python::arg("fromAtoms") = python::list(),
       python::arg("useChirality") = false, python::arg("useBondTypes") = true,
       python::arg("useFeatures") = false, python::arg("useCounts") = true,
       python::arg("bitInfo") = python::object(),
       python::arg("includeRedundantEnvironments") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  docString = "Returns a hashed Morgan fingerprint for a molecule";
  python::def(
      "GetHashedMorganFingerprint", GetHashedMorganFingerprint,
      (python::arg("mol"), python::arg("radius"), python::arg("nBits") = 2048,
       python::arg("invariants") = python::list(),
       python::arg("fromAtoms") = python::list(),
       python::arg("useChirality") = false, python::arg("useBondTypes") = true,
       python::arg("useFeatures") = false,
       python::arg("bitInfo") = python::object(),
       python::arg("includeRedundantEnvironments") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  docString = "Returns a Morgan fingerprint for a molecule as a bit vector";
  python::def(
      "GetMorganFingerprintAsBitVect", GetMorganFingerprintBV,
      (python::arg("mol"), python::arg("radius"), python::arg("nBits") = 2048,
       python::arg("invariants") = python::list(),
       python::arg("fromAtoms") = python::list(),
       python::arg("useChirality") = false, python::arg("useBondTypes") = true,
       python::arg("useFeatures") = false,
       python::arg("bitInfo") = python::object(),
       python::arg("includeRedundantEnvironments") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  python::scope().attr("_MorganFingerprint_version") =
      RDKit::MorganFingerprints::morganFingerprintVersion;
  docString = "Returns connectivity invariants (ECFP-like) for a molecule.";
  python::def("GetConnectivityInvariants", GetConnectivityInvariants,
              (python::arg("mol"), python::arg("includeRingMembership") = true),
              docString.c_str());
  python::scope().attr("_ConnectivityInvariants_version") =
      RDKit::MorganFingerprints::morganConnectivityInvariantVersion;

  docString = "Returns feature invariants (FCFP-like) for a molecule.";
  python::def("GetFeatureInvariants", GetFeatureInvariants,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_FeatureInvariants_version") =
      RDKit::MorganFingerprints::morganFeatureInvariantVersion;

  // USR descriptor
  docString = "Returns a USR descriptor for one conformer of a molecule";
  python::def("GetUSR", GetUSR,
              (python::arg("mol"), python::arg("confId") = -1),
              docString.c_str());
  docString =
      "Returns the four USR distance distributions for a set of coordinates";
  python::def("GetUSRDistributions", GetUSRDistributions,
              (python::arg("coords"), python::arg("points") = python::object()),
              docString.c_str());
  docString =
      "Returns the USR distance distributions for a set of coordinates and "
      "points";
  python::def("GetUSRDistributionsFromPoints", GetUSRDistributionsFromPoints,
              (python::arg("coords"), python::arg("points")),
              docString.c_str());
  docString = "Returns the USR descriptor from a set of distance distributions";
  python::def("GetUSRFromDistributions", GetUSRFromDistributions,
              (python::arg("distances")), docString.c_str());
  docString = "Returns the USR score for two USR or USRCAT descriptors";
  python::def("GetUSRScore", GetUSRScore,
              (python::arg("descriptor1"), python::arg("descriptor2"),
               python::arg("weights") = python::list()),
              docString.c_str());
  docString = "Returns a USRCAT descriptor for one conformer of a molecule";
  python::def(
      "GetUSRCAT", GetUSRCAT,
      (python::arg("mol"), python::arg("atomSelections") = python::object(),
       python::arg("confId") = -1),
      docString.c_str());

  docString =
      "returns (as a list of 2-tuples) the contributions of each atom to\n"
      "the Wildman-Cripppen logp and mr value";
  python::def("_CalcCrippenContribs", computeCrippenContribs,
              (python::arg("mol"), python::arg("force") = false,
               python::arg("atomTypes") = python::list(),
               python::arg("atomTypeLabels") = python::list()),
              docString.c_str());
  docString = "returns a 2-tuple with the Wildman-Crippen logp,mr values";
  python::def("CalcCrippenDescriptors", calcCrippenDescriptors,
              (python::arg("mol"), python::arg("includeHs") = true,
               python::arg("force") = false),
              docString.c_str());
  python::scope().attr("_CalcCrippenDescriptors_version") =
      RDKit::Descriptors::crippenVersion;

  docString = "returns the Labute ASA value for a molecule";
  python::def("CalcLabuteASA", RDKit::Descriptors::calcLabuteASA,
              (python::arg("mol"), python::arg("includeHs") = true,
               python::arg("force") = false),
              docString.c_str());
  python::scope().attr("_CalcLabuteASA_version") =
      RDKit::Descriptors::labuteASAVersion;

  docString = "returns a list of atomic contributions to the Labute ASA";
  python::def("_CalcLabuteASAContribs", computeASAContribs,
              (python::arg("mol"), python::arg("includeHs") = true,
               python::arg("force") = false),
              docString.c_str());

  docString = "returns the TPSA value for a molecule";
  python::def("CalcTPSA", RDKit::Descriptors::calcTPSA,
              (python::arg("mol"), python::arg("force") = false,
               python::arg("includeSandP") = false),
              docString.c_str());
  python::scope().attr("_CalcTPSA_version") = RDKit::Descriptors::tpsaVersion;

  docString = "returns a list of atomic contributions to the TPSA";
  python::def("_CalcTPSAContribs", computeTPSAContribs,
              (python::arg("mol"), python::arg("force") = false,
               python::arg("includeSandP") = false),
              docString.c_str());

  docString = "returns the molecule's molecular weight";
  python::def("_CalcMolWt", RDKit::Descriptors::calcAMW,
              (python::arg("mol"), python::arg("onlyHeavy") = false),
              docString.c_str());
  python::scope().attr("_CalcMolWt_version") = "1.0.0";

  docString = "returns the molecule's exact molecular weight";
  python::def("CalcExactMolWt", RDKit::Descriptors::calcExactMW,
              (python::arg("mol"), python::arg("onlyHeavy") = false),
              docString.c_str());
  python::scope().attr("_CalcExactMolWt_version") = "1.0.0";

  docString = "returns the molecule's formula";
  python::def("CalcMolFormula", RDKit::Descriptors::calcMolFormula,
              (python::arg("mol"), python::arg("separateIsotopes") = false,
               python::arg("abbreviateHIsotopes") = true),
              docString.c_str());
  python::scope().attr("_CalcMolFormula_version") = "1.3.0";

  docString = "returns the number of Lipinski H-bond donors for a molecule";
  python::def("CalcNumLipinskiHBD", RDKit::Descriptors::calcLipinskiHBD,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumLipinskiHBD_version") =
      RDKit::Descriptors::lipinskiHBDVersion;
  docString = "returns the number of Lipinski H-bond acceptors for a molecule";
  python::def("CalcNumLipinskiHBA", RDKit::Descriptors::calcLipinskiHBA,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumLipinskiHBA_version") =
      RDKit::Descriptors::lipinskiHBAVersion;
  docString = "returns the number of H-bond donors for a molecule";
  python::def("CalcNumHBD", RDKit::Descriptors::calcNumHBD,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumHBD_version") =
      RDKit::Descriptors::NumHBDVersion;
  docString = "returns the number of H-bond acceptors for a molecule";
  python::def("CalcNumHBA", RDKit::Descriptors::calcNumHBA,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumHBA_version") =
      RDKit::Descriptors::NumHBAVersion;

  // exposes calcNumRotatableBondOptions (must be a better way!)

  docString =
      "Options for generating rotatable bonds\n\
 NonStrict - standard loose definitions\n\
 Strict - stricter definition excluding amides, esters, etc\n\
 StrictLinkages - adds rotors between rotatable bonds\n\
 Default - Current RDKit default\n";

  python::enum_<RDKit::Descriptors::NumRotatableBondsOptions>(
      "NumRotatableBondsOptions", docString.c_str())
      .value("NonStrict", RDKit::Descriptors::NonStrict)
      .value("Strict", RDKit::Descriptors::Strict)
      .value("StrictLinkages", RDKit::Descriptors::StrictLinkages)
      .value("Default", RDKit::Descriptors::Default);

#ifdef RDK_USE_STRICT_ROTOR_DEFINITION
  docString =
      "returns the number of rotatable bonds for a molecule.\n\
   strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.\n\
   strict = NumRotatableBondsOptions.Strict - (default) does not count things like\n\
            amide or ester bonds\n\
   strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring\n\
      systems.\n\
      - Single bonds between aliphatic ring Cs are always rotatable. This\n\
        means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now \n\
        considered rotatable; it was not before\n\
      - Heteroatoms in the linked rings no longer affect whether or not\n\
        the linking bond is rotatable\n\
      - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now\n\
         considered non-rotatable";
#else
  docString =
      "returns the number of rotatable bonds for a molecule.\n\
   strict = NumRotatableBondsOptions.NonStrict - (default) Simple rotatable bond definition.\n\
   strict = NumRotatableBondsOptions.Strict - does not count things like\n\
            amide or ester bonds\n\
   strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring\n\
      systems.\n\
      - Single bonds between aliphatic ring Cs are always rotatable. This\n\
        means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now \n\
        considered rotatable; it was not before\n\
      - Heteroatoms in the linked rings no longer affect whether or not\n\
        the linking bond is rotatable\n\
      - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now\n\
         considered non-rotatable";
#endif

  python::def("CalcNumRotatableBonds",
              (unsigned int (*)(const RDKit::ROMol &,
                                bool))RDKit::Descriptors::calcNumRotatableBonds,
              (python::arg("mol"), python::arg("strict")), docString.c_str());

  python::def(
      "CalcNumRotatableBonds",
      (unsigned int (*)(const RDKit::ROMol &,
                        RDKit::Descriptors::NumRotatableBondsOptions))
          RDKit::Descriptors::calcNumRotatableBonds,
      (python::arg("mol"), python::arg("strict") = RDKit::Descriptors::Default),
      docString.c_str());
  python::scope().attr("_CalcNumRotatableBonds_version") =
      RDKit::Descriptors::NumRotatableBondsVersion;

  docString = "returns the number of rings for a molecule";
  python::def("CalcNumRings", RDKit::Descriptors::calcNumRings,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumRings_version") =
      RDKit::Descriptors::NumRingsVersion;

  docString = "returns the number of aromatic rings for a molecule";
  python::def("CalcNumAromaticRings", RDKit::Descriptors::calcNumAromaticRings,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAromaticRings_version") =
      RDKit::Descriptors::NumAromaticRingsVersion;

  docString = "returns the number of saturated rings for a molecule";
  python::def("CalcNumSaturatedRings",
              RDKit::Descriptors::calcNumSaturatedRings, (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumSaturatedRings_version") =
      RDKit::Descriptors::NumSaturatedRingsVersion;

  docString = "returns the number of heterocycles for a molecule";
  python::def("CalcNumHeterocycles", RDKit::Descriptors::calcNumHeterocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumHeterocycles_version") =
      RDKit::Descriptors::NumHeterocyclesVersion;

  docString = "returns the number of aromatic heterocycles for a molecule";
  python::def("CalcNumAromaticHeterocycles",
              RDKit::Descriptors::calcNumAromaticHeterocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAromaticHeterocycles_version") =
      RDKit::Descriptors::NumAromaticHeterocyclesVersion;

  docString = "returns the number of aromatic carbocycles for a molecule";
  python::def("CalcNumAromaticCarbocycles",
              RDKit::Descriptors::calcNumAromaticCarbocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAromaticCarbocycles_version") =
      RDKit::Descriptors::NumAromaticCarbocyclesVersion;

  docString = "returns the number of saturated heterocycles for a molecule";
  python::def("CalcNumSaturatedHeterocycles",
              RDKit::Descriptors::calcNumSaturatedHeterocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumSaturatedHeterocycles_version") =
      RDKit::Descriptors::NumSaturatedHeterocyclesVersion;

  docString = "returns the number of saturated carbocycles for a molecule";
  python::def("CalcNumSaturatedCarbocycles",
              RDKit::Descriptors::calcNumSaturatedCarbocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumSaturatedCarbocycles_version") =
      RDKit::Descriptors::NumSaturatedCarbocyclesVersion;

  docString =
      "returns the number of aliphatic (containing at least one non-aromatic "
      "bond) rings for a molecule";
  python::def("CalcNumAliphaticRings",
              RDKit::Descriptors::calcNumAliphaticRings, (python::arg("mol")),
              docString.c_str());
  python::scope().attr("_CalcNumAliphaticRings_version") =
      RDKit::Descriptors::NumAliphaticRingsVersion;

  docString =
      "returns the number of aliphatic (containing at least one non-aromatic "
      "bond) heterocycles for a molecule";
  python::def("CalcNumAliphaticHeterocycles",
              RDKit::Descriptors::calcNumAliphaticHeterocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAliphaticHeterocycles_version") =
      RDKit::Descriptors::NumAliphaticHeterocyclesVersion;

  docString =
      "returns the number of aliphatic (containing at least one non-aromatic "
      "bond) carbocycles for a molecule";
  python::def("CalcNumAliphaticCarbocycles",
              RDKit::Descriptors::calcNumAliphaticCarbocycles,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAliphaticCarbocycles_version") =
      RDKit::Descriptors::NumAliphaticCarbocyclesVersion;

  docString = "returns the number of heavy atoms for a molecule";
  python::def("CalcNumHeavyAtoms", RDKit::Descriptors::calcNumHeavyAtoms,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumHeavyAtoms_version") =
      RDKit::Descriptors::NumHeavyAtomsVersion;
  docString = "returns the total number of atoms for a molecule";
  python::def("CalcNumAtoms", RDKit::Descriptors::calcNumAtoms,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAtoms_version") =
      RDKit::Descriptors::NumAtomsVersion;
  docString = "returns the number of heteroatoms for a molecule";
  python::def("CalcNumHeteroatoms", RDKit::Descriptors::calcNumHeteroatoms,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumHeteroatoms_version") =
      RDKit::Descriptors::NumHeteroatomsVersion;
  docString = "returns the number of amide bonds in a molecule";
  python::def("CalcNumAmideBonds", RDKit::Descriptors::calcNumAmideBonds,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcNumAmideBonds_version") =
      RDKit::Descriptors::NumAmideBondsVersion;

  docString = "returns the fraction of C atoms that are SP3 hybridized";
  python::def("CalcFractionCSP3", RDKit::Descriptors::calcFractionCSP3,
              (python::arg("mol")), docString.c_str());
  python::scope().attr("_CalcFractionCSP3_version") =
      RDKit::Descriptors::FractionCSP3Version;

  docString = "returns the SlogP VSA contributions for a molecule";
  python::def("SlogP_VSA_", CalcSlogPVSA,
              (python::arg("mol"), python::arg("bins") = python::list(),
               python::arg("force") = false));
  docString = "returns the SMR VSA contributions for a molecule";
  python::def("SMR_VSA_", CalcSMRVSA,
              (python::arg("mol"), python::arg("bins") = python::list(),
               python::arg("force") = false));
  docString = "returns the PEOE VSA contributions for a molecule";
  python::def("PEOE_VSA_", CalcPEOEVSA,
              (python::arg("mol"), python::arg("bins") = python::list(),
               python::arg("force") = false));
  docString =
      "returns the VSA contributions based on a custom property for a molecule";
  python::def("CustomProp_VSA_", CalcCustomPropVSA,
              (python::arg("mol"), python::arg("customPropName"),
               python::arg("bins"), python::arg("force") = false));
  docString = "returns the MQN descriptors for a molecule";
  python::def("MQNs_", CalcMQNs,
              (python::arg("mol"), python::arg("force") = false));

  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def(
      "CalcChiNv", RDKit::Descriptors::calcChiNv,
      (python::arg("mol"), python::arg("n"), python::arg("force") = false));
  python::scope().attr("_CalcChiNv_version") = RDKit::Descriptors::chiNvVersion;
  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def("CalcChi0v", RDKit::Descriptors::calcChi0v,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi0v_version") = RDKit::Descriptors::chi0vVersion;
  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def("CalcChi1v", RDKit::Descriptors::calcChi1v,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi1v_version") = RDKit::Descriptors::chi1vVersion;
  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def("CalcChi2v", RDKit::Descriptors::calcChi2v,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi2v_version") = RDKit::Descriptors::chi2vVersion;
  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def("CalcChi3v", RDKit::Descriptors::calcChi3v,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi3v_version") = RDKit::Descriptors::chi3vVersion;
  docString =
      "From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991)";
  python::def("CalcChi4v", RDKit::Descriptors::calcChi4v,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi4v_version") = RDKit::Descriptors::chi4vVersion;
  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def(
      "CalcChiNn", RDKit::Descriptors::calcChiNn,
      (python::arg("mol"), python::arg("n"), python::arg("force") = false));
  python::scope().attr("_CalcChiNn_version") = RDKit::Descriptors::chiNnVersion;

  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def("CalcChi0n", RDKit::Descriptors::calcChi0n,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi0n_version") = RDKit::Descriptors::chi0nVersion;

  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def("CalcChi1n", RDKit::Descriptors::calcChi1n,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi1n_version") = RDKit::Descriptors::chi1nVersion;

  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def("CalcChi2n", RDKit::Descriptors::calcChi2n,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi2n_version") = RDKit::Descriptors::chi2nVersion;

  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def("CalcChi3n", RDKit::Descriptors::calcChi3n,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi3n_version") = RDKit::Descriptors::chi3nVersion;

  docString =
      "Similar to ChiXv, but uses uses nVal instead of valence. This makes a "
      "big difference after we get out of the first row.";
  python::def("CalcChi4n", RDKit::Descriptors::calcChi4n,
              (python::arg("mol"), python::arg("force") = false));
  python::scope().attr("_CalcChi4n_version") = RDKit::Descriptors::chi4nVersion;

  docString =
      "From equation (58) of Rev. Comp. Chem. vol 2, 367-422, (1991). "
      "NOTE: Because hybridization is used to calculate this, results may "
      "differ from other implementations which have different conventions for "
      "assigning hybridization";
  python::def(
      "CalcHallKierAlpha", hkAlphaHelper,
      (python::arg("mol"), python::arg("atomContribs") = python::object()));
  python::scope().attr("_CalcHallKierAlpha_version") =
      RDKit::Descriptors::hallKierAlphaVersion;

  docString =
      "From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991) "
      "NOTE: Because hybridization is used to calculate this, results may "
      "differ from other implementations which have different conventions for "
      "assigning hybridization";
  python::def("CalcKappa1", RDKit::Descriptors::calcKappa1,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa1_version") =
      RDKit::Descriptors::kappa1Version;
  docString =
      "From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991) "
      "NOTE: Because hybridization is used to calculate this, results may "
      "differ from other implementations which have different conventions for "
      "assigning hybridization";
  python::def("CalcKappa2", RDKit::Descriptors::calcKappa2,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa2_version") =
      RDKit::Descriptors::kappa2Version;
  docString =
      "From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, "
      "(1991) "
      "NOTE: Because hybridization is used to calculate this, results may "
      "differ from other implementations which have different conventions for "
      "assigning hybridization";
  python::def("CalcKappa3", RDKit::Descriptors::calcKappa3,
              (python::arg("mol")));
  python::scope().attr("_CalcKappa3_version") =
      RDKit::Descriptors::kappa3Version;
  docString =
      "From Quantitative Structure-Activity Relationships 8, 221–224 (1989). "
      "NOTE: Because hybridization is used to calculate this, results may "
      "differ from other implementations which have different conventions for "
      "assigning hybridization";
  python::def("CalcPhi", RDKit::Descriptors::calcPhi, (python::arg("mol")));
  python::scope().attr("_CalcPhi_version") = RDKit::Descriptors::PhiVersion;

  docString = "Returns the MACCS keys for a molecule as an ExplicitBitVect";
  python::def("GetMACCSKeysFingerprint",
              RDKit::MACCSFingerprints::getFingerprintAsBitVect,
              (python::arg("mol")), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  python::scope().attr("_GetAtomFeatures_version") =
      RDKit::Descriptors::AtomFeatVersion;
  docString = "Returns the Atom Features vector";
  python::def("GetAtomFeatures", GetAtomFeatures,
              (python::arg("mol"), python::arg("atomid"),
               python::arg("addchiral") = false),
              docString.c_str());

  python::scope().attr("_CalcNumSpiroAtoms_version") =
      RDKit::Descriptors::NumSpiroAtomsVersion;
  docString =
      "Returns the number of spiro atoms (atoms shared between rings that "
      "share exactly one atom)";
  python::def("CalcNumSpiroAtoms", numSpiroAtoms,
              (python::arg("mol"), python::arg("atoms") = python::object()),
              docString.c_str());

  python::scope().attr("_CalcNumBridgeheadAtoms_version") =
      RDKit::Descriptors::NumBridgeheadAtomsVersion;
  docString =
      "Returns the number of bridgehead atoms (atoms shared between rings that "
      "share at least two bonds)";
  python::def("CalcNumBridgeheadAtoms", numBridgeheadAtoms,
              (python::arg("mol"), python::arg("atoms") = python::object()),
              docString.c_str());

  python::scope().attr("_CalcNumAtomStereoCenters_version") =
      RDKit::Descriptors::NumAtomStereoCentersVersion;
  docString =
      "Returns the total number of atomic stereocenters (specified and "
      "unspecified)";
  python::def("CalcNumAtomStereoCenters",
              RDKit::Descriptors::numAtomStereoCenters, (python::arg("mol")),
              docString.c_str());

  python::scope().attr("_CalcNumUnspecifiedAtomStereoCenters_version") =
      RDKit::Descriptors::NumUnspecifiedAtomStereoCentersVersion;
  docString = "Returns the number of unspecified atomic stereocenters";
  python::def("CalcNumUnspecifiedAtomStereoCenters",
              RDKit::Descriptors::numUnspecifiedAtomStereoCenters,
              (python::arg("mol")), docString.c_str());

  docString =
      "Property computation class stored in the property registry.\n"
      "See rdkit.Chem.rdMolDescriptor.Properties.GetProperty and \n"
      "rdkit.Chem.Descriptor.Properties.PropertyFunctor for creating new ones";
  python::class_<RDKit::Descriptors::PropertyFunctor,
                 boost::shared_ptr<RDKit::Descriptors::PropertyFunctor>,
                 boost::noncopyable>("PropertyFunctor", docString.c_str(),
                                     python::no_init)
      .def("__call__", &RDKit::Descriptors::PropertyFunctor::operator(),
           python::args("self", "mol"),
           "Compute the property for the specified molecule")
      .def("GetName", &RDKit::Descriptors::PropertyFunctor::getName,
           python::args("self"), "Return the name of the property to calculate")
      .def("GetVersion", &RDKit::Descriptors::PropertyFunctor::getVersion,
           python::args("self"),
           "Return the version of the calculated property");

  iterable_converter().from_python<std::vector<std::string>>();

  docString =
      "Property computation and registry system.  To compute all registered "
      "properties:\n"
      "mol = Chem.MolFromSmiles('c1ccccc1')\n"
      "properties = rdMolDescriptors.Properties()\n"
      "for name, value in zip(properties.GetPropertyNames(), "
      "properties.ComputeProperties(mol)):\n"
      "  print(name, value)\n\n"
      "To compute a subset\n"
      "properties = rdMolDescriptors.Properties(['exactmw', 'lipinskiHBA'])\n"
      "for name, value in zip(properties.GetPropertyNames(), "
      "properties.ComputeProperties(mol)):\n"
      "  print(name, value)\n\n"
      "";

  python::class_<RDKit::Descriptors::Properties,
                 RDKit::Descriptors::Properties *>(
      "Properties", docString.c_str(), python::init<>(python::args("self")))
      .def(python::init<const std::vector<std::string> &>(
          python::args("self", "propNames")))
      .def("GetPropertyNames",
           &RDKit::Descriptors::Properties::getPropertyNames,
           python::args("self"),
           "Return the property names computed by this instance")
      .def("ComputeProperties",
           &RDKit::Descriptors::Properties::computeProperties,
           ((python::arg("self"), python::arg("mol")),
            python::arg("annotateMol") = false),
           "Return a list of computed properties, if annotateMol==True, "
           "annotate the molecule with "
           "the computed properties.")
      .def("AnnotateProperties",
           &RDKit::Descriptors::Properties::annotateProperties,
           (python::arg("self"), python::arg("mol")),
           "Annotate the molecule with the computed properties.  These "
           "properties will be available "
           "as SDData or from mol.GetProp(prop)")
      .def("GetAvailableProperties",
           &RDKit::Descriptors::Properties::getAvailableProperties,
           "Return all available property names that can be computed")
      .staticmethod("GetAvailableProperties")
      .def("GetProperty", &RDKit::Descriptors::Properties::getProperty,
           python::arg("propName"), "Return the named property if it exists")
      .staticmethod("GetProperty")
      .def("RegisterProperty",
           &RDKit::Descriptors::Properties::registerProperty,
           python::arg("propertyFunctor"),
           "Register a new property object (not thread safe)")
      .staticmethod("RegisterProperty");

  python::class_<PythonPropertyFunctor, boost::noncopyable,
                 python::bases<RDKit::Descriptors::PropertyFunctor>>(
      "PythonPropertyFunctor", "",
      python::init<PyObject *, const std::string &, const std::string &>(
          python::args("self", "callback", "name", "version")))
      .def("__call__", &PythonPropertyFunctor::operator(),
           python::args("self", "mol"),
           "Compute the property for the specified molecule");

  docString =
      "Property Range Query for a molecule.  Match(mol) -> true if in range";
  python::class_<Queries::RangeQuery<double, RDKit::ROMol const &, true>,
                 Queries::RangeQuery<double, RDKit::ROMol const &, true> *,
                 boost::noncopyable>("PropertyRangeQuery", docString.c_str(),
                                     python::no_init)
      .def("Match",
           &Queries::RangeQuery<double, RDKit::ROMol const &, true>::Match,
           python::args("self", "what"));

  docString =
      "Generates a Range property for the specified property, between min and "
      "max\n"
      "query = MakePropertyRangeQuery('exactmw', 0, 500)\n"
      "query.Match( mol )";

  python::def("MakePropertyRangeQuery",
              RDKit::Descriptors::makePropertyRangeQuery,
              (python::arg("name"), python::arg("min"), python::arg("max")),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      R"DOC(ARGUMENTS:
      "   - mol: molecule or protein under consideration
      "   - isProtein: flag to indicate if the input is a protein (default=False, free ligand).
      "   - includeLigand: flag to include or exclude a bound ligand when input is a protein (default=True)
      "   - probeRadius: radius of the solvent probe (default=1.2)
      "   - depth: control of number of dots per atom (default=4)
      "   - dotDensity: control of accuracy (default=0)
      ")DOC";
  python::class_<RDKit::Descriptors::DoubleCubicLatticeVolume>(
      "DoubleCubicLatticeVolume",
      "Class for the Double Cubic Lattice Volume method",
      python::init<const RDKit::ROMol &,
                   python::optional<bool, bool, double, int, int>>(
          (python::args("self", "mol"), python::args("isProtein") = false,
           python::args("includeLigand") = true,
           python::args("probeRadius") = 1.2, python::args("depth") = 4,
           python::args("dotDensity") = 0),
          docString.c_str()))
      .def("GetSurfaceArea",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getSurfaceArea,
           "Get the Surface Area of the Molecule or Protein")
      .def("GetVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getVolume,
           "Get the Total Volume of the Molecule or Protein")
      .def("GetVDWVolume",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getVDWVolume,
           "Get the van der Waals Volume of the Molecule or Protein")
      .def("GetCompactness",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getCompactness,
           "Get the Compactness of the Protein")
      .def("GetPackingDensity",
           &RDKit::Descriptors::DoubleCubicLatticeVolume::getPackingDensity,
           "Get the PackingDensity of the Protein");

#ifdef RDK_BUILD_DESCRIPTORS3D
  python::scope().attr("_CalcCoulombMat_version") =
      RDKit::Descriptors::CoulombMatVersion;
  docString = "Returns severals Coulomb randomized matrices";
  python::def("CalcCoulombMat", calcCoulombMat,
              (python::arg("mol"), python::arg("confId") = -1),
              docString.c_str());

  python::scope().attr("_CalcEMMcharges_version") =
      RDKit::Descriptors::EEMVersion;
  docString = "Returns EEM atomic partial charges";
  python::def("CalcEEMcharges", calcEEMcharges,
              (python::arg("mol"), python::arg("confId") = -1),
              docString.c_str());

  python::scope().attr("_CalcWHIM_version") = RDKit::Descriptors::WHIMVersion;
  docString = "Returns the WHIM descriptors vector";
  python::def(
      "CalcWHIM", calcWHIMs,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("thresh") = 0.001, python::arg("CustomAtomProperty") = ""),
      docString.c_str());

  python::scope().attr("_CalcGETAWAY_version") =
      RDKit::Descriptors::GETAWAYVersion;
  docString = "Returns the GETAWAY descriptors vector";
  python::def(
      "CalcGETAWAY", calcGETAWAYs,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("precision") = 2, python::arg("CustomAtomProperty") = ""),
      docString.c_str());

  python::scope().attr("_CalcRDF_version") = RDKit::Descriptors::RDFVersion;
  docString = "Returns radial distribution fonction descriptors (RDF)";
  python::def("CalcRDF", calcRDFs,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("CustomAtomProperty") = ""),
              docString.c_str());

  python::scope().attr("_CalcMORSE_version") = RDKit::Descriptors::MORSEVersion;
  docString =
      "Returns Molecule Representation of Structures based on Electron "
      "diffraction descriptors";
  python::def("CalcMORSE", calcMORSEs,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("CustomAtomProperty") = ""),
              docString.c_str());

  python::scope().attr("_CalcAUTOCORR3D_version") =
      RDKit::Descriptors::AUTOCORR3DVersion;
  docString = "Returns 3D Autocorrelation descriptors vector";
  python::def("CalcAUTOCORR3D", calcAUTOCORR3Ds,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("CustomAtomProperty") = ""),
              docString.c_str());

  python::scope().attr("_CalcPBF_version") = RDKit::Descriptors::PBFVersion;
  docString =
      "Returns the PBF (plane of best fit) descriptor "
      "(https://doi.org/10.1021/ci300293f)";
  python::def("CalcPBF", RDKit::Descriptors::PBF,
              (python::arg("mol"), python::arg("confId") = -1),
              docString.c_str());
  python::scope().attr("_CalcNPR1_version") = RDKit::Descriptors::NPR1Version;
  docString = "";
  python::def(
      "CalcNPR1", RDKit::Descriptors::NPR1,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcNPR2_version") = RDKit::Descriptors::NPR2Version;
  docString = "";
  python::def(
      "CalcNPR2", RDKit::Descriptors::NPR2,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcPMI1_version") = RDKit::Descriptors::PMI1Version;
  docString = "";
  python::def(
      "CalcPMI1", RDKit::Descriptors::PMI1,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcPMI2_version") = RDKit::Descriptors::PMI2Version;
  docString = "";
  python::def(
      "CalcPMI2", RDKit::Descriptors::PMI2,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcPMI3_version") = RDKit::Descriptors::PMI3Version;
  docString = "";
  python::def(
      "CalcPMI3", RDKit::Descriptors::PMI3,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());

  python::scope().attr("_CalcRadiusOfGyration_version") =
      RDKit::Descriptors::radiusOfGyrationVersion;
  docString = "";
  python::def(
      "CalcRadiusOfGyration", RDKit::Descriptors::radiusOfGyration,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcInertialShapeFactor_version") =
      RDKit::Descriptors::inertialShapeFactorVersion;
  docString = "";
  python::def(
      "CalcInertialShapeFactor", RDKit::Descriptors::inertialShapeFactor,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());

  python::scope().attr("_CalcEccentricity_version") =
      RDKit::Descriptors::eccentricityVersion;
  docString = "";
  python::def(
      "CalcEccentricity", RDKit::Descriptors::eccentricity,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcAsphericity_version") =
      RDKit::Descriptors::asphericityVersion;
  docString = "";
  python::def(
      "CalcAsphericity", RDKit::Descriptors::asphericity,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("useAtomicMasses") = true, python::arg("force") = true),
      docString.c_str());
  python::scope().attr("_CalcSpherocityIndex_version") =
      RDKit::Descriptors::spherocityIndexVersion;
  docString = "";
  python::def("CalcSpherocityIndex", RDKit::Descriptors::spherocityIndex,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("force") = true),
              docString.c_str());

  python::scope().attr("_CalcAUTOCORR2D_version") =
      RDKit::Descriptors::AUTOCORR2DVersion;
  docString = "Returns 2D Autocorrelation descriptors vector";
  python::def("CalcAUTOCORR2D", calcAUTOCORR2Ds,
              (python::arg("mol"), python::arg("CustomAtomProperty") = ""),
              docString.c_str());

#endif

#ifdef RDK_HAS_EIGEN3
  python::scope().attr("_BCUT2D_version") = RDKit::Descriptors::BCUT2DVersion;
  std::pair<double, double> (*BCUT_atomprops)(
      const RDKit::ROMol &, const std::string &) = &RDKit::Descriptors::BCUT2D;
  docString =
      "Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, "
      "No. 1, 1999"
      "Diagonal elements are (currently) atomic mass, gasteiger charge,"
      "crippen logP and crippen MRReturns the 2D BCUT2D descriptors vector as "
      "described in\n"
      "returns [mass eigen value high, mass eigen value low,\n"
      "         gasteiger charge eigenvalue high, gasteiger charge low,\n"
      "         crippen lowgp  eigenvalue high, crippen lowgp  low,\n"
      "         crippen mr eigenvalue high, crippen mr low]\n"
      "";

  python::def("BCUT2D", BCUT, (python::arg("mol")), docString.c_str());

  std_pair_to_python_converter<double, double>();
  docString =
      "Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule "
      "and the specified atom props\n"
      " atom_props must be a list or tuple of floats equal in "
      "size to the number of atoms in mol";

  python::def("BCUT2D", BCUT2D_list,
              (python::arg("mol"), python::arg("atom_props")),
              docString.c_str());
  python::def("BCUT2D", BCUT2D_tuple,
              (python::arg("mol"), python::arg("atom_props")),
              docString.c_str());

  docString =
      "Returns a 2D BCUT (eigen value high, eigen value low) given the "
      "molecule and the specified atom prop name\n"
      "atom_propname must exist on each atom and be convertible to a float";
  python::def("BCUT2D", BCUT_atomprops,
              (python::arg("mol"), python::arg("atom_propname")),
              docString.c_str());

  docString =
      "Adds the oxidation number/state to the atoms of a molecule as"
      " property OxidationNumber on each atom.  Use Pauling"
      " electronegativities."
      "  This is experimental code, still under development.";
  python::def("CalcOxidationNumbers", RDKit::Descriptors::calcOxidationNumbers,
              (python::arg("mol")), docString.c_str());
#endif
}
