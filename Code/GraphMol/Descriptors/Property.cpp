//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "Property.h"
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/Descriptors/BCUT.h>
#include <RDGeneral/RDThreads.h>

#include <algorithm>
#include <exception>
#include <atomic>
#include <limits>
#include <map>
#include <new>
#include <unordered_map>
#include <utility>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#include <mutex>
#endif

namespace RDKit {
namespace Descriptors {

namespace {
void registerFunc(const char *name, const std::string &version,
                  double (*func)(const ROMol &)) {
  Properties::registerProperty(new PropertyFunctor(name, version, func));
}

std::vector<std::string> numberedNames(const std::string &prefix,
                                       unsigned int count) {
  std::vector<std::string> res;
  for (unsigned int i = 1; i <= count; ++i) {
    res.push_back(prefix + std::to_string(i));
  }
  return res;
}

#ifdef RDK_HAS_EIGEN3
std::vector<double> bcut2D(const ROMol &m) { return BCUT2D(m); }
#endif
std::vector<double> peoeVSA(const ROMol &m) { return calcPEOE_VSA(m); }
std::vector<double> smrVSA(const ROMol &m) { return calcSMR_VSA(m); }
std::vector<double> slogpVSA(const ROMol &m) { return calcSlogP_VSA(m); }
std::vector<double> estateVSA(const ROMol &m) { return calcEState_VSA(m); }
std::vector<double> vsaEState(const ROMol &m) { return calcVSA_EState(m); }
std::vector<double> fragments(const ROMol &m) {
  auto counts = calcFragmentDescriptors(m);
  return std::vector<double>(counts.begin(), counts.end());
}
constexpr size_t numFragments = 85;

template <VectorElementPropertyFunctor::VectorFunc F, size_t I>
double vectorElement(const ROMol &m) {
  return F(m).at(I);
}

// registers element I of F as names[I], for each I
template <VectorElementPropertyFunctor::VectorFunc F, size_t... I>
void registerVector(const std::vector<std::string> &names,
                    const std::string &version, std::index_sequence<I...>) {
  (Properties::registerProperty(new VectorElementPropertyFunctor(
       names[I], version, F, I, &vectorElement<F, I>)),
   ...);
}

void _registerDescriptors() {
  REGISTER_DESCRIPTOR(exactmw, calcExactMW);
  REGISTER_DESCRIPTOR(amw, calcAMW);
  REGISTER_DESCRIPTOR(lipinskiHBA, calcLipinskiHBA);
  REGISTER_DESCRIPTOR(lipinskiHBD, calcLipinskiHBD);
  REGISTER_DESCRIPTOR(NumRotatableBonds, calcNumRotatableBonds);
  REGISTER_DESCRIPTOR(NumHBD, calcNumHBD);
  REGISTER_DESCRIPTOR(NumHBA, calcNumHBA);
  REGISTER_DESCRIPTOR(NumHeavyAtoms, calcNumHeavyAtoms);
  REGISTER_DESCRIPTOR(NumAtoms, calcNumAtoms);
  REGISTER_DESCRIPTOR(NumHeteroatoms, calcNumHeteroatoms);
  REGISTER_DESCRIPTOR(NumAmideBonds, calcNumAmideBonds);
  REGISTER_DESCRIPTOR(FractionCSP3, calcFractionCSP3);
  REGISTER_DESCRIPTOR(NumRings, calcNumRings);
  REGISTER_DESCRIPTOR(NumAromaticRings, calcNumAromaticRings);
  REGISTER_DESCRIPTOR(NumAliphaticRings, calcNumAliphaticRings);
  REGISTER_DESCRIPTOR(NumSaturatedRings, calcNumSaturatedRings);
  REGISTER_DESCRIPTOR(NumHeterocycles, calcNumHeterocycles);
  REGISTER_DESCRIPTOR(NumAromaticHeterocycles, calcNumAromaticHeterocycles);
  REGISTER_DESCRIPTOR(NumSaturatedHeterocycles, calcNumSaturatedHeterocycles);
  REGISTER_DESCRIPTOR(NumAliphaticHeterocycles, calcNumAliphaticHeterocycles);
  REGISTER_DESCRIPTOR(NumSpiroAtoms, calcNumSpiroAtoms);
  REGISTER_DESCRIPTOR(NumBridgeheadAtoms, calcNumBridgeheadAtoms);
  REGISTER_DESCRIPTOR(NumAtomStereoCenters, numAtomStereoCenters);
  REGISTER_DESCRIPTOR(NumUnspecifiedAtomStereoCenters,
                      numUnspecifiedAtomStereoCenters);
  REGISTER_DESCRIPTOR(labuteASA, calcLabuteASA);
  REGISTER_DESCRIPTOR(tpsa, calcTPSA);
  REGISTER_DESCRIPTOR(CrippenClogP, calcClogP);
  REGISTER_DESCRIPTOR(CrippenMR, calcMR);
  REGISTER_DESCRIPTOR(chi0v, calcChi0v);
  REGISTER_DESCRIPTOR(chi1v, calcChi1v);
  REGISTER_DESCRIPTOR(chi2v, calcChi2v);
  REGISTER_DESCRIPTOR(chi3v, calcChi3v);
  REGISTER_DESCRIPTOR(chi4v, calcChi4v);
  REGISTER_DESCRIPTOR(chi0n, calcChi0n);
  REGISTER_DESCRIPTOR(chi1n, calcChi1n);
  REGISTER_DESCRIPTOR(chi2n, calcChi2n);
  REGISTER_DESCRIPTOR(chi3n, calcChi3n);
  REGISTER_DESCRIPTOR(chi4n, calcChi4n);
  REGISTER_DESCRIPTOR(hallKierAlpha, calcHallKierAlpha);
  REGISTER_DESCRIPTOR(kappa1, calcKappa1);
  REGISTER_DESCRIPTOR(kappa2, calcKappa2);
  REGISTER_DESCRIPTOR(kappa3, calcKappa3);
  REGISTER_DESCRIPTOR(Phi, calcPhi);

  // the rest of the descriptors in rdkit.Chem.Descriptors._descList, under
  // the names used there
  registerFunc("MaxAbsEStateIndex", MaxAbsEStateIndexVersion,
               calcMaxAbsEStateIndex);
  registerFunc("MaxEStateIndex", MaxEStateIndexVersion, calcMaxEStateIndex);
  registerFunc("MinAbsEStateIndex", MinAbsEStateIndexVersion,
               calcMinAbsEStateIndex);
  registerFunc("MinEStateIndex", MinEStateIndexVersion, calcMinEStateIndex);
  registerFunc("qed", QEDVersion, [](const ROMol &m) { return calcQED(m); });
  registerFunc("SPS", SPSVersion, [](const ROMol &m) { return calcSPS(m); });
  registerFunc("HeavyAtomMolWt", HeavyAtomMolWtVersion, calcHeavyAtomMolWt);
  registerFunc("NumValenceElectrons", NumValenceElectronsVersion,
               [](const ROMol &m) {
                 return static_cast<double>(calcNumValenceElectrons(m));
               });
  registerFunc("NumRadicalElectrons", NumRadicalElectronsVersion,
               [](const ROMol &m) {
                 return static_cast<double>(calcNumRadicalElectrons(m));
               });
  registerFunc("MaxPartialCharge", MaxPartialChargeVersion,
               calcMaxPartialCharge);
  registerFunc("MinPartialCharge", MinPartialChargeVersion,
               calcMinPartialCharge);
  registerFunc("MaxAbsPartialCharge", MaxAbsPartialChargeVersion,
               calcMaxAbsPartialCharge);
  registerFunc("MinAbsPartialCharge", MinAbsPartialChargeVersion,
               calcMinAbsPartialCharge);
  registerFunc("FpDensityMorgan1", FpDensityMorganVersion,
               [](const ROMol &m) { return calcFpDensityMorgan(m, 1); });
  registerFunc("FpDensityMorgan2", FpDensityMorganVersion,
               [](const ROMol &m) { return calcFpDensityMorgan(m, 2); });
  registerFunc("FpDensityMorgan3", FpDensityMorganVersion,
               [](const ROMol &m) { return calcFpDensityMorgan(m, 3); });
#ifdef RDK_HAS_EIGEN3
  registerVector<bcut2D>(
      {"BCUT2D_MWHI", "BCUT2D_MWLOW", "BCUT2D_CHGHI", "BCUT2D_CHGLO",
       "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW", "BCUT2D_MRHI", "BCUT2D_MRLOW"},
      BCUT2DVersion, std::make_index_sequence<8>());
#endif
#ifdef RDK_HAS_EIGEN3
  registerFunc("AvgIpc", avgIpcVersion,
               [](const ROMol &m) { return calcAvgIpc(m); });
#endif
  registerFunc("BalabanJ", balabanJVersion,
               [](const ROMol &m) { return calcBalabanJ(m); });
  registerFunc("BertzCT", bertzCTVersion,
               [](const ROMol &m) { return calcBertzCT(m); });
  registerFunc("Chi0", chi0Version, calcChi0);
  registerFunc("Chi1", chi1Version, calcChi1);
#ifdef RDK_HAS_EIGEN3
  registerFunc("Ipc", ipcVersion, [](const ROMol &m) { return calcIpc(m); });
#endif
  registerVector<peoeVSA>(numberedNames("PEOE_VSA", 14), PEOE_VSAVersion,
                          std::make_index_sequence<14>());
  registerVector<smrVSA>(numberedNames("SMR_VSA", 10), SMR_VSAVersion,
                         std::make_index_sequence<10>());
  registerVector<slogpVSA>(numberedNames("SlogP_VSA", 12), SlogP_VSAVersion,
                           std::make_index_sequence<12>());
  registerVector<estateVSA>(numberedNames("EState_VSA", 11), EState_VSAVersion,
                            std::make_index_sequence<11>());
  registerVector<vsaEState>(numberedNames("VSA_EState", 10), VSA_EStateVersion,
                            std::make_index_sequence<10>());
  REGISTER_DESCRIPTOR(NumAliphaticCarbocycles, calcNumAliphaticCarbocycles);
  REGISTER_DESCRIPTOR(NumAromaticCarbocycles, calcNumAromaticCarbocycles);
  REGISTER_DESCRIPTOR(NumSaturatedCarbocycles, calcNumSaturatedCarbocycles);
  CHECK_INVARIANT(getFragmentDescriptorNames().size() == numFragments,
                  "unexpected number of fragment descriptors");
  registerVector<fragments>(getFragmentDescriptorNames(),
                            FragmentDescriptorsVersion,
                            std::make_index_sequence<numFragments>());
}
}  // namespace

void registerDescriptors() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static std::once_flag once;
  std::call_once(once, _registerDescriptors);
#else
  static bool initialized = false;
  if (!initialized) {
    _registerDescriptors();
    initialized = true;
  }
#endif
}

std::vector<boost::shared_ptr<PropertyFunctor>> Properties::registry;
int Properties::registerProperty(boost::shared_ptr<PropertyFunctor> prop) {
  for (size_t i = 0; i < Properties::registry.size(); ++i) {
    if (registry[i]->getName() == prop->getName()) {
      Properties::registry[i] = std::move(prop);
      return i;
    }
  }
  // XXX Add mutex?
  Properties::registry.push_back(std::move(prop));
  return Properties::registry.size() - 1;
}

int Properties::registerProperty(PropertyFunctor *prop) {
  return Properties::registerProperty(boost::shared_ptr<PropertyFunctor>(prop));
}

std::vector<std::string> Properties::getAvailableProperties() {
  registerDescriptors();
  std::vector<std::string> names;
  for (auto prop : Properties::registry) {
    names.push_back(prop->getName());
  }
  return names;
}

boost::shared_ptr<PropertyFunctor> Properties::getProperty(
    const std::string &name) {
  registerDescriptors();
  for (auto prop : Properties::registry) {
    if (prop.get() && prop->getName() == name) {
      return prop;
    }
  }
  throw KeyErrorException(name);
}

Properties::Properties() : m_properties() {
  registerDescriptors();
  for (auto prop : Properties::registry) {
    m_properties.push_back(prop);
  }
}

Properties::Properties(const std::vector<std::string> &propNames) {
  registerDescriptors();
  for (const auto &name : propNames) {
    m_properties.push_back(Properties::getProperty(name));
  }
}

std::vector<std::string> Properties::getPropertyNames() const {
  std::vector<std::string> names;
  for (auto prop : m_properties) {
    names.push_back(prop->getName());
  }
  return names;
}

namespace {
// what a property gives when its calculation throws, per its failureMode()
double failedValue(const PropertyFunctor &prop) {
  switch (prop.failureMode()) {
    case PropertyFailureMode::THROW:
      throw;
    case PropertyFailureMode::DEFAULT_ZERO:
      return 0.0;
    case PropertyFailureMode::DEFAULT_NAN:
      return std::numeric_limits<double>::quiet_NaN();
  }
}

// computes properties, sharing each descriptor vector between the properties
// that are elements of it
class PropertyCalculator {
 public:
  PropertyCalculator(const ROMol &mol)
    : d_mol(mol) {}

  double operator()(const PropertyFunctor &prop) {
    auto vecProp = dynamic_cast<const VectorElementPropertyFunctor *>(&prop);
    if (!vecProp) {
      try {
        return prop(d_mol);
      } catch (const std::bad_alloc &) {
        throw;
      } catch (const std::exception &) {
        return failedValue(prop);
      }
    }
    auto it = d_vectors.find(vecProp->d_vectorFunc);
    if (it == d_vectors.end()) {
      VectorResult vr;
      try {
        vr.vals = vecProp->d_vectorFunc(d_mol);
      } catch (const std::bad_alloc &) {
        throw;
      } catch (const std::exception &) {
        vr.error = std::current_exception();
      }
      it = d_vectors.emplace(vecProp->d_vectorFunc, std::move(vr)).first;
    }
    const auto &vr = it->second;
    if (vr.error) {
      try {
        std::rethrow_exception(vr.error);
      } catch (const std::exception &) {
        return failedValue(prop);
      }
    }
    if (vecProp->d_index >= vr.vals.size()) {
      try {
        throw ValueErrorException("vector property index out of range");
      } catch (const std::exception &) {
        return failedValue(prop);
      }
    }
    return vr.vals[vecProp->d_index];
  }

 private:
  struct VectorResult {
    std::vector<double> vals;
    std::exception_ptr error;
  };
  const ROMol &d_mol;
  std::map<VectorElementPropertyFunctor::VectorFunc, VectorResult> d_vectors;
};
}  // namespace

std::vector<double> Properties::computeValues(const RDKit::ROMol &mol) const {

  PropertyCalculator calc(mol);
  std::vector<double> res;
  res.reserve(m_properties.size());
  for (const auto &prop : m_properties) {
    res.push_back(calc(*prop));
  }
  return res;
}

std::vector<double> Properties::computeProperties(const RDKit::ROMol &mol,
                                                  bool annotate) const {
  auto res = computeValues(mol);
  if (annotate) {
    for (size_t i = 0; i < m_properties.size(); ++i) {
      mol.setProp<double>(m_properties[i]->getName(), res[i]);
    }
  }
  return res;
}

std::vector<std::vector<double>> Properties::computePropertiesForMols(
    const std::vector<const RDKit::ROMol *> &mols, int numThreads) const {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::vector<double>> res(
      mols.size(), std::vector<double>(m_properties.size(), nan));

  // Each distinct molecule is computed exactly once: property calculators
  // may cache results on the molecule (or initialize its ring info), so the
  // same molecule must never be processed by more than one thread.
  std::vector<size_t> todo;
  std::vector<std::pair<size_t, size_t>> repeats;  // (repeat, first)
  {
    std::unordered_map<const ROMol *, size_t> firstSeen;
    firstSeen.reserve(mols.size());
    for (size_t i = 0; i < mols.size(); ++i) {
      if (!mols[i]) {
        continue;
      }
      auto [it, inserted] = firstSeen.emplace(mols[i], i);
      if (inserted) {
        todo.push_back(i);
      } else {
        repeats.emplace_back(i, it->second);
      }
    }
  }

  // Workers pull the next unprocessed entry of `todo` from a shared atomic
  // index, so uneven molecule costs balance across threads.
  std::atomic<size_t> nextIdx{0};
  auto worker = [&]() {
    for (size_t t = nextIdx.fetch_add(1); t < todo.size();
         t = nextIdx.fetch_add(1)) {
      // virtual, so that subclasses (e.g. normalizers) apply here too
      res[todo[t]] = computeProperties(*mols[todo[t]]);
    }
  };

  auto nThreads = std::min<size_t>(todo.size(), getNumThreadsToUse(numThreads));
#ifdef RDK_BUILD_THREADSAFE_SSS
  if (nThreads <= 1) {
    worker();
  } else {
    // declared after res so that, if an exception unwinds the stack, the
    // futures are destroyed (waiting for their tasks) before res is
    std::vector<std::future<void>> tasks;
    tasks.reserve(nThreads);
    for (size_t ti = 0; ti < nThreads; ++ti) {
      tasks.emplace_back(std::async(std::launch::async, worker));
    }
    for (auto &task : tasks) {
      task.get();
    }
  }
#else
  // threads are compiled out: getNumThreadsToUse() is 1 here anyway
  RDUNUSED_PARAM(nThreads);
  worker();
#endif

  for (const auto &[repeat, first] : repeats) {
    res[repeat] = res[first];
  }
  return res;
}

void Properties::annotateProperties(RDKit::ROMol &mol) const {
  computeProperties(mol, true);
}

PROP_RANGE_QUERY *makePropertyRangeQuery(const std::string &name, double min,
                                         double max) {
  auto *filter = new PROP_RANGE_QUERY(min, max);
  filter->setDataFunc(Properties::getProperty(name)->d_dataFunc);
  return filter;
}
}  // namespace Descriptors
}  // namespace RDKit
