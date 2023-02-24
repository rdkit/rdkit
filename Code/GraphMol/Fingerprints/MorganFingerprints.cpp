//
//  Copyright (c) 2009-2022, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
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
//  Created by Greg Landrum, July 2008
//
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

namespace RDKit {
namespace MorganFingerprints {

SparseIntVect<uint32_t> *getFingerprint(
    const ROMol &mol, unsigned int radius, std::vector<uint32_t> *invariants,
    const std::vector<uint32_t> *fromAtoms, bool useChirality,
    bool useBondTypes, bool useCounts, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  bool countSimulation = false;
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
      MorganFingerprint::getMorganGenerator<std::uint32_t>(
          radius, countSimulation, useChirality, useBondTypes,
          onlyNonzeroInvariants, includeRedundantEnvironments));
  RDKit::FingerprintFuncArguments args;
  args.fromAtoms = fromAtoms;
  args.customAtomInvariants = invariants;
  AdditionalOutput ao;
  if (atomsSettingBits) {
    args.additionalOutput = &ao;
    ao.allocateBitInfoMap();
  }

  SparseIntVect<uint32_t> *res;
  if (!useCounts) {
    auto tmp = fpgen->getSparseFingerprint(mol, args);
    res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());
    for (auto idx : *(tmp->dp_bits)) {
      res->setVal(idx, 1);
    }
  } else {
    res = fpgen->getSparseCountFingerprint(mol, args).release();
  }

  if (atomsSettingBits) {
    atomsSettingBits->clear();
    for (const auto &pr : *(ao.bitInfoMap)) {
      (*atomsSettingBits)[pr.first] = pr.second;
    }
  }

  return res;
}

SparseIntVect<uint32_t> *getHashedFingerprint(
    const ROMol &mol, unsigned int radius, unsigned int nBits,
    std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  if (nBits == 0) {
    throw ValueErrorException("nBits can not be zero");
  }

  bool countSimulation = false;
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
      MorganFingerprint::getMorganGenerator<std::uint32_t>(
          radius, countSimulation, useChirality, useBondTypes,
          onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
          nBits));
  RDKit::FingerprintFuncArguments args;
  args.fromAtoms = fromAtoms;
  args.customAtomInvariants = invariants;
  AdditionalOutput ao;
  if (atomsSettingBits) {
    args.additionalOutput = &ao;
    ao.allocateBitInfoMap();
  }

  auto res = fpgen->getCountFingerprint(mol, args).release();

  if (atomsSettingBits) {
    atomsSettingBits->clear();
    for (const auto &pr : *(ao.bitInfoMap)) {
      (*atomsSettingBits)[pr.first] = pr.second;
    }
  }

  return res;
}

ExplicitBitVect *getFingerprintAsBitVect(
    const ROMol &mol, unsigned int radius, unsigned int nBits,
    std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
  if (nBits == 0) {
    throw ValueErrorException("nBits can not be zero");
  }

  bool countSimulation = false;
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
      MorganFingerprint::getMorganGenerator<std::uint32_t>(
          radius, countSimulation, useChirality, useBondTypes,
          onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
          nBits));
  RDKit::FingerprintFuncArguments args;
  args.fromAtoms = fromAtoms;
  args.customAtomInvariants = invariants;
  AdditionalOutput ao;
  if (atomsSettingBits) {
    args.additionalOutput = &ao;
    ao.allocateBitInfoMap();
  }

  auto res = fpgen->getFingerprint(mol, args).release();

  if (atomsSettingBits) {
    atomsSettingBits->clear();
    for (const auto &pr : *(ao.bitInfoMap)) {
      (*atomsSettingBits)[pr.first] = pr.second;
    }
  }

  return res;
}

}  // end of namespace MorganFingerprints
}  // end of namespace RDKit
