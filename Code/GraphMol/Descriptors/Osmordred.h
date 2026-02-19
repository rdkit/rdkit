//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC’s and others
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
#include <RDGeneral/export.h>

#ifndef _OSMORDRED_H
#define _OSMORDRED_H

#include <GraphMol/ROMol.h>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// v2.0: Control function to check if Gasteiger parameters exist for all atoms
// Returns true if all atoms have parameters for their specific environment,
// false otherwise Use this BEFORE calling any function that uses Gasteiger
// charges to avoid crashes
RDKIT_DESCRIPTORS_EXPORT bool checkGasteigerParameters(const ROMol &mol);

// v2.0: Filter function to check if a molecule is too large (will cause hangs)
// Returns true if molecule has >10 rings OR >200 heavy atoms
// Use this to filter out overly complex molecules before descriptor calculation
RDKIT_DESCRIPTORS_EXPORT bool isMoleculeTooLarge(const ROMol &mol);

// Group 1 + 2: Basic physchem + counts/rules
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcABCIndex(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAcidBase(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAromatic(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAtomCounts(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcBondCounts(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWeight(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVdwVolumeABC(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopoPSA(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSLogP(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHydrogenBond(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcLogS(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcLipinskiGhose(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMcGowanVolume(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPolarizability(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRotatableBond(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFragmentComplexity(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcConstitutional(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAddFeatures(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT int calcAcidicGroupCount(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT int calcBasicGroupCount(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticAtoms(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticBonds(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahams(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPol(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMR(const ROMol &mol);

// Group 3: Topological indices/connectivity/shape
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBalabanJ(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBertzCT(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVertexAdjacencyInformation(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcWienerIndex(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalIndex(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcCarbonTypes(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEccentricConnectivityIndex(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcZagrebIndex(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularDistanceEdgeDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWalkCounts(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalChargeDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAllChiDescriptors(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPathCount(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcKappaShapeIndex(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAlphaKappaShapeIndex(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFlexibility(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcODT(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSchultz(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRNCG_RPCG(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcRingDescriptors(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularId(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFramework(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcExtendedTopochemicalAtom(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calculateETADescriptors(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipath(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChichain(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChicluster(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipathcluster(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcInformationContent(
    const ROMol &mol, int maxradius = 5);  // Inspired by 1984 Basak paper

// Group 4: Matrix/autocorr/EState/fragments
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescsL(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescsL(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescsL(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescsL(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescs(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMoeType(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEStateDescs(
    const ROMol &mol, bool extended = false);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBEStateDescs(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHEStateDescs(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBCUTs(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAutoCorrelation(
    const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANMat(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZMat(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASMat(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSMat(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Mat(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZV(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASV(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSV(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZS(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASZ(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2S(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2I(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASI(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSI(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASN(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSN(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2N(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANS(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANV(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZN(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANZ(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANI(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSZ(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANN(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Z(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFrags(const ROMol &mol);

// Aggregated fast path that calls all Osmordred descriptors in C++
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordred(const ROMol &mol);

// v2.0: Single molecule with timeout protection (default 60 seconds)
// Returns NaN vector (3585 NaN values) if computation exceeds timeout
// This is the RECOMMENDED function for production use to prevent hanging
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordredWithTimeout(
    const ROMol &mol, int timeout_seconds = 60);

// v2.0: Batch from SMILES: parses each SMILES with SmilesToMol() -> NEW mol.
// Tautomer canonical LOST.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatch(
    const std::vector<std::string> &smiles_list, int n_jobs = 0);

// v2.0: Batch from mol objects (Python Mol via ToBinary/MolPickler). PRESERVES
// tautomer canonical.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>>
calcOsmordredBatchFromMols(const std::vector<const ROMol *> &mols,
                           int n_jobs = 0);

// v2.0: Get descriptor names in the same order as calcOsmordred returns values
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getOsmordredDescriptorNames();
  
}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

#endif  //_DESCRIPTORS_H
