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

RDKIT_DESCRIPTORS_EXPORT  bool hasOsmordredSupport();
  
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcABCIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAcidBase(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAromatic(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAtomCounts(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBalabanJ(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBertzCT(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcBondCounts(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVertexAdjacencyInformation(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWeight(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcWienerIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVdwVolumeABC(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopoPSA(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSLogP(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHydrogenBond(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcLogS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcLipinskiGhose(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMcGowanVolume(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPolarizability(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRotatableBond(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFragmentComplexity(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcConstitutional(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalIndex(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcCarbonTypes(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEccentricConnectivityIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcZagrebIndex(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMoeType(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularDistanceEdgeDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEStateDescs(const RDKit::ROMol& mol, bool extended = false);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWalkCounts(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalChargeDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAllChiDescriptors(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPathCount(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcKappaShapeIndex(const RDKit::ROMol& mol); // closer "missing" k3 path count not correct on few cases
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcRingCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularId(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBCUTs(const RDKit::ROMol& mol); // 10x faster the
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAutoCorrelation(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFramework(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcExtendedTopochemicalAtom(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calculateETADescriptors(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipath(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChichain(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChicluster(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipathcluster(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int calcAcidicGroupCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int calcBasicGroupCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticAtoms(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticBonds(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBEStateDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHEStateDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAlphaKappaShapeIndex(const RDKit::ROMol& mol); // closer "missing" k3 path count not correct on few cases
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahams(const RDKit::ROMol& mol); // Platts, Butina, Abraham, Hersey  paper J Chem Inf Comput Sci. 1999 30/8/01;39(5):835-45
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPol(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMR(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFlexibility(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcODT(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSchultz(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2S(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2I(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2N(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Z(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Mat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFrags(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAddFeatures(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcInformationContent(const RDKit::ROMol& mol, int maxradius=5); // Inspired by 1984 Basak paper

// Aggregated fast path that calls all Osmordred descriptors in C++
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordred(const RDKit::ROMol& mol);
} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit

#endif //_DESCRIPTORS_H
