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

//! Which equivalence key the InformationContent descriptors use.
/*!
  Measured against POLLY -- the program Basak's own group wrote and used --
  on 411 molecules x r=0..5 (2466 values), tolerance 0.0015:

    BASAK     88.0%      (this is what osmordred v2 did)
    EXTENDED  66.1%      (osmordred v3; added the neighbour degree)
    mordred   49.8%      (for reference; not reachable from here)

  So v2's key was substantially closer to Basak than the mordred package is,
  and v3 moved away from Basak rather than towards it. The difference is the
  neighbour DEGREE: Basak labels a vertex by (element, valency), and in a
  hydrogen-filled graph an sp2 carbon has degree 3 but valency 4, so keying on
  degree splits atoms Basak keeps together. On Basak's own worked example
  (2-butenol, Roy/Basak/Harriss/Magnuson 1983, Table 1) BASAK reproduces
  IC1 = 2.0349 with the partition [1,1,1,1,2,7]; EXTENDED gives 2.4997 and
  [1,1,1,1,2,2,5], which is mordred's answer.
*/
enum class ICKeyFlavor {
  BASAK,     //!< neighbour degree excluded; reproduces Basak/POLLY (default)
  EXTENDED,  //!< neighbour degree included (osmordred v3 behaviour)
};

//! How aromatic bonds are encoded in the equivalence key.
/*!
  DISTINCT scores 88.0% against POLLY, KEKULIZED 82.9%, so DISTINCT is the
  default. KEKULIZED is exposed because POLLY predates aromatic bond types and
  it is a reasonable thing to want to test, not because it scores better.
*/
enum class ICAromaticHandling {
  DISTINCT,   //!< aromatic bonds get their own bond code (default)
  KEKULIZED,  //!< kekulize first and use integer bond orders
};

//! Which per-atom label the equivalence key is built from.
/*!
  Basak's paper defines a vertex as (element, valency) (p.746), so VALENCY is
  the literal reading of the text. Measured against POLLY with this key the two
  are near-equivalent: DEGREE 88.0%, VALENCY 87.8% of 2466 values. DEGREE is
  the default only because it is marginally ahead and is the existing
  behaviour, not because valency is wrong.

  (An earlier revision of this comment quoted 83.2% for VALENCY. That figure
  came from a path-code reconstruction of the paper's condition (iii), which is
  a different algorithm from this key, and did not belong here.)
*/
enum class ICVertexLabel {
  DEGREE,   //!< graph degree in the hydrogen-filled graph (default)
  VALENCY,  //!< Basak's (element, valency) as written in the paper
};

//! Options for the InformationContent family (IC/TIC/SIC/BIC/CIC/MIC/ZMIC).
/*!
  Default-constructed, this reproduces Basak/POLLY as closely as osmordred
  currently can (88.0% of 2466 reference values).

  The residual disagreement is concentrated where RDKit's graph breaks a
  symmetry the molecule actually has:

    plain molecules                                      92.3%
    tautomer-ambiguous                                   88.2%
    resonance-asymmetric (nitro, carboxylate, sulfonate) 33.3%

  RDKit writes a nitro group as [N+](=O)[O-], so its two chemically equivalent
  oxygens take different canonical ranks and get different keys; every nitro
  molecule in the reference set disagrees with POLLY at one or more orders,
  always by over-splitting. equalizeDelocalizedBonds is a partial mitigation --
  it gives one bond code to every bond inside a nitro/carboxylate/sulfonate
  group, which lifts those molecules from 33.3% to 40.3%. It is OFF by default
  because it does not close the gap and it changes the descriptor's meaning;
  note also that every nitro molecule in the reference set is also aromatic, so
  the two effects cannot be separated on the data available.
*/
struct RDKIT_DESCRIPTORS_EXPORT InformationContentOptions {
  ICKeyFlavor keyFlavor = ICKeyFlavor::BASAK;
  ICAromaticHandling aromaticHandling = ICAromaticHandling::DISTINCT;
  ICVertexLabel vertexLabel = ICVertexLabel::DEGREE;
  //! give one bond code to every bond inside a delocalized group; see above
  bool equalizeDelocalizedBonds = false;
};


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
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipath(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChichain(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChicluster(const ROMol &mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipathcluster(
    const ROMol &mol);
//! InformationContent (IC/TIC/SIC/BIC/CIC/MIC/ZMIC), Basak's neighbourhood
//! complexity indices.
/*!
  \param mol        the molecule of interest
  \param maxradius  highest neighbourhood order; must be >= 0
  \param options    see InformationContentOptions; the default reproduces
                    Basak/POLLY

  \return 7 * (maxradius + 1) values, in blocks: IC, TIC, SIC, BIC, CIC, MIC,
          ZMIC, each running r = 0..maxradius. Hydrogens are added internally,
          so pass the molecule without explicit Hs.
*/
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcInformationContent(
    const ROMol &mol, int maxradius,
    const InformationContentOptions &options);

//! \overload  uses default (Basak) options
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcInformationContent(
    const ROMol &mol, int maxradius = 5);

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

// v2.0: Timeout constant for Osmordred computation (60 seconds = 1 minute)
constexpr int OSMORDRED_TIMEOUT_SECONDS = 60;

// v2.0: Single molecule with timeout protection (default 60 seconds)
// Returns NaN vector (3585 NaN values) if computation exceeds timeout
// This is the RECOMMENDED function for production use to prevent hanging
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordredWithTimeout(
    const ROMol &mol, int timeout_seconds = OSMORDRED_TIMEOUT_SECONDS);

// v2.0: Batch from SMILES: parses each SMILES with SmilesToMol() -> NEW mol.
// Tautomer canonical LOST.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatch(
    const std::vector<std::string> &smiles_list, int n_jobs = 0,
    int timeout_seconds = OSMORDRED_TIMEOUT_SECONDS);

// v2.0: Batch from mol objects (Python Mol via ToBinary/MolPickler). PRESERVES
// tautomer canonical.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>>
calcOsmordredBatchFromMols(const std::vector<const ROMol *> &mols,
                           int n_jobs = 0,
			   int timeout_seconds = OSMORDRED_TIMEOUT_SECONDS);
// v2.0: Get descriptor names in the same order as calcOsmordred returns values
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getOsmordredDescriptorNames();

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

#endif  //_DESCRIPTORS_H
