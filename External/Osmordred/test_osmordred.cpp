#include "Osmordred.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>


using namespace RDKit;

TEST_CASE("Test Osmordred") {
  SECTION("Calling all functions") {
    auto mol = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"_smiles;

    {
      std::vector<double> v = calcABCIndex(*mol); 
    }
    {
      std::vector<int> v = calcAcidBase(*mol); 
    }
    {
      std::vector<int> v = calcAromatic(*mol); 
    }
    {
      std::vector<int> v = calcAtomCounts(*mol);
    }
    {
      std::vector<double> v = calcBalabanJ(*mol); 
    }
    {
      std::vector<double> v = calcBertzCT(*mol);
    }
    {
      std::vector<int> v = calcBondCounts(*mol);
    }
    {
      std::vector<double> v = calcVertexAdjacencyInformation(*mol); 
    }
    {
      std::vector<double> v = calcWeight(*mol); 
    }
    {
      std::vector<int> v = calcWienerIndex(*mol); 
    }
    {
      std::vector<double> v = calcVdwVolumeABC(*mol); 
    }
    {
      std::vector<double> v = calcTopoPSA(*mol); 
    }
    {
      std::vector<double> v = calcSLogP(*mol); 
    }
    {
      std::vector<double> v = calcHydrogenBond(*mol); 
    }
    {
      std::vector<double> v = calcLogS(*mol);
    }
    {
      std::vector<int> v = calcLipinskiGhose(*mol);
    }
    {
      std::vector<double> v = calcMcGowanVolume(*mol);
    }
    {
      std::vector<double> v = calcPolarizability(*mol);
    }
    {
      std::vector<double> v = calcRotatableBond(*mol); 
    }
    {
      std::vector<double> v = calcFragmentComplexity(*mol); 
    }
    {
      std::vector<double> v = calcConstitutional(*mol); 
    }
    {
      std::vector<double> v = calcTopologicalIndex(*mol);
    }
    {
      std::vector<double> v = calcDetourMatrixDescs(*mol);
    }
    {
      std::vector<double> v = calcDetourMatrixDescsL(*mol);
    }
    {
      std::vector<double> v = calcDistMatrixDescs(*mol);
    }
    {
      std::vector<double> v = calcDistMatrixDescsL(*mol);
    }
    {
      std::vector<double> v = calcAdjMatrixDescs(*mol);
    }
    {
      std::vector<double> v = calcAdjMatrixDescsL(*mol);
    }
    {
      std::vector<double> v = calcCarbonTypes(*mol);
    }
    {
      std::vector<double> v = calcEccentricConnectivityIndex(*mol); 
    }
    {
      std::vector<double> v = calcBaryszMatrixDescsL(*mol);
    }
    {
      std::vector<double> v = calcBaryszMatrixDescs(*mol);
    }
    {
      std::vector<double> v = calcZagrebIndex(*mol);
    }
    {
      std::vector<double> v = calcMoeType(*mol); 
    }
    {
      std::vector<double> v = calcMolecularDistanceEdgeDescs(*mol);
    }
    {
      std::vector<double> v = calcEStateDescs(*mol);
      std::vector<double> v2 = calcEStateDescs(*mol, true);
    }
    {
      std::vector<double> v = calcWalkCounts(*mol);
    }
    {
      std::vector<double> v = calcTopologicalChargeDescs(*mol);
    }
    {
      std::vector<double> v = calcAllChiDescriptors(*mol);
    }
    {
      std::vector<double> v = calcPathCount(*mol);
    }
    {
      std::vector<double> v = calcKappaShapeIndex(*mol);
    }
    {
      std::vector<int> v = calcRingCount(*mol); 
    }
    {
      std::vector<double> v = calcMolecularId(*mol);
    }
    {
      std::vector<double> v = calcBCUTs(*mol); // 10x faster the
    }
    {
      std::vector<double> v = calcAutoCorrelation(*mol);
    }
    {
      std::vector<double> v = calcFramework(*mol); 
    }
    {
      std::vector<double> v = calcExtendedTopochemicalAtom(*mol);
    }
    {
      std::vector<double> v = calculateETADescriptors(*mol);
    }
    {
      std::vector<double> v = calcChipath(*mol);
    }
    {
      std::vector<double> v = calcChichain(*mol);
    }
    {
      std::vector<double> v = calcChicluster(*mol);
    }
    {
      std::vector<double> v = calcChipathcluster(*mol);
    }
    {
      int v = calcAcidicGroupCount(*mol); (void) v;
    }
    {
      int v = calcBasicGroupCount(*mol); (void) v;
    }
    {
      int v = countAromaticAtoms(*mol); (void) v;
    }
    {
      int v = countAromaticBonds(*mol); (void) v;
    }
    {
      std::vector<double> v = calcBEStateDescs(*mol);
    }
    {
      std::vector<double> v = calcHEStateDescs(*mol);
    }
    {
      std::vector<double> v = calcAlphaKappaShapeIndex(*mol); // closer "missing" k3 path count not correct on few cases
    }
    {
      std::vector<double> v = calcAbrahams(*mol); // Platts, Butina, Abraham, Hersey  paper J Chem Inf Comput Sci. 1999 30/8/01;39(5):835-45
    }
    {
      std::vector<double> v = calcPol(*mol);
    }
    {
      std::vector<double> v = calcMR(*mol);
    }
    {
      std::vector<double> v = calcFlexibility(*mol);
    }
    {
      std::vector<double> v = calcODT(*mol);
    }
    {
      std::vector<double> v = calcSchultz(*mol);
    }
    {
      std::vector<double> v = calcRNCG_RPCG(*mol);
    }
    {
      std::vector<double> v = calcAZV(*mol);
    }
    {
      std::vector<double> v = calcASV(*mol);
    }
    {
      std::vector<double> v = calcDSV(*mol);
    }
    {
      std::vector<double> v = calcAZS(*mol);
    }
    {
      std::vector<double> v = calcASZ(*mol);
    }
    {
      std::vector<double> v = calcDN2S(*mol);
    }
    {
      std::vector<double> v = calcDN2I(*mol);
    }
    {
      std::vector<double> v = calcASI(*mol);
    }
    {
      std::vector<double> v = calcDSI(*mol);
    }
    {
      std::vector<double> v = calcASN(*mol);
    }
    {
      std::vector<double> v = calcDSN(*mol);
    }
    {
      std::vector<double> v = calcDN2N(*mol);
    }
    {
      std::vector<double> v = calcANS(*mol);
    }
    {
      std::vector<double> v = calcANV(*mol);
    }
    {
      std::vector<double> v = calcAZN(*mol);
    }
    {
      std::vector<double> v = calcANZ(*mol);
    }
    {
      std::vector<double> v = calcANI(*mol);
    }
    {
      std::vector<double> v = calcDSZ(*mol);
    }
    {
      std::vector<double> v = calcANN(*mol);
    }
    {
      std::vector<double> v = calcDN2Z(*mol);
    }
    {
      std::vector<double> v = calcANMat(*mol);
    }
    {
      std::vector<double> v = calcAZMat(*mol);
    }
    {
      std::vector<double> v = calcASMat(*mol);
    }
    {
      std::vector<double> v = calcDSMat(*mol);
    }
    {
      std::vector<double> v = calcDN2Mat(*mol);
    }
    {
      std::vector<double> v = calcFrags(*mol);
    }
    {
      std::vector<double> v = calcAddFeatures(*mol);
    }
    {
      std::vector<double> v = calcInformationContent(*mol);
    }
  
  }

}
