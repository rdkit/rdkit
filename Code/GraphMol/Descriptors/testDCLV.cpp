#include <RDGeneral/test.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/RDKitBase.h>

#include "DCLV.h"

using namespace RDKit;
using namespace std;
using namespace Descriptors;

void testPDB() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing DoubleCubicLatticeVolume from PDB"
      << endl;

  string rdbase = getenv("RDBASE");
  string fName = rdbase += "/Code/GraphMol/Descriptors/test_data/1mup.pdb";

  ROMol* m;
  {
    const bool sanitize = false;
    const bool removeHs = false;
    m = PDBFileToMol(fName, sanitize, removeHs);
  }

  // test default params
  DoubleCubicLatticeVolume vol_default(m);

  TEST_ASSERT(fabs(vol_default.getSurfaceArea() - 8330.59) < 0.05);
  TEST_ASSERT(fabs(vol_default.getVolume() - 31789.6) < 0.05);
  TEST_ASSERT(fabs(vol_default.getVDWVolume() - 15355.3) < 0.05);
  TEST_ASSERT(fabs(vol_default.getCompactness() - 1.7166) < 0.05);
  TEST_ASSERT(fabs(vol_default.getPackingDensity() - 0.48303) < 0.05);

  // test set depth and radius
  DoubleCubicLatticeVolume vol_depthrad(m, true, false, 1.6, 6);

  TEST_ASSERT(fabs(vol_depthrad.getSurfaceArea() - 8186.06) < 0.05);
  TEST_ASSERT(fabs(vol_depthrad.getVolume() - 33464.5) < 0.05);
  TEST_ASSERT(fabs(vol_depthrad.getVDWVolume() - 15350.7) < 0.05);
  TEST_ASSERT(fabs(vol_depthrad.getCompactness() - 1.63005) < 0.05);
  TEST_ASSERT(fabs(vol_depthrad.getPackingDensity() - 0.458717) < 0.05);

  // test include ligand
  DoubleCubicLatticeVolume vol_withlig(m, true, true);

  TEST_ASSERT(fabs(vol_withlig.getSurfaceArea() - 8010.56) < 0.05);
  TEST_ASSERT(fabs(vol_withlig.getVolume() - 31228.4) < 0.05);
  TEST_ASSERT(fabs(vol_withlig.getVDWVolume() - 15155.7) < 0.05);
  TEST_ASSERT(fabs(vol_withlig.getCompactness() - 1.67037) < 0.05);
  TEST_ASSERT(fabs(vol_withlig.getPackingDensity() - 0.48532) < 0.05);
}

void testSDF() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing DoubleCubicLatticeVolume from PDB"
      << endl;

  string rdbase = getenv("RDBASE");
  string fName = rdbase +=
      "/Code/GraphMol/Descriptors/test_data/TZL_model.sdf";  // ligand
                                                             // from
                                                             // 1MU

  SDMolSupplier reader(fName, false, false, false);
  while (!reader.atEnd()) {
    ROMol* mol = reader.next();
    TEST_ASSERT(mol);
    DoubleCubicLatticeVolume vol_sdf(mol, false);
  
  

  // NOTE - expected values generated from Roger's original C code
  // Original did not return surface area for Ligand only
  // so no check for Surface Area or Compactness

  TEST_ASSERT(fabs(vol_sdf.getVolume() - 1048.53) < 0.05);
  TEST_ASSERT(fabs(vol_sdf.getVDWVolume() - 231.971) < 0.05);
  TEST_ASSERT(fabs(vol_sdf.getPackingDensity() - 0.221234) < 0.05);
  }
}

int main() {
  RDLog::InitLogs();
  testPDB();
  testSDF();
}