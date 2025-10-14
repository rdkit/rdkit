#include <catch2/catch_all.hpp>
#include <string>

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("Chirality::findPotentialStereo", "[stereo]") {
  auto cases = {
      "COC1/C=C/OC2(C)Oc3c(C)c(O)c4c(O)c(c(/C=N/OC(c5ccccc5)c5ccccc5)cc4c3C2=O)NC(=O)/C(C)=C\\C=C\\C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C1C",
      "O=C1NC(=O)C2=C1c1cn(c3ccccc13)CCO[C@H](CNCc1ccccc1)CCn1cc2c2ccccc21",
      "CN[C@@H](C)C(=O)N[C@H](C(=O)N1C[C@@H]2C[C@H]1C(=O)N[C@@H](Cc1ccc3ccccc3c1)C(=O)N[C@@H](C(=O)O)Cc1ccc(cc1)OCc1cn2nn1)C1CCCCC1",
      "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
      "C12(CCCCC1)CCCC3CCCCC23",
      "c1ccc2c(c1)c3ccccc3c4ccccc24",
      "F[C@@H](Cl)[C@H](Br)[C@](I)(O)[C@@](N)(C#N)C(=O)O",
      "C[S+](C)(C)[O-]",
      "C(C(C(C(C(C(CO)O)O)O)O)O)O",
  };

  for (auto smiles : cases) {
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmilesToMol(smiles)};
    REQUIRE(mol);

    BENCHMARK("Chirality::findPotentialStereo: " + std::string(smiles)) {
      return Chirality::findPotentialStereo(*mol);
    };
  }
}

TEST_CASE("CIPLabeler::CIPLabeler", "[stereo]") {
  auto cases = {
      "COC1/C=C/OC2(C)Oc3c(C)c(O)c4c(O)c(c(/C=N/OC(c5ccccc5)c5ccccc5)cc4c3C2=O)NC(=O)/C(C)=C\\C=C\\C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C1C",
      "O=C1NC(=O)C2=C1c1cn(c3ccccc13)CCO[C@H](CNCc1ccccc1)CCn1cc2c2ccccc21",
      "CN[C@@H](C)C(=O)N[C@H](C(=O)N1C[C@@H]2C[C@H]1C(=O)N[C@@H](Cc1ccc3ccccc3c1)C(=O)N[C@@H](C(=O)O)Cc1ccc(cc1)OCc1cn2nn1)C1CCCCC1",
      "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
      "C12(CCCCC1)CCCC3CCCCC23",
      "c1ccc2c(c1)c3ccccc3c4ccccc24",
      "F[C@@H](Cl)[C@H](Br)[C@](I)(O)[C@@](N)(C#N)C(=O)O",
      "C[S+](C)(C)[O-]",
      "C(C(C(C(C(C(CO)O)O)O)O)O)O",
  };

  for (auto smiles : cases) {
    std::unique_ptr<RDKit::ROMol> mol{RDKit::SmilesToMol(smiles)};
    REQUIRE(mol);

    BENCHMARK("CIPLabeler::assignCIPLabels: " + std::string(smiles)) {
      return Chirality::findPotentialStereo(*mol);
    };
  }
}
