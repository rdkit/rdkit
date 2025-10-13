#include <catch2/catch_all.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("smiles") {
  const std::vector<std::string> testSmiles = {
      "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
      "C12(CCCCC1)CCCC3CCCCC23",
      "c1ccc2c(c1)c3ccccc3c4ccccc24",
      "F[C@@H](Cl)[C@H](Br)[C@](I)(O)[C@@](N)(C#N)C(=O)O",
      "C[S+](C)(C)[O-]",
      "C(C(C(C(C(C(CO)O)O)O)O)O)O"};

  for (const auto &smiles : testSmiles) {
    std::unique_ptr<ROMol> mol{SmilesToMol(smiles)};
    REQUIRE(mol);

    BENCHMARK("SmilesToMol: " + smiles) {
      std::unique_ptr<ROMol> temp{SmilesToMol(smiles)};
      return temp;
    };

    BENCHMARK("MolToSmiles: " + smiles) { return MolToSmiles(*mol); };
  }
}
