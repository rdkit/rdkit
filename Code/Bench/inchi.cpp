#include <catch2/catch_all.hpp>
#include <string>
#include <utility>
#include <vector>

#include "bench_common.hpp"

#include <GraphMol/Conformer.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <INCHI-API/inchi.h>

using namespace RDKit;

namespace {
std::vector<ROMol> loadConformerSamples(
    const std::vector<std::string> &smiles, bool is3D) {
  std::vector<ROMol> ret;
  ret.reserve(smiles.size());
  for (const auto &smi : smiles) {
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(mol);
    RWMol withConformer(*mol);
    auto *conf = new Conformer(withConformer.getNumAtoms());
    conf->set3D(is3D);
    for (auto atomIdx = 0u; atomIdx < withConformer.getNumAtoms();
         ++atomIdx) {
      const auto idx = static_cast<double>(atomIdx);
      conf->setAtomPos(
          atomIdx,
          RDGeom::Point3D(idx, static_cast<double>(atomIdx % 3) * 0.5,
                          is3D ? idx * 0.1 : 0.0));
    }
    withConformer.addConformer(conf, true);
    ret.push_back(std::move(withConformer));
  }
  return ret;
}

std::vector<std::string> molsToInchis(const std::vector<ROMol> &mols) {
  std::vector<std::string> inchis;
  inchis.reserve(mols.size());
  for (const auto &mol : mols) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  return inchis;
}
}  // namespace

TEST_CASE("MolToInchi", "[inchi]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToInchi") {
    return molsToInchis(samples);
  };
}

TEST_CASE("MolToInchi conformer paths", "[inchi]") {
  auto noCandidate2D =
      loadConformerSamples({"CCCC", "C1CCCCC1", "CCOCC"}, false);
  auto asymmetric2D =
      loadConformerSamples({"CC=CC", "CC=C(C)C", "CCC=CC"}, false);

  std::vector<std::string> symmetricSmiles = {
      "COC(=O)C=C1CCC(C)CC1"};
  for (auto ringSize : {6u, 12u, 24u, 48u}) {
    symmetricSmiles.emplace_back("COC(=O)C=C1" +
                                 std::string(ringSize - 1, 'C') + "1");
  }
  auto symmetric2D = loadConformerSamples(symmetricSmiles, false);
  auto symmetric3D = loadConformerSamples(symmetricSmiles, true);

  BENCHMARK("MolToInchi/2D-no-candidate") {
    return molsToInchis(noCandidate2D);
  };
  BENCHMARK("MolToInchi/2D-asymmetric-candidate") {
    return molsToInchis(asymmetric2D);
  };
  BENCHMARK("MolToInchi/2D-symmetric-candidate") {
    return molsToInchis(symmetric2D);
  };
  BENCHMARK("MolToInchi/3D-symmetric-candidate") {
    return molsToInchis(symmetric3D);
  };
}

TEST_CASE("InchiToInchiKey", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToInchiKey") {
    std::vector<std::string> inchikeys;
    for (auto &inchi : inchis) {
      inchikeys.push_back(InchiToInchiKey(inchi));
    }
    return inchikeys;
  };
}

TEST_CASE("InchiToMol", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToMol") {
    std::vector<std::unique_ptr<ROMol>> mols;
    for (auto &inchi : inchis) {
      ExtraInchiReturnValues rv_inner;
      mols.emplace_back(InchiToMol(inchi, rv_inner));
    }
    return mols;
  };
}
