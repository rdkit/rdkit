#include <catch2/catch_all.hpp>

#include <set>
#include <memory>
#include <string>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

namespace {

std::vector<ROMOL_SPTR> make_reagent_pool(unsigned int count,
                                         const std::string &suffix) {
  std::vector<ROMOL_SPTR> pool;
  pool.reserve(count);
  for (unsigned int i = 0; i < count; ++i) {
    const std::string smiles = std::string(i + 1, 'C') + suffix;
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    pool.push_back(ROMOL_SPTR(mol.release()));
  }
  return pool;
}

std::vector<ROMOL_SPTR> make_acid_pool(unsigned int count) {
  return make_reagent_pool(count, "(=O)O");
}

std::vector<ROMOL_SPTR> make_amine_pool(unsigned int count) {
  return make_reagent_pool(count, "N");
}

constexpr unsigned int BENCH_REPEATS = 8;

std::size_t run_uncached(const ChemicalReaction &rxn,
                         const std::vector<ROMOL_SPTR> &acids,
                         const std::vector<ROMOL_SPTR> &amines) {
  std::size_t total_products = 0;
  for (unsigned int repeat = 0; repeat < BENCH_REPEATS; ++repeat) {
    for (const auto &acid : acids) {
      for (const auto &amine : amines) {
        MOL_SPTR_VECT reactants{acid, amine};
        total_products += run_Reactants(rxn, reactants).size();
      }
    }
  }
  return total_products;
}

std::size_t run_cached(const ChemicalReaction &rxn,
                       const std::vector<ROMOL_SPTR> &acids,
                       const std::vector<ROMOL_SPTR> &amines) {
  ReactantMatchCache cache;
  std::size_t total_products = 0;
  for (unsigned int repeat = 0; repeat < BENCH_REPEATS; ++repeat) {
    for (const auto &acid : acids) {
      for (const auto &amine : amines) {
        MOL_SPTR_VECT reactants{acid, amine};
        total_products += run_Reactants(rxn, reactants, cache).size();
      }
    }
  }
  return total_products;
}

}  // namespace

TEST_CASE("reaction enumeration matching", "[reaction][enumerate]") {
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
      "[C:1](=[O:2])[O;H1].[N:3]>>[C:1](=[O:2])[N:3]"));
  REQUIRE(rxn);
  rxn->initReactantMatchers();

  const auto acids = make_acid_pool(64);
  const auto amines = make_amine_pool(64);

  const std::size_t uncached_total = run_uncached(*rxn, acids, amines);
  const std::size_t cached_total = run_cached(*rxn, acids, amines);
  CHECK(cached_total == uncached_total);

  BENCHMARK("runReactants uncached") {
    std::size_t total_products = 0;
    for (unsigned int repeat = 0; repeat < BENCH_REPEATS; ++repeat) {
      for (const auto &acid : acids) {
        for (const auto &amine : amines) {
          MOL_SPTR_VECT reactants{acid, amine};
          total_products += run_Reactants(*rxn, reactants).size();
        }
      }
    }
    return total_products;
  };

  BENCHMARK("runReactants cached") {
    ReactantMatchCache cache;
    std::size_t total_products = 0;
    for (unsigned int repeat = 0; repeat < BENCH_REPEATS; ++repeat) {
      for (const auto &acid : acids) {
        for (const auto &amine : amines) {
          MOL_SPTR_VECT reactants{acid, amine};
          total_products += run_Reactants(*rxn, reactants, cache).size();
        }
      }
    }
    return total_products;
  };
}

TEST_CASE("reaction enumeration matching (recursive SMARTS)",
          "[reaction][enumerate]") {
  // Both reactant templates use recursive SMARTS ($(...)), so the
  // substructure matching is expensive relative to product assembly. This is
  // where caching repeated per-reagent matches pays off the most. Smaller
  // pools are used here because the recursive matching is costly.
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
      "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
      "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>"
      "[N:3]-[C:1]-[N+0:2]"));
  REQUIRE(rxn);
  rxn->initReactantMatchers();

  constexpr unsigned int POOL = 24;
  // isothiocyanates of varying chain length match the first template
  auto isothiocyanates = make_reagent_pool(POOL, "N=C=S");
  // simple primary amines match the second (restrictive) template
  std::vector<ROMOL_SPTR> amines;
  amines.reserve(POOL);
  for (unsigned int i = 0; i < POOL; ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles("N" + std::string(i + 1, 'C'));
    REQUIRE(mol);
    amines.push_back(ROMOL_SPTR(mol.release()));
  }

  // correctness: cached and uncached produce the same number of products
  std::size_t uncached_total = 0;
  std::size_t cached_total = 0;
  ReactantMatchCache checkCache;
  for (const auto &itc : isothiocyanates) {
    for (const auto &amine : amines) {
      MOL_SPTR_VECT reactants{itc, amine};
      uncached_total += run_Reactants(*rxn, reactants).size();
      cached_total += run_Reactants(*rxn, reactants, checkCache).size();
    }
  }
  CHECK(cached_total == uncached_total);

  BENCHMARK("runReactants uncached (recursive SMARTS)") {
    std::size_t total_products = 0;
    for (const auto &itc : isothiocyanates) {
      for (const auto &amine : amines) {
        MOL_SPTR_VECT reactants{itc, amine};
        total_products += run_Reactants(*rxn, reactants).size();
      }
    }
    return total_products;
  };

  BENCHMARK("runReactants cached (recursive SMARTS)") {
    ReactantMatchCache cache;
    std::size_t total_products = 0;
    for (const auto &itc : isothiocyanates) {
      for (const auto &amine : amines) {
        MOL_SPTR_VECT reactants{itc, amine};
        total_products += run_Reactants(*rxn, reactants, cache).size();
      }
    }
    return total_products;
  };
}

TEST_CASE("reaction enumeration matching (symmetric dedup)",
          "[reaction][enumerate]") {
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
      "[N;H2:1].[C:2](=[O:3])[O;H1]>>[N:1][C:2]=[O:3]"));
  REQUIRE(rxn);
  rxn->initReactantMatchers();

  constexpr unsigned int POOL = 16;
  std::vector<ROMOL_SPTR> diamines;
  diamines.reserve(POOL);
  for (unsigned int i = 0; i < POOL; ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles("N" + std::string(i + 1, 'C') +
                                              "N");
    REQUIRE(mol);
    diamines.push_back(ROMOL_SPTR(mol.release()));
  }
  const auto acids = make_acid_pool(POOL);

  auto collect_unique_product_smiles =
      [&rxn, &diamines, &acids](bool dedupeSymmetricMatches) {
        std::set<std::string> product_smiles;
        ReactantMatchCache cache;
        for (const auto &diamine : diamines) {
          for (const auto &acid : acids) {
            MOL_SPTR_VECT reactants{diamine, acid};
            const auto products = run_Reactants(*rxn, reactants, cache,
                                                dedupeSymmetricMatches);
            for (const auto &product_set : products) {
              for (const auto &product : product_set) {
                product_smiles.insert(MolToSmiles(*product));
              }
            }
          }
        }
        return product_smiles;
      };

  std::size_t dedupOff_total = 0;
  std::size_t dedupOn_total = 0;
  // the dedup flag is part of the cache key, so false/true entries never
  // collide even when the same cache is shared between both modes here
  ReactantMatchCache checkCache;
  for (const auto &diamine : diamines) {
    for (const auto &acid : acids) {
      MOL_SPTR_VECT reactants{diamine, acid};
      dedupOff_total += run_Reactants(*rxn, reactants, checkCache, false).size();
      dedupOn_total += run_Reactants(*rxn, reactants, checkCache, true).size();
    }
  }
  CHECK(dedupOn_total < dedupOff_total);
  CHECK(collect_unique_product_smiles(false) ==
        collect_unique_product_smiles(true));

  BENCHMARK("runReactants dedupe off") {
    ReactantMatchCache cache;
    std::size_t total_products = 0;
    for (const auto &diamine : diamines) {
      for (const auto &acid : acids) {
        MOL_SPTR_VECT reactants{diamine, acid};
        total_products += run_Reactants(*rxn, reactants, cache, false).size();
      }
    }
    return total_products;
  };

  BENCHMARK("runReactants dedupe on") {
    ReactantMatchCache cache;
    std::size_t total_products = 0;
    for (const auto &diamine : diamines) {
      for (const auto &acid : acids) {
        MOL_SPTR_VECT reactants{diamine, acid};
        total_products += run_Reactants(*rxn, reactants, cache, true).size();
      }
    }
    return total_products;
  };
}