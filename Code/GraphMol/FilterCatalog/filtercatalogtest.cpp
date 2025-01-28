//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>

using namespace RDKit;
using namespace std;

struct IntPair {
  int first;
  int second;
};

template <class T>
void dump(std::string name, const T &v) {
  std::cerr << name << " = { " << std::endl;
  for (size_t i = 0; i < v.size(); ++i) {
    std::cerr << "\t" << v[i].first << "," << v[i].second << "}," << std::endl;
    ;
  }
  std::cerr << "}" << std::endl;
}

bool check(MatchVectType v, MatchVectType match) {
  dump("v", v);
  dump("match", match);
  for (size_t i = 0; i < v.size(); ++i) {
    if (v[i].first != match[i].first) {
      return false;
    }
    if (v[i].second != match[i].second) {
      return false;
    }
  }
  return true;
}

void testFilterCatalog() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing the filter catalog " << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    SmilesMolSupplier suppl(pathName + "pains.smi");

    FilterCatalogParams params;
    params.addCatalog(FilterCatalogParams::PAINS_A);
    params.addCatalog(FilterCatalogParams::PAINS_B);
    params.addCatalog(FilterCatalogParams::PAINS_C);

    FilterCatalog catalog(params);
    boost::scoped_ptr<ROMol> mol;
    const IntPair match1[10] = {{0, 23}, {1, 22}, {2, 20}, {3, 19}, {4, 25},
                                {5, 24}, {6, 18}, {7, 17}, {8, 16}, {9, 21}};
    MatchVectType matchvec1;
    for (auto i : match1) {
      matchvec1.push_back(std::make_pair(i.first, i.second));
    }

    const IntPair match2[13] = {{0, 11}, {1, 12},  {2, 13}, {3, 14}, {4, 15},
                                {5, 10}, {6, 9},   {7, 8},  {8, 7},  {9, 6},
                                {10, 5}, {11, 17}, {12, 16}};
    MatchVectType matchvec2;
    for (auto i : match2) {
      matchvec2.push_back(std::make_pair(i.first, i.second));
    }

    const IntPair match3[12] = {{0, 0}, {1, 1},  {2, 2},   {3, 4},
                                {4, 5}, {5, 6},  {6, 7},   {7, 8},
                                {8, 9}, {9, 14}, {10, 15}, {11, 16}};
    MatchVectType matchvec3;
    for (auto i : match3) {
      matchvec3.push_back(std::make_pair(i.first, i.second));
    }
    int count = 0;
    while (!suppl.atEnd()) {
      mol.reset(suppl.next());
      std::string name = mol->getProp<std::string>(common_properties::_Name);

      TEST_ASSERT(mol.get());
      if (catalog.hasMatch(*mol)) {
        std::cerr << "Warning: molecule failed filter " << std::endl;
      }
      // More detailed
      FilterCatalog::CONST_SENTRY entry = catalog.getFirstMatch(*mol);
      TEST_ASSERT(entry);
      if (entry) {
        std::cerr << "Warning: molecule failed filter: reason "
                  << entry->getDescription() << std::endl;
        switch (count) {
          case 0:
            TEST_ASSERT(entry->getDescription() == "hzone_phenol_A(479)");
            break;
          case 1:
            TEST_ASSERT(entry->getDescription() == "cyano_imine_B(17)");
            break;
          case 2:
            TEST_ASSERT(entry->getDescription() == "keto_keto_gamma(5)");
            break;
        }
        TEST_ASSERT(entry->getDescription() == name);

        // get the substructure atoms for visualization
        std::vector<FilterMatch> matches;
        if (entry->getFilterMatches(*mol, matches)) {
          for (std::vector<FilterMatch>::const_iterator it = matches.begin();
               it != matches.end(); ++it) {
            // Get the FilterMatcherBase that matched
            const FilterMatch &fm = (*it);
            boost::shared_ptr<FilterMatcherBase> matchingFilter =
                fm.filterMatch;

            // Get the matching atom indices
            const MatchVectType &vect = fm.atomPairs;
            switch (count) {
              case 0:
                TEST_ASSERT(check(vect, matchvec1));
                break;
              case 1:
                TEST_ASSERT(check(vect, matchvec2));
                break;
              case 2:
                TEST_ASSERT(check(vect, matchvec3));
                break;
            }

            // do something with these...
          }
        }
      }
      count++;
    }  // end while
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testFilterCatalogEntry() {
  SmartsMatcher *sm = new SmartsMatcher("Aromatic carbon chain");
  boost::shared_ptr<FilterMatcherBase> matcher(sm);
  TEST_ASSERT(!matcher->isValid());
  const int debugParse = 0;
  const bool mergeHs = true;
  ROMOL_SPTR pattern(SmartsToMol("c:c:c:c:c", debugParse, mergeHs));
  TEST_ASSERT(pattern.get() != nullptr);
  sm->setPattern(pattern);
  sm->setMinCount(1);
  FilterCatalogEntry entry("Bar", matcher);
  TEST_ASSERT(entry.getDescription() == "Bar");
  TEST_ASSERT(sm->getMinCount() == 1);
  TEST_ASSERT(sm->getMaxCount() == (unsigned int)-1);

  entry.setDescription("Foo");
  TEST_ASSERT(entry.getDescription() == "Foo");

  entry.setProp("foo", "foo");
  TEST_ASSERT(entry.getProp<std::string>("foo") == "foo");
  entry.setProp(std::string("bar"), "bar");
  TEST_ASSERT(entry.getProp<std::string>("bar") == "bar");

  RWMol *newM = SmilesToMol("c1ccccc1", 0, true);
  TEST_ASSERT(entry.hasFilterMatch(*newM));
  delete newM;
}

void testFilterCatalogThreadedRunner() {
  FilterCatalogParams params;
  params.addCatalog(FilterCatalogParams::PAINS_A);
  params.addCatalog(FilterCatalogParams::PAINS_B);
  params.addCatalog(FilterCatalogParams::PAINS_C);

  FilterCatalog catalog(params);

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/pains.smi";

  std::ifstream infile(pathName);
  std::vector<std::string> smiles;

  std::string line;
  int count = 0;
  while (std::getline(infile, line)) {
    if (count) {
      smiles.push_back(line);
    }
    count += 1;
  }
  TEST_ASSERT(smiles.size() == 3);

  int numThreads = 3;  // one per entry
  auto results = RunFilterCatalog(catalog, smiles, numThreads);
  TEST_ASSERT(results.size() == smiles.size());
  count = 0;
  for (auto &entries : results) {
    TEST_ASSERT(entries.size() > 0);
    switch (count) {
      case 0:
        TEST_ASSERT(entries[0]->getDescription() == "hzone_phenol_A(479)");
        break;
      case 1:
        TEST_ASSERT(entries[0]->getDescription() == "cyano_imine_B(17)");
        break;
      case 2:
        TEST_ASSERT(entries[0]->getDescription() == "keto_keto_gamma(5)");
        break;
    }
    count += 1;
  }
}

void testFilterCatalogCHEMBL() {
  FilterCatalogParams params;
  auto catalogs = {FilterCatalogParams::CHEMBL_BMS,
                   FilterCatalogParams::CHEMBL_LINT,
                   FilterCatalogParams::CHEMBL_Glaxo,
                   FilterCatalogParams::CHEMBL_MLSMR,
                   FilterCatalogParams::CHEMBL_Dundee,
                   FilterCatalogParams::CHEMBL_Inpharmatica,
                   FilterCatalogParams::CHEMBL_SureChEMBL};
  std::vector<unsigned> entries = {180, 57, 55, 116, 105, 91, 166};
  int i = 0;
  for (auto catalog : catalogs) {
    FilterCatalogParams params;
    params.addCatalog(catalog);
    FilterCatalog filtercat(params);
    TEST_ASSERT(entries[i++] == filtercat.getNumEntries());
  }

  {  // Test inpharmatica merge hydrogens
    const std::vector<std::pair<std::string, std::string>> tests = {
        {"CN=NC", "Filter5_azo"},
        {"CN=N", "Filter5_azo"},
        {"N=N", "Filter5_azo"},
        {"NN=N", ""},
        {"CC(=O)CCBr",
         "Filter26_alkyl_halide|Filter30_beta_halo_carbonyl|Filter75_alkyl_Br_I"},
        {"CC(=O)CCI",
         "Filter26_alkyl_halide|Filter30_beta_halo_carbonyl|Filter75_alkyl_Br_I"},
        {"C(=O)CCI",
         "Filter26_alkyl_halide|Filter30_beta_halo_carbonyl|Filter38_aldehyde|Filter75_alkyl_Br_I"},
        {"NC(=O)CCI", "Filter26_alkyl_halide|Filter75_alkyl_Br_I"},
        {"CC=NOS(=O)N", "Filter18_oxime_ester|Filter89_hydroxylamine"},
        {"NC=NOP(=O)N", "Filter89_hydroxylamine"},
        {"C=NOP(=O)N", "Filter18_oxime_ester|Filter89_hydroxylamine"},
        {"NC=NO", "Filter89_hydroxylamine"},
        {"CC(=O)NOC", ""},
        {"[Fe]C", "Filter9_metal"},
        {"[FeH]", "Filter9_metal"},
        {"[Fe]", ""},
        {"CN=O", "Filter12_nitroso"},
        {"N=O", "Filter12_nitroso"},
        {"S=P",
         "Filter13_PS_double_bond|Filter61_phosphor_halide_and_P_S_bond"},
        {"S=PC",
         "Filter13_PS_double_bond|Filter61_phosphor_halide_and_P_S_bond"},
        {"CC(=O)SC", "Filter29_thioester"},
        {"CC(=S)SC", "Filter29_thioester"},
        {"CC(=N)SC", "Filter29_thioester"},
        {"CC(=O)S", "Filter29_thioester|Filter74_thiol"},
        {"CC(=S)S", "Filter29_thioester|Filter74_thiol"},
        {"CC(=N)S", "Filter29_thioester|Filter74_thiol"},
        {"SO", "Filter31_so_bond|Filter74_thiol"},
        {"SOC", "Filter31_so_bond|Filter74_thiol"},
        {"c1csocc1", ""},
        {"OO", "Filter32_oo_bond"},
        {"COO", "Filter32_oo_bond"},
        {"c1coocc1", ""},
        {"CC(=O)C=C", "Filter44_michael_acceptor2"},
        {"C(=O)C=C", "Filter38_aldehyde|Filter44_michael_acceptor2"},
        {"OC(=O)C=C", "Filter44_michael_acceptor2"},
        {"C1(=O)C=CC(=O)C=C1", "Filter53_para_quinones"},
        {"CIC", "Filter49_halogen|Filter75_alkyl_Br_I"},
        {"CI(C)C", "Filter49_halogen|Filter75_alkyl_Br_I"}};
    FilterCatalogParams params;
    params.addCatalog(FilterCatalogParams::CHEMBL_Inpharmatica);
    FilterCatalog catalog(params);
    for (auto &test : tests) {
      std::unique_ptr<RWMol> mol(SmilesToMol(test.first));
      std::string matches;
      for (auto &match : catalog.getMatches(*mol)) {
        if (matches.size()) matches += "|";
        matches += match->getDescription();
      }

      TEST_ASSERT(matches == test.second);
    }
  }
  {
    // Test Lint merge hydrogens
    const std::vector<std::pair<std::string, std::string>> tests = {
        {"CN=C", ""},
        {"CN=CC", "acyclic imines"},
        {"CN=C(C)", "acyclic imines"},
        {"C1N=CC1", ""}};
    FilterCatalogParams params;
    params.addCatalog(FilterCatalogParams::CHEMBL_LINT);
    FilterCatalog catalog(params);
    for (auto &test : tests) {
      std::unique_ptr<RWMol> mol(SmilesToMol(test.first));
      std::string matches;
      for (auto &match : catalog.getMatches(*mol)) {
        if (matches.size()) matches += "|";
        matches += match->getDescription();
      }
      TEST_ASSERT(matches == test.second);
    }
  }
}

int main() {
  RDLog::InitLogs();
  // boost::logging::enable_logs("rdApp.debug");

  testFilterCatalog();
  testFilterCatalogEntry();
  testFilterCatalogThreadedRunner();
  testFilterCatalogCHEMBL();
  return 0;
}
