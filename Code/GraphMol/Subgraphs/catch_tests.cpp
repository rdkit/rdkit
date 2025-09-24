//
//  Copyright (C) 2023 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <iostream>
#include <algorithm>

using namespace RDKit;

TEST_CASE("shortestPathsOnly") {
  SECTION("basics, findAllPathsOfLengthN") {
    auto m = "CC1CCC1"_smiles;
    REQUIRE(m);
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;
    auto ps = findAllPathsOfLengthN(*m, 4);
    CHECK(ps.size() == 3);

    ps = findAllPathsOfLengthN(*m, 4, useBonds, useHs, rootedAt,
                               onlyShortestPaths);
    CHECK(ps.empty());

    ps = findAllPathsOfLengthN(*m, 3);
    CHECK(ps.size() == 6);
    std::cerr << "---------" << std::endl;
    ps = findAllPathsOfLengthN(*m, 3, useBonds, useHs, rootedAt,
                               onlyShortestPaths);
    CHECK(ps.size() == 2);
  }
  SECTION("basics, findAllPathsOfLengthsMtoN") {
    auto m = "CC1CCC1"_smiles;
    REQUIRE(m);
    bool useBonds = true;
    bool useHs = false;
    int rootedAt = -1;
    bool onlyShortestPaths = true;
    auto ps = findAllPathsOfLengthsMtoN(*m, 3, 4);
    CHECK(ps.size() == 2);
    CHECK(ps[3].size() == 6);
    CHECK(ps[4].size() == 3);

    ps = findAllPathsOfLengthsMtoN(*m, 3, 4, useBonds, useHs, rootedAt,
                                   onlyShortestPaths);
    CHECK(ps.size() == 1);
    CHECK(ps[3].size() == 2);
  }
}

// helper api to get test data for copyMolSubset
struct SelectedComponents {
  std::vector<bool> selected_atoms;
  std::vector<bool> selected_bonds;
};

// helper api to get test mol for copyMolSubset api.
[[nodiscard]] static std::unique_ptr<RDKit::RWMol> get_test_mol()
{
  std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol("CCCCCCCCCCCCCCC")};
  for (auto& atom : mol->atoms()) {
    atom->setProp("orig_idx", atom->getIdx());
  }

  for (auto& bond : mol->bonds()) {
    bond->setProp("orig_idx", bond->getIdx());
  }

  return mol;
}

// Helper api to get the included atoms and bonds from test atom indices.
[[nodiscard]] static SelectedComponents
get_selected_components(::RDKit::RWMol& mol,
		                        const std::vector<unsigned int>& atom_ids)
{
  const auto num_atoms = mol.getNumAtoms();
  std::vector<bool> selected_atoms(num_atoms);

  for (auto& atom_idx : atom_ids) {
    if (atom_idx < num_atoms) {
      selected_atoms[atom_idx] = true;
    }
  }

  std::vector<bool> selected_bonds(mol.getNumBonds());
  for (auto& bond : mol.bonds()) {
    if (selected_atoms[bond->getBeginAtomIdx()] &&
      selected_atoms[bond->getEndAtomIdx()]) {
      selected_bonds[bond->getIdx()] = true;
    }
  }

  return {std::move(selected_atoms), std::move(selected_bonds)};
}

// This test makes sure we correctly extract atoms
TEST_CASE("test_extract_atoms", "[copyMolSubset]") {
  auto selected_atoms = GENERATE(
    // unique values
    std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12},
    // duplicate values
    std::vector<unsigned int>{0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
    // values outside of atom indices
    std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12, 100, 200, 300}
  );

  std::vector<unsigned int> expected_atoms{0, 2, 4, 6, 8, 10, 12};

  auto mol = get_test_mol();
  auto extracted_mol = copyMolSubset(*mol, selected_atoms);
  REQUIRE(extracted_mol->getNumAtoms() == expected_atoms.size());

  std::vector<unsigned int> extracted_atoms;
  for (auto& atom : extracted_mol->atoms()) {
    extracted_atoms.push_back(atom->template getProp<unsigned int>("orig_idx"));
  }

  CHECK(extracted_atoms == expected_atoms);
}

// This test makes sure we correctly extract atoms
TEST_CASE("test_extract_bonds", "[copyMolSubset]")
{
  auto test_mol = get_test_mol();

  for (auto& bond : test_mol->bonds()) {
    bond->setProp("test_prop", true);
  }

  for (auto& bond : test_mol->bonds()) {
    auto begin_idx = bond->getBeginAtomIdx();
    auto end_idx = bond->getEndAtomIdx();
    auto m = copyMolSubset(*test_mol, {begin_idx, end_idx});

    REQUIRE(m->getNumBonds() == 1);
    CHECK(m->getBondWithIdx(0)->getProp<bool>("test_prop") == true);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(0)->getProp<unsigned int>("orig_idx") == begin_idx);
    CHECK(m->getAtomWithIdx(1)->getProp<unsigned int>("orig_idx") == end_idx);
  }
}

// This test makes sure we correctly extract substance groups
TEST_CASE("test_extract_substance_groups", "[copyMolSubset]") {
  auto mol = get_test_mol();
  ::RDKit::SubstanceGroup sgroup{mol.get(), "COP"};

  auto test_sgroup_atoms = GENERATE(
     std::vector<unsigned int>{},
     std::vector<unsigned int>{0, 1, 2, 3, 4},
     std::vector<unsigned int>{9, 10, 11}
  );
  sgroup.setAtoms(test_sgroup_atoms);

  auto test_sgroup_bonds = GENERATE(
     std::vector<unsigned int>{},
     std::vector<unsigned int>{0, 1, 2},
     std::vector<unsigned int>{3, 4, 5}
  );
  sgroup.setBonds(test_sgroup_bonds);

  auto test_sgroup_patoms = GENERATE(
     std::vector<unsigned int>{},
     std::vector<unsigned int>{3, 4},
     std::vector<unsigned int>{5, 6}
  );
  sgroup.setParentAtoms(test_sgroup_patoms);

  ::RDKit::addSubstanceGroup(*mol, std::move(sgroup));

  auto has_selected_components = [&](auto& components, auto& ref_bitset) {
   return components.empty() ||
     std::ranges::all_of(components, [&](auto& idx) {
         return idx < ref_bitset.size() && ref_bitset[idx];
       });
  };

  auto test_selected_atoms = GENERATE(
     std::vector<unsigned int>{0, 1, 2, 3, 4, 5},
     std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12},
     std::vector<unsigned int>{0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
     std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12, 100, 200, 300}
  );

  auto extracted_mol = copyMolSubset(*mol, test_selected_atoms);

  auto [selected_atoms, selected_bonds] = get_selected_components(*mol, test_selected_atoms);
  // sgroup should only be copied if all components are selected
  auto flag = ::RDKit::getSubstanceGroups(*extracted_mol).size() == 1;
  REQUIRE(flag ==
             (has_selected_components(test_sgroup_atoms, selected_atoms) &&
              has_selected_components(test_sgroup_patoms, selected_atoms) &&
              has_selected_components(test_sgroup_bonds, selected_bonds)));

  // now make sure we copied the correct components
  if (flag) {
      auto& extracted_sgroup = ::RDKit::getSubstanceGroups(*extracted_mol)[0];
      for (auto& idx : extracted_sgroup.getAtoms()) {
          auto atom = extracted_mol->getAtomWithIdx(idx);
          CHECK(
              selected_atoms[atom->template getProp<unsigned int>("orig_idx")] ==
              true);
      }

      for (auto& idx : extracted_sgroup.getParentAtoms()) {
          auto atom = extracted_mol->getAtomWithIdx(idx);
          CHECK(
              selected_atoms[atom->template getProp<unsigned int>("orig_idx")] ==
              true);
      }

      for (auto& idx : extracted_sgroup.getBonds()) {
          auto bond = extracted_mol->getBondWithIdx(idx);
          CHECK(
              selected_bonds[bond->template getProp<unsigned int>("orig_idx")] ==
              true);
      }
  }
}

// This test makes sure we correctly extract stereo groups
TEST_CASE("test_extract_stereo_groups", "[copyMolSubset]") {
  auto mol = get_test_mol();

  auto test_stereo_group_atoms = GENERATE(
     std::vector<unsigned int>{},
     std::vector<unsigned int>{0, 1, 2, 3, 4},
     std::vector<unsigned int>{9, 10, 11}
  );

  std::vector<::RDKit::Atom*> sg_atoms;
  for (auto& idx : test_stereo_group_atoms) {
      sg_atoms.push_back(mol->getAtomWithIdx(idx));
  }

  auto test_stereo_group_bonds = GENERATE(
     std::vector<unsigned int>{},
     std::vector<unsigned int>{0, 1, 2},
     std::vector<unsigned int>{3, 4, 5}
  );

  std::vector<::RDKit::Bond*> sg_bonds;
  for (auto& idx : test_stereo_group_bonds) {
      sg_bonds.push_back(mol->getBondWithIdx(idx));
  }

  ::RDKit::StereoGroup stereo_group{::RDKit::StereoGroupType::STEREO_ABSOLUTE,
                                    std::move(sg_atoms), std::move(sg_bonds)};
  mol->setStereoGroups({std::move(stereo_group)});

  auto has_selected_components = [&](auto& components, auto& ref_bitset) {
      return components.empty() ||
             std::ranges::all_of(components, [&](auto& idx) {
                     return idx < ref_bitset.size() && ref_bitset[idx];
                 });
  };

  auto test_selected_atoms = GENERATE(
     std::vector<unsigned int>{0, 1, 2, 3, 4, 5},
     std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12},
     std::vector<unsigned int>{0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
     std::vector<unsigned int>{0, 2, 4, 6, 8, 10, 12, 100, 200, 300}
  );

  auto extracted_mol = copyMolSubset(*mol, test_selected_atoms);

  auto [selected_atoms, selected_bonds] =
      get_selected_components(*mol, test_selected_atoms);

  // stereo group should only be copied if all components are selected
  auto flag = extracted_mol->getStereoGroups().size() == 1;
  REQUIRE(flag ==
             (has_selected_components(test_stereo_group_atoms, selected_atoms) &&
              has_selected_components(test_stereo_group_bonds, selected_bonds)));

  // now make sure we copied the correct components
  if (flag) {
      auto& extracted_stereo_group = extracted_mol->getStereoGroups()[0];
      for (auto& atom : extracted_stereo_group.getAtoms()) {
          CHECK(
              selected_atoms[atom->template getProp<int>("orig_idx")] ==
              true);
      }

      for (auto& bond : extracted_stereo_group.getBonds()) {
          CHECK(
              selected_bonds[bond->template getProp<int>("orig_idx")] ==
              true);
      }
  }
}
