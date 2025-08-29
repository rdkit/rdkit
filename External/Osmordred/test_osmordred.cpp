#include "Osmordred.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <iostream>

using namespace RDKit;

TEST_CASE("Osmordred Basic Functionality") {
  SECTION("Simple molecule tests") {
    auto mol = "CCO"_smiles; // ethanol
    REQUIRE(mol != nullptr);
    
    // Test basic atom counts
    auto atom_counts = calcAtomCounts(*mol);
    REQUIRE(!atom_counts.empty());
    REQUIRE(atom_counts[0] == 9); // Total atoms
    REQUIRE(atom_counts[1] == 3); // Heavy atoms
    REQUIRE(atom_counts[2] == 0); // Spiro atoms
    REQUIRE(atom_counts[3] == 0); // Bridgehead atoms
    
    REQUIRE(atom_counts[5] == 6); // Hydrogen count
    REQUIRE(atom_counts[6] == 0); // Boron count
    REQUIRE(atom_counts[7] == 2); // Carbon count
    REQUIRE(atom_counts[8] == 0); // Nitrogen count
    REQUIRE(atom_counts[9] == 1); // Oxygen count
    
    // Test bond counts
    auto bond_counts = calcBondCounts(*mol);
    REQUIRE(!bond_counts.empty());
    REQUIRE(bond_counts[0] == 8); // Total bonds
    REQUIRE(bond_counts[1] == 2); // Heavy bonds
    REQUIRE(bond_counts[2] == 8); // Single bonds (including H bonds)
    REQUIRE(bond_counts[3] == 0); // Double bonds
    REQUIRE(bond_counts[4] == 0); // Triple bonds
  }
  
  SECTION("Aromatic molecule tests") {
    auto mol = "c1ccccc1"_smiles; // benzene
    REQUIRE(mol != nullptr);
    
    // Test aromatic counts
    auto aromatic_counts = calcAromatic(*mol);
    REQUIRE(!aromatic_counts.empty());
    REQUIRE(aromatic_counts[0] == 6); // Aromatic atoms
    REQUIRE(aromatic_counts[1] == 6); // Aromatic bonds
    
    // Test ring count
    auto ring_counts = calcRingCount(*mol);
    REQUIRE(!ring_counts.empty());
    REQUIRE(ring_counts[0] == 1); // Number of rings
  }
  
  SECTION("Complex molecule tests") {
    auto mol = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"_smiles;
    REQUIRE(mol != nullptr);
    
    // Test that all functions return non-empty results
    REQUIRE(!calcABCIndex(*mol).empty());
    REQUIRE(!calcAcidBase(*mol).empty());
    REQUIRE(!calcAromatic(*mol).empty());
    REQUIRE(!calcAtomCounts(*mol).empty());
    REQUIRE(!calcBalabanJ(*mol).empty());
    REQUIRE(!calcBertzCT(*mol).empty());
    REQUIRE(!calcBondCounts(*mol).empty());
    REQUIRE(!calcVertexAdjacencyInformation(*mol).empty());
    REQUIRE(!calcWeight(*mol).empty());
    REQUIRE(!calcWienerIndex(*mol).empty());
    REQUIRE(!calcVdwVolumeABC(*mol).empty());
    REQUIRE(!calcTopoPSA(*mol).empty());
    REQUIRE(!calcSLogP(*mol).empty());
    REQUIRE(!calcHydrogenBond(*mol).empty());
    REQUIRE(!calcLogS(*mol).empty());
    REQUIRE(!calcLipinskiGhose(*mol).empty());
    REQUIRE(!calcMcGowanVolume(*mol).empty());
    REQUIRE(!calcPolarizability(*mol).empty());
    REQUIRE(!calcRotatableBond(*mol).empty());
    REQUIRE(!calcFragmentComplexity(*mol).empty());
    REQUIRE(!calcConstitutional(*mol).empty());
    REQUIRE(!calcTopologicalIndex(*mol).empty());
  }
}

TEST_CASE("Osmordred Matrix Descriptors") {
  SECTION("Distance matrix descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto dist_descs = calcDistMatrixDescs(*mol);
    REQUIRE(!dist_descs.empty());
    
    auto dist_descs_l = calcDistMatrixDescsL(*mol);
    REQUIRE(!dist_descs_l.empty());    
  }
  
  SECTION("Adjacency matrix descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto adj_descs = calcAdjMatrixDescs(*mol);
    REQUIRE(!adj_descs.empty());
    
    auto adj_descs_l = calcAdjMatrixDescsL(*mol);
    REQUIRE(!adj_descs_l.empty());
    
  }
  
  SECTION("Detour matrix descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto detour_descs = calcDetourMatrixDescs(*mol);
    REQUIRE(!detour_descs.empty());
    
    auto detour_descs_l = calcDetourMatrixDescsL(*mol);
    REQUIRE(!detour_descs_l.empty());
  }
  
  SECTION("Barysz matrix descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto barysz_descs = calcBaryszMatrixDescs(*mol);
    REQUIRE(!barysz_descs.empty());
    
    auto barysz_descs_l = calcBaryszMatrixDescsL(*mol);
    REQUIRE(!barysz_descs_l.empty());
  }
}

TEST_CASE("Osmordred Carbon Types") {
  SECTION("Carbon type descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto carbon_types = calcCarbonTypes(*mol);
    REQUIRE(!carbon_types.empty());    
  }
}

TEST_CASE("Osmordred EState Descriptors") {
  SECTION("EState descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto estate_descs = calcEStateDescs(*mol);
    REQUIRE(!estate_descs.empty());
    
    // Test extended version
    auto estate_descs_ext = calcEStateDescs(*mol, true);
    REQUIRE(!estate_descs_ext.empty());
  }
}

TEST_CASE("Osmordred Chi Descriptors") {
  SECTION("Chi descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto chi_descs = calcAllChiDescriptors(*mol);
    REQUIRE(!chi_descs.empty());
    
    auto chipath = calcChipath(*mol);
    REQUIRE(!chipath.empty());
    
    auto chichain = calcChichain(*mol);
    REQUIRE(!chichain.empty());
    
    auto chicluster = calcChicluster(*mol);
    REQUIRE(!chicluster.empty());
    
    auto chipathcluster = calcChipathcluster(*mol);
    REQUIRE(!chipathcluster.empty());
  }
}

TEST_CASE("Osmordred Count Functions") {
  SECTION("Count functions") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto acidic_count = calcAcidicGroupCount(*mol);
    REQUIRE(acidic_count >= 0);
    
    auto basic_count = calcBasicGroupCount(*mol);
    REQUIRE(basic_count >= 0);
    
    auto aromatic_atoms = countAromaticAtoms(*mol);
    REQUIRE(aromatic_atoms >= 0);
    
    auto aromatic_bonds = countAromaticBonds(*mol);
    REQUIRE(aromatic_bonds >= 0);
  }
}

TEST_CASE("Osmordred Information Content") {
  SECTION("Information content descriptors") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto info_content = calcInformationContent(*mol);
    REQUIRE(!info_content.empty());
    
    // Test with different max radius
    auto info_content_r3 = calcInformationContent(*mol, 3);
    REQUIRE(!info_content_r3.empty());
    
    auto info_content_r7 = calcInformationContent(*mol, 7);
    REQUIRE(!info_content_r7.empty());
  }
}

TEST_CASE("Osmordred Edge Cases") {
  SECTION("Empty molecule") {
    auto mol = ""_smiles;
    // This should return nullptr for invalid SMILES
    if (mol != nullptr) {
      // If molecule is valid, test that functions handle it gracefully
      auto atom_counts = calcAtomCounts(*mol);
      REQUIRE((atom_counts.empty() || atom_counts[0] == 0));
    }
  }
  
  SECTION("Single atom molecule") {
    auto mol = "[He]"_smiles; // Helium atom
    if (mol != nullptr) {
      auto atom_counts = calcAtomCounts(*mol);
      REQUIRE(!atom_counts.empty());
    }
  }
  
  SECTION("Large molecule") {
    auto mol = "CC(C)(C)c1ccc(C(C)(C)C)cc1"_smiles; // t-butyl benzene
    REQUIRE(mol != nullptr);
    
    // Test that all functions work with larger molecules
    REQUIRE(!calcABCIndex(*mol).empty());
    REQUIRE(!calcBertzCT(*mol).empty());
    REQUIRE(!calcWienerIndex(*mol).empty());
    REQUIRE(!calcVdwVolumeABC(*mol).empty());
    REQUIRE(!calcTopoPSA(*mol).empty());
    REQUIRE(!calcSLogP(*mol).empty());
    REQUIRE(!calcLogS(*mol).empty());
    REQUIRE(!calcMcGowanVolume(*mol).empty());
    REQUIRE(!calcPolarizability(*mol).empty());
    REQUIRE(!calcRotatableBond(*mol).empty());
    REQUIRE(!calcFragmentComplexity(*mol).empty());
    REQUIRE(!calcConstitutional(*mol).empty());
    REQUIRE(!calcTopologicalIndex(*mol).empty());
  }
}

TEST_CASE("Osmordred Validation Tests") {
  SECTION("Consistency checks") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    // Test that multiple calls return consistent results
    auto atom_counts1 = calcAtomCounts(*mol);
    auto atom_counts2 = calcAtomCounts(*mol);
    REQUIRE(atom_counts1 == atom_counts2);
    
    auto bond_counts1 = calcBondCounts(*mol);
    auto bond_counts2 = calcBondCounts(*mol);
    REQUIRE(bond_counts1 == bond_counts2);
    
    auto wiener1 = calcWienerIndex(*mol);
    auto wiener2 = calcWienerIndex(*mol);
    REQUIRE(wiener1 == wiener2);
  }
  
  SECTION("Physical property validation") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    // Test that physical properties are reasonable
    auto mw = calcWeight(*mol);
    REQUIRE(!mw.empty());
    REQUIRE(mw[0] > 0.0); // Molecular weight should be positive
    
    auto vdw_vol = calcVdwVolumeABC(*mol);
    REQUIRE(!vdw_vol.empty());
    REQUIRE(vdw_vol[0] > 0.0); // Volume should be positive
    
    auto mcgowan_vol = calcMcGowanVolume(*mol);
    REQUIRE(!mcgowan_vol.empty());
    REQUIRE(mcgowan_vol[0] > 0.0); // McGowan volume should be positive
  }
}
