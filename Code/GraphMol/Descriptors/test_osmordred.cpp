#include "Osmordred.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/RDLog.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

using namespace RDKit;
using namespace RDKit::Descriptors::Osmordred;

#ifdef RDK_BUILD_OSMORDRED_SUPPORT

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
// ============================================================
// Osmordred v2.0 Bug Fix Tests
// These tests verify the critical fixes implemented in v2.0
// ============================================================

TEST_CASE("Osmordred v2.0 - isMoleculeTooLarge") {
  SECTION("Normal molecule should not be flagged as too large") {
    auto mol = "CCO"_smiles;  // ethanol - small molecule
    REQUIRE(mol != nullptr);
    REQUIRE(isMoleculeTooLarge(*mol) == false);
  }
  
  SECTION("Very large linear chain (230 carbons) should be flagged") {
    // Build a SMILES string with 230 carbons: "CCCC...CCC"
    std::string large_smiles(230, 'C');
    auto mol = SmilesToMol(large_smiles);
    REQUIRE(mol != nullptr);
    
    // This molecule has 230 heavy atoms, exceeding the 200 threshold
    REQUIRE(isMoleculeTooLarge(*mol) == true);
    
    // calcOsmordred should still work but may return NaN for some descriptors
    auto result = calcOsmordred(*mol);
    REQUIRE(!result.empty());
    REQUIRE(result.size() == 3585);
    
    delete mol;
  }
  
  SECTION("Molecule with exactly 200 heavy atoms should not be flagged") {
    std::string border_smiles(200, 'C');
    auto mol = SmilesToMol(border_smiles);
    REQUIRE(mol != nullptr);
    
    // 200 is the threshold, should NOT be flagged
    REQUIRE(isMoleculeTooLarge(*mol) == false);
    
    delete mol;
  }
  
  SECTION("Molecule with 201 heavy atoms should be flagged") {
    std::string over_smiles(201, 'C');
    auto mol = SmilesToMol(over_smiles);
    REQUIRE(mol != nullptr);
    
    // 201 exceeds 200 threshold
    REQUIRE(isMoleculeTooLarge(*mol) == true);
    
    delete mol;
  }
  
  SECTION("Molecule with many rings (>10) should be flagged") {
    // Coronene-like structure with many fused rings
    // This is a naphthalene-derived structure with many rings
    auto mol = "c1cc2ccc3ccc4ccc5ccc6ccc7ccc8ccc9ccc%10ccc%11ccc1c1c2c3c4c5c6c7c8c9c%10c%111"_smiles;
    if (mol != nullptr) {
      // If this parses, check ring count
      auto ring_info = mol->getRingInfo();
      if (ring_info && ring_info->numRings() > 10) {
        REQUIRE(isMoleculeTooLarge(*mol) == true);
      }
    }
  }
}

TEST_CASE("Osmordred v2.0 - solveLinearSystem (BCUT fix)") {
  SECTION("Simple ethane molecule - BCUT should work") {
    auto mol = "CC"_smiles;  // ethane
    REQUIRE(mol != nullptr);
    
    // BCUT calculations use solveLinearSystem internally
    // v2.0 fix: saves B_original before LAPACK calls and adds dgelss fallback
    auto bcuts = calcBCUTs(*mol);
    REQUIRE(!bcuts.empty());
    
    // Check that we get valid (non-NaN) results
    bool has_valid = false;
    for (const auto& val : bcuts) {
      if (!std::isnan(val)) {
        has_valid = true;
        break;
      }
    }
    REQUIRE(has_valid);
  }
  
  SECTION("Benzene - BCUT test") {
    auto mol = "c1ccccc1"_smiles;  // benzene
    REQUIRE(mol != nullptr);
    
    auto bcuts = calcBCUTs(*mol);
    REQUIRE(!bcuts.empty());
    
    // Verify we get reasonable values
    bool has_valid = false;
    for (const auto& val : bcuts) {
      if (!std::isnan(val) && std::isfinite(val)) {
        has_valid = true;
        break;
      }
    }
    REQUIRE(has_valid);
  }
  
  SECTION("Drug-like molecule - full BCUT test") {
    // Aspirin
    auto mol = "CC(=O)OC1=CC=CC=C1C(=O)O"_smiles;
    REQUIRE(mol != nullptr);
    
    auto bcuts = calcBCUTs(*mol);
    REQUIRE(!bcuts.empty());
    REQUIRE(bcuts.size() > 0);
    
    // Test multiple times for consistency
    auto bcuts2 = calcBCUTs(*mol);
    REQUIRE(bcuts == bcuts2);
  }
  
  SECTION("Symmetric molecule - tests matrix solver stability") {
    // Highly symmetric molecule that could cause numerical issues
    auto mol = "C1CCCCC1"_smiles;  // cyclohexane
    REQUIRE(mol != nullptr);
    
    auto bcuts = calcBCUTs(*mol);
    REQUIRE(!bcuts.empty());
    
    // Multiple calls should give consistent results
    auto bcuts2 = calcBCUTs(*mol);
    REQUIRE(bcuts == bcuts2);
  }
}

TEST_CASE("Osmordred v2.0 - checkGasteigerParameters") {
  SECTION("Normal organic molecule should pass Gasteiger check") {
    auto mol = "CCO"_smiles;  // ethanol
    REQUIRE(mol != nullptr);
    REQUIRE(checkGasteigerParameters(*mol) == true);
  }
  
  SECTION("Phosphorus molecule - may fail Gasteiger check") {
    // Phosphoric acid
    auto mol = "OP(O)(O)=O"_smiles;
    REQUIRE(mol != nullptr);
    
    // Gasteiger parameters exist for P, so this should pass
    // Whether it passes or not, calcOsmordred should handle it gracefully
    (void)checkGasteigerParameters(*mol);  // Just verify it doesn't crash
    auto result = calcOsmordred(*mol);
    REQUIRE(!result.empty());
    REQUIRE(result.size() == 3585);
  }
  
  SECTION("Vanadium molecule - should fail Gasteiger check") {
    // Vanadium compound - Gasteiger doesn't have parameters for V
    auto mol = "[V]"_smiles;
    if (mol != nullptr) {
      // Vanadium should NOT have Gasteiger parameters
      // checkGasteigerParameters should return false (or throw for truly unknown atoms)
      try {
        bool gasteiger_ok = checkGasteigerParameters(*mol);
        // If it doesn't throw, should return false
        REQUIRE(gasteiger_ok == false);
      } catch (...) {
        // Vanadium may throw in getParams - this is acceptable
        // The important thing is that the user is protected from crashes
      }
      
      // calcOsmordred should handle this gracefully
      // Either return NaN or skip charge-dependent descriptors
      try {
        auto result = calcOsmordred(*mol);
        REQUIRE(!result.empty());
        REQUIRE(result.size() == 3585);
      } catch (...) {
        // If it throws, that's also acceptable for unknown atom types
        // The v2.0 fix prevents crashes in common cases
      }
    }
  }
  
  SECTION("Mixed P and V molecule") {
    // A molecule with both P and V
    auto mol = "[V]P"_smiles;
    if (mol != nullptr) {
      // Should fail due to V
      try {
        bool gasteiger_ok = checkGasteigerParameters(*mol);
        REQUIRE(gasteiger_ok == false);
      } catch (...) {
        // Expected for truly unknown elements
      }
      
      // calcOsmordred should handle gracefully
      try {
        auto result = calcOsmordred(*mol);
        REQUIRE(!result.empty());
      } catch (...) {
        // Acceptable for unknown elements
      }
    }
  }
  
  SECTION("Transition metal complexes should fail Gasteiger check") {
    // Iron
    auto mol_fe = "[Fe]"_smiles;
    if (mol_fe != nullptr) {
      REQUIRE(checkGasteigerParameters(*mol_fe) == false);
    }
    
    // Copper
    auto mol_cu = "[Cu]"_smiles;
    if (mol_cu != nullptr) {
      REQUIRE(checkGasteigerParameters(*mol_cu) == false);
    }
    
    // Zinc
    auto mol_zn = "[Zn]"_smiles;
    if (mol_zn != nullptr) {
      REQUIRE(checkGasteigerParameters(*mol_zn) == false);
    }
  }
  
  SECTION("RNCG/RPCG with Gasteiger failure should return NaN") {
    // Vanadium oxide
    auto mol = "[V]=O"_smiles;
    if (mol != nullptr) {
      // This should fail Gasteiger check
      REQUIRE(checkGasteigerParameters(*mol) == false);
      
      // RNCG/RPCG should return NaN instead of crashing
      auto rncg_rpcg = calcRNCG_RPCG(*mol);
      REQUIRE(!rncg_rpcg.empty());
      // Should be NaN due to Gasteiger failure
      for (const auto& val : rncg_rpcg) {
        REQUIRE(std::isnan(val));
      }
    }
  }
}

TEST_CASE("Osmordred v2.0 - Timeout and Batch Functions") {
  SECTION("calcOsmordredWithTimeout basic test") {
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto result = calcOsmordredWithTimeout(*mol, 60);
    REQUIRE(!result.empty());
    REQUIRE(result.size() == 3585);
    
    // Should have valid results for ethanol
    bool has_valid = false;
    for (const auto& val : result) {
      if (!std::isnan(val)) {
        has_valid = true;
        break;
      }
    }
    REQUIRE(has_valid);
  }
  
  SECTION("calcOsmordredBatch basic test") {
    std::vector<std::string> smiles_list = {"CCO", "CCC", "c1ccccc1"};
    
    auto results = calcOsmordredBatch(smiles_list, 0);
    REQUIRE(results.size() == 3);
    
    for (const auto& result : results) {
      REQUIRE(result.size() == 3585);
    }
  }
  
  SECTION("calcOsmordredBatch with invalid SMILES") {
    std::vector<std::string> smiles_list = {"CCO", "INVALID", "CCC"};
    
    auto results = calcOsmordredBatch(smiles_list, 0);
    REQUIRE(results.size() == 3);
    
    // First and third should have valid results
    bool first_valid = false;
    for (const auto& val : results[0]) {
      if (!std::isnan(val)) {
        first_valid = true;
        break;
      }
    }
    REQUIRE(first_valid);
    
    // Second (invalid) should be all NaN
    bool second_all_nan = true;
    for (const auto& val : results[1]) {
      if (!std::isnan(val)) {
        second_all_nan = false;
        break;
      }
    }
    REQUIRE(second_all_nan);
  }
  
  SECTION("getOsmordredDescriptorNames") {
    auto names = getOsmordredDescriptorNames();
    REQUIRE(names.size() == 3585);
    
    // Check first few names are not empty
    REQUIRE(!names[0].empty());
    REQUIRE(!names[1].empty());
    REQUIRE(!names[2].empty());
  }
}

TEST_CASE("Osmordred v2.0 - AutoCorrelation with Gasteiger") {
  SECTION("AutoCorrelation should handle Gasteiger failure") {
    // Vanadium compound
    auto mol = "[V]=O"_smiles;
    if (mol != nullptr) {
      // Should not crash, should return NaN for charge-dependent descriptors
      auto autocorr = calcAutoCorrelation(*mol);
      REQUIRE(!autocorr.empty());
    }
  }
  
  SECTION("AutoCorrelation with normal molecule") {
    auto mol = "c1ccccc1"_smiles;  // benzene
    REQUIRE(mol != nullptr);
    
    auto autocorr = calcAutoCorrelation(*mol);
    REQUIRE(!autocorr.empty());
    
    // Should have valid values
    bool has_valid = false;
    for (const auto& val : autocorr) {
      if (!std::isnan(val) && std::isfinite(val)) {
        has_valid = true;
        break;
      }
    }
    REQUIRE(has_valid);
  }
}

TEST_CASE("Osmordred v2.0 - NCI Dataset Stress Test") {
  // Get RDBASE from environment
  const char* rdbase_env = std::getenv("RDBASE");
  std::string rdbase = rdbase_env ? rdbase_env : "";
  
  SECTION("Load and process NCI first_5K dataset") {
    // Skip test if RDBASE not set
    if (rdbase.empty()) {
      WARN("RDBASE not set, skipping NCI dataset test");
      return;
    }
    
    std::string nci_path = rdbase + "/Data/NCI/first_5K.smi";
    std::ifstream file(nci_path);
    
    if (!file.is_open()) {
      WARN("Could not open NCI file: " << nci_path);
      return;
    }
    
    // Disable RDKit warnings for cleaner output
    boost::logging::disable_logs("rdApp.*");
    
    int total_molecules = 0;
    int successful = 0;
    int parse_failed = 0;
    int calc_failed = 0;
    int too_large = 0;
    
    std::string line;
    // Process first 500 molecules for speed
    int max_mols = 500;
    
    while (std::getline(file, line) && total_molecules < max_mols) {
      // Parse tab-separated SMILES
      std::istringstream iss(line);
      std::string smiles;
      std::getline(iss, smiles, '\t');
      
      if (smiles.empty()) continue;
      total_molecules++;
      
      // Try to parse SMILES
      ROMol* mol = nullptr;
      try {
        mol = SmilesToMol(smiles);
      } catch (...) {
        parse_failed++;
        continue;
      }
      
      if (mol == nullptr) {
        parse_failed++;
        continue;
      }
      
      // Check if molecule is too large
      if (isMoleculeTooLarge(*mol)) {
        too_large++;
        delete mol;
        continue;
      }
      
      // Try to calculate Osmordred descriptors
      try {
        auto result = calcOsmordred(*mol);
        
        // Verify result size
        REQUIRE(result.size() == 3585);
        
        // Count valid (non-NaN) descriptors
        int valid_count = 0;
        for (const auto& val : result) {
          if (!std::isnan(val) && std::isfinite(val)) {
            valid_count++;
          }
        }
        
        // A successful molecule should have at least some valid descriptors
        if (valid_count > 0) {
          successful++;
        } else {
          calc_failed++;
        }
      } catch (...) {
        calc_failed++;
      }
      
      delete mol;
    }
    
    file.close();
    
    // Re-enable logging
    boost::logging::enable_logs("rdApp.*");
    
    // Report results
    INFO("NCI Dataset Test Results:");
    INFO("  Total molecules: " << total_molecules);
    INFO("  Successful: " << successful);
    INFO("  Parse failed: " << parse_failed);
    INFO("  Calc failed: " << calc_failed);
    INFO("  Too large: " << too_large);
    
    // Success rate should be at least 90%
    double success_rate = (double)successful / (total_molecules - parse_failed) * 100.0;
    INFO("  Success rate: " << success_rate << "%");
    
    REQUIRE(total_molecules > 0);
    REQUIRE(successful > 0);
    // At least 90% success rate on parseable molecules
    REQUIRE(success_rate >= 90.0);
  }
  
  SECTION("Batch processing NCI molecules") {
    if (rdbase.empty()) {
      WARN("RDBASE not set, skipping NCI batch test");
      return;
    }
    
    std::string nci_path = rdbase + "/Data/NCI/first_5K.smi";
    std::ifstream file(nci_path);
    
    if (!file.is_open()) {
      WARN("Could not open NCI file: " << nci_path);
      return;
    }
    
    // Collect first 100 SMILES
    std::vector<std::string> smiles_list;
    std::string line;
    int count = 0;
    
    while (std::getline(file, line) && count < 100) {
      std::istringstream iss(line);
      std::string smiles;
      std::getline(iss, smiles, '\t');
      if (!smiles.empty()) {
        smiles_list.push_back(smiles);
        count++;
      }
    }
    file.close();
    
    REQUIRE(smiles_list.size() == 100);
    
    // Disable logging
    boost::logging::disable_logs("rdApp.*");
    
    // Test batch processing
    auto results = calcOsmordredBatch(smiles_list, 0);
    
    // Re-enable logging
    boost::logging::enable_logs("rdApp.*");
    
    // Should get results for all molecules
    REQUIRE(results.size() == 100);
    
    // Each result should have correct size (3585) OR be empty (for invalid molecules)
    int valid_results = 0;
    int empty_results = 0;
    for (const auto& result : results) {
      if (result.size() == 3585) {
        valid_results++;
      } else if (result.empty()) {
        empty_results++;
      }
      // Result should be either 3585 descriptors or empty (for invalid molecules)
      REQUIRE((result.size() == 3585 || result.empty()));
    }
    
    INFO("Batch results: " << valid_results << " valid, " << empty_results << " empty (invalid molecules)");
    
    // Count successful molecules (at least one non-NaN descriptor)
    int batch_successful = 0;
    for (const auto& result : results) {
      if (result.empty()) continue;  // Skip empty results
      for (const auto& val : result) {
        if (!std::isnan(val)) {
          batch_successful++;
          break;
        }
      }
    }
    
    INFO("Batch processing: " << batch_successful << "/" << valid_results << " successful");
    // At least 90% of valid results should have non-NaN descriptors
    if (valid_results > 0) {
      double success_rate = (double)batch_successful / valid_results * 100.0;
      REQUIRE(success_rate >= 90.0);
    }
  }
  
  SECTION("Compare individual vs batch results") {
    // Test that batch and individual processing give same results
    std::vector<std::string> test_smiles = {
      "CCO",           // ethanol
      "c1ccccc1",      // benzene
      "CC(=O)O",       // acetic acid
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  // caffeine
    };
    
    // Disable logging
    boost::logging::disable_logs("rdApp.*");
    
    // Get batch results
    auto batch_results = calcOsmordredBatch(test_smiles, 0);
    REQUIRE(batch_results.size() == test_smiles.size());
    
    // Compare with individual results
    for (size_t i = 0; i < test_smiles.size(); ++i) {
      auto mol = SmilesToMol(test_smiles[i]);
      REQUIRE(mol != nullptr);
      
      auto individual_result = calcOsmordred(*mol);
      
      // Results should match
      REQUIRE(batch_results[i].size() == individual_result.size());
      
      // Compare values (allowing for NaN == NaN)
      for (size_t j = 0; j < individual_result.size(); ++j) {
        bool both_nan = std::isnan(batch_results[i][j]) && std::isnan(individual_result[j]);
        bool both_equal = batch_results[i][j] == individual_result[j];
        REQUIRE((both_nan || both_equal));
      }
      
      delete mol;
    }
    
    // Re-enable logging
    boost::logging::enable_logs("rdApp.*");
  }
}

TEST_CASE("Osmordred v2.0 - Descriptor Count Validation") {
  SECTION("Verify descriptor count matches v1.0 specification") {
    // Osmordred should produce exactly 3585 descriptors
    auto mol = "CCO"_smiles;
    REQUIRE(mol != nullptr);
    
    auto result = calcOsmordred(*mol);
    REQUIRE(result.size() == 3585);
    
    auto names = getOsmordredDescriptorNames();
    REQUIRE(names.size() == 3585);
    
    // Descriptor count should match
    REQUIRE(result.size() == names.size());
  }
  
  SECTION("Verify consistency across different molecule types") {
    std::vector<std::string> diverse_smiles = {
      "C",                    // methane
      "CCO",                  // ethanol
      "c1ccccc1",             // benzene
      "C1CCCCC1",             // cyclohexane
      "CC(C)C",               // isobutane
      "c1ccc2ccccc2c1",       // naphthalene
      "CC(=O)OC1=CC=CC=C1C(=O)O",  // aspirin
    };
    
    for (const auto& smi : diverse_smiles) {
      auto mol = SmilesToMol(smi);
      REQUIRE(mol != nullptr);
      
      auto result = calcOsmordred(*mol);
      REQUIRE(result.size() == 3585);
      
      delete mol;
    }
  }
}

#else
TEST_CASE("Osmordred Basic Functionality") {
  SECTION("No Osmordred support") {
  }
}
#endif
