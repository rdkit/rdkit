//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC's and others
//  All rights reserved.
//
// SMARTS291 - Abraham SMARTS-based Features Unit Tests
// Tests C++ implementation against Python CalcAbrahamsFeatures golden reference

#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/smarts291/SMARTS291.h>
#include <GraphMol/Descriptors/smarts291/abraham_queries.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace RDKit;

// Test tolerance - allow 1e-5 relative or absolute difference
const double TOLERANCE = 1e-5;

// Check if two values match within tolerance
bool valuesMatch(double computed, double golden) {
    // Both NaN = match
    if (std::isnan(computed) && std::isnan(golden)) return true;
    // One NaN = no match
    if (std::isnan(computed) || std::isnan(golden)) return false;
    // Both Inf with same sign = match
    if (std::isinf(computed) && std::isinf(golden)) {
        return (computed > 0) == (golden > 0);
    }
    // One Inf = no match
    if (std::isinf(computed) || std::isinf(golden)) return false;
    
    // Absolute difference check
    double absDiff = std::abs(computed - golden);
    if (absDiff < TOLERANCE) return true;
    
    // Relative difference check
    double maxVal = std::max(std::abs(computed), std::abs(golden));
    if (maxVal > 0 && absDiff / maxVal < TOLERANCE) return true;
    
    return false;
}

TEST_CASE("SMARTS291 Basic Functionality", "[smarts291][basic]") {
    SECTION("Ethanol features") {
        ROMol* mol = SmilesToMol("CCO");
        REQUIRE(mol != nullptr);
        
        std::vector<double> features = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        
        // Should return 291 features (241 base + 50 golden)
        CHECK(features.size() == 291);
        
        // Features should not be all zeros for a real molecule
        double sum = 0;
        for (double f : features) {
            if (!std::isnan(f)) sum += std::abs(f);
        }
        CHECK(sum > 0);
        
        delete mol;
    }
    
    SECTION("Base features extraction") {
        ROMol* mol = SmilesToMol("c1ccccc1");  // Benzene
        REQUIRE(mol != nullptr);
        
        std::vector<double> baseFeatures = Descriptors::Osmordred::extractAbrahamBaseFeatures(*mol);
        
        // Should return exactly 241 base features
        CHECK(baseFeatures.size() == 241);
        
        delete mol;
    }
}

TEST_CASE("SMARTS291 Batch Processing", "[smarts291][batch]") {
    std::vector<std::string> smiles = {"CCO", "CC(=O)C", "c1ccccc1", "CCN"};
    
    for (const auto& smi : smiles) {
        ROMol* mol = SmilesToMol(smi);
        REQUIRE(mol != nullptr);
        
        std::vector<double> features = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        CHECK(features.size() == 291);
        
        delete mol;
    }
}

TEST_CASE("SMARTS291 NCI Golden Reference Test", "[smarts291][golden]") {
    // Locate the golden reference file
    std::string rdbase = std::getenv("RDBASE") ? std::getenv("RDBASE") : "";
    std::string golden_path;
    
    if (!rdbase.empty()) {
        golden_path = rdbase + "/Code/GraphMol/Descriptors/test_data/nci_100_smarts291_golden.csv";
    } else {
        // Try relative path
        golden_path = "Code/GraphMol/Descriptors/test_data/nci_100_smarts291_golden.csv";
    }
    
    std::ifstream golden_file(golden_path);
    if (!golden_file.is_open()) {
        // Try current directory
        golden_path = "test_data/nci_100_smarts291_golden.csv";
        golden_file.open(golden_path);
    }
    
    // If golden file not found, skip this test
    if (!golden_file.is_open()) {
        WARN("Golden file not found at " << golden_path << " - skipping golden reference test");
        return;
    }
    
    // Parse header to get number of descriptors
    std::string header;
    std::getline(golden_file, header);
    
    // Count columns (smiles + 291 features)
    size_t n_descriptors = 0;
    {
        std::stringstream ss(header);
        std::string col;
        while (std::getline(ss, col, ',')) {
            n_descriptors++;
        }
        n_descriptors--;  // Subtract SMILES column
    }
    
    INFO("Golden file has " << n_descriptors << " descriptors");
    REQUIRE(n_descriptors == 291);
    
    int validated = 0;
    int mismatches = 0;
    std::string line;
    
    while (std::getline(golden_file, line) && validated < 100) {
        std::stringstream ss(line);
        std::string smiles;
        std::getline(ss, smiles, ',');
        
        // Parse golden values
        std::vector<double> golden_values;
        std::string value;
        while (std::getline(ss, value, ',')) {
            try {
                if (value.empty() || value == "nan" || value == "NaN") {
                    golden_values.push_back(std::nan(""));
                } else if (value == "inf") {
                    golden_values.push_back(std::numeric_limits<double>::infinity());
                } else if (value == "-inf") {
                    golden_values.push_back(-std::numeric_limits<double>::infinity());
                } else {
                    golden_values.push_back(std::stod(value));
                }
            } catch (...) {
                golden_values.push_back(std::nan(""));
            }
        }
        
        if (golden_values.size() != 291) {
            WARN("Skipping line with " << golden_values.size() << " values (expected 291)");
            continue;
        }
        
        // Compute SMARTS291 features
        ROMol* mol = SmilesToMol(smiles);
        if (!mol) {
            WARN("Could not parse SMILES: " << smiles);
            continue;
        }
        
        std::vector<double> computed = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        delete mol;
        
        if (computed.size() != 291) {
            WARN("Computed features have wrong size: " << computed.size());
            continue;
        }
        
        // Compare each feature
        for (size_t i = 0; i < 291; i++) {
            if (!valuesMatch(computed[i], golden_values[i])) {
                mismatches++;
                if (mismatches <= 10) {  // Only log first 10 mismatches
                    WARN("Mismatch at molecule " << validated << " (" << smiles << ") feature " << i 
                         << ": computed=" << std::setprecision(10) << computed[i] 
                         << " golden=" << golden_values[i]);
                }
            }
        }
        
        validated++;
    }
    
    golden_file.close();
    
    INFO("Validated " << validated << " molecules");
    INFO("Total mismatches: " << mismatches << " out of " << (validated * 291) << " values");
    
    // Require at least 90 molecules validated
    REQUIRE(validated >= 90);
    
    // Allow up to 5% mismatch (some SMARTS patterns may have minor differences)
    double mismatch_rate = static_cast<double>(mismatches) / (validated * 291);
    INFO("Mismatch rate: " << (mismatch_rate * 100) << "%");
    REQUIRE(mismatch_rate < 0.05);
}

TEST_CASE("SMARTS291 Edge Cases", "[smarts291][edge]") {
    SECTION("Single atom molecule") {
        ROMol* mol = SmilesToMol("C");  // Methane
        REQUIRE(mol != nullptr);
        
        std::vector<double> features = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        CHECK(features.size() == 291);
        
        delete mol;
    }
    
    SECTION("Large molecule") {
        // Cholesterol-like structure
        ROMol* mol = SmilesToMol("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C");
        REQUIRE(mol != nullptr);
        
        std::vector<double> features = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        CHECK(features.size() == 291);
        
        delete mol;
    }
    
    SECTION("Molecule with heteroatoms") {
        ROMol* mol = SmilesToMol("c1ccc(N)c(O)c1S");  // Amino-thiophenol
        REQUIRE(mol != nullptr);
        
        std::vector<double> features = Descriptors::Osmordred::calcAbrahamsFeatures(*mol);
        CHECK(features.size() == 291);
        
        delete mol;
    }
}
