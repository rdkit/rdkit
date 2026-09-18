// CLEAN Abraham Integration - Uses new query system
// Generated to match trained model EXACTLY
// EACH MODEL TYPE USES ITS OWN GOLDEN FEATURES!

#include "SMARTS291.h"
#include "abraham_queries.h"
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/RDThreads.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <future>
#include <memory>
#include <limits>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

std::vector<double> extractAbrahamBaseFeatures(const RDKit::ROMol& mol) {
    std::vector<double> features;
    features.reserve(241);
    
    auto queries = GetQueriesAbrahamBaseFeatures();
    size_t max_queries = (queries.size() > 241) ? 241 : queries.size();
    
    RDKit::RWMol mol_rw(mol);
    
    for (size_t i = 0; i < max_queries; ++i) {
        if (queries[i]) {
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol_rw, *queries[i], matches, true);
            features.push_back(static_cast<double>(matches.size()));
        } else {
            features.push_back(0.0);
        }
    }
    
    while (features.size() < 241) {
        features.push_back(0.0);
    }

    return features;
}

std::vector<double> generateGoldenFeaturesA(const std::vector<double>& baseFeatures) {
    std::vector<double> goldenFeatures;
    goldenFeatures.reserve(50);

    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[89] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[29] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[91] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[108] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[29] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[39] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[54] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[0] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[54] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[34] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[18] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[90] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[90] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[18] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[89] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[55] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[108] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[2] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[39] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[91] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[51] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[9] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[107] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[51] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[13] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[93] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[9] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[55] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[93] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[11] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[52] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[92] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[92] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[53] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[2] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[53] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[52] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[56] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[34] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[56] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[59] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[85] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[85] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[38] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[110] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[110] / baseFeatures[127] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[60] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[59] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[3] / baseFeatures[113] : 0.0);
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[84] / baseFeatures[127] : 0.0);

    return goldenFeatures;
}

std::vector<double> calcAbrahamsFeatures(const RDKit::ROMol& mol) {
    std::vector<double> baseFeatures = extractAbrahamBaseFeatures(mol);
    std::vector<double> goldenA = generateGoldenFeaturesA(baseFeatures);
    
    std::vector<double> allFeatures = baseFeatures;
    allFeatures.insert(allFeatures.end(), goldenA.begin(), goldenA.end());
    
    return allFeatures;
}

}  // namespace Osmordred

namespace SMARTS291 {

bool hasSMARTS291Support() {
    return true;
}

std::vector<std::string> getBaseFeatureNames() {
    std::vector<std::string> names;
    names.reserve(241);
    for (int i = 0; i < 241; ++i) {
        names.push_back("base_" + std::to_string(i));
    }
    return names;
}

std::vector<std::string> getGoldenFeatureNames(char param) {
    std::vector<std::string> names;
    names.reserve(50);
    for (int i = 0; i < 50; ++i) {
        names.push_back(std::string("golden_") + param + "_" + std::to_string(i));
    }
    return names;
}

std::vector<std::string> getSMARTS291FeatureNames(char param) {
    std::vector<std::string> names = getBaseFeatureNames();
    std::vector<std::string> golden = getGoldenFeatureNames(param);
    names.insert(names.end(), golden.begin(), golden.end());
    return names;
}

std::vector<std::vector<double>> extractSMARTS291Batch(
    const std::vector<std::string>& smiles_list, char /*param*/, int n_jobs) {
    
    std::vector<std::vector<double>> results;
    results.reserve(smiles_list.size());
    
    unsigned int nThreads = getNumThreadsToUse(n_jobs);
    
    if (nThreads <= 1 || smiles_list.size() < 10) {
        for (const auto& smi : smiles_list) {
            ROMol* mol = SmilesToMol(smi);
            if (mol) {
                results.push_back(Osmordred::calcAbrahamsFeatures(*mol));
                delete mol;
            } else {
                results.push_back(std::vector<double>(291, 0.0));
            }
        }
        return results;
    }
    
    std::vector<std::future<std::vector<double>>> futures;
    futures.reserve(smiles_list.size());
    
    for (const auto& smi : smiles_list) {
        futures.emplace_back(std::async(std::launch::async, [smi]() {
            try {
                ROMol* mol = SmilesToMol(smi);
                if (mol) {
                    std::vector<double> feats = Osmordred::calcAbrahamsFeatures(*mol);
                    delete mol;
                    return feats;
                }
            } catch (...) {}
            return std::vector<double>(291, 0.0);
        }));
    }
    
    for (auto& f : futures) {
        results.push_back(f.get());
    }
    
    return results;
}

std::vector<std::vector<double>> extractSMARTS291FromMolsBatch(
    const std::vector<const ROMol*>& mols, char /*param*/, int n_jobs) {
    
    const double kNaN = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> nanRow(291, kNaN);
    
    std::vector<std::vector<double>> results;
    results.reserve(mols.size());
    
    unsigned int nThreads = getNumThreadsToUse(n_jobs);
    
    if (nThreads <= 1 || mols.size() < 10) {
        for (const auto* mol : mols) {
            if (mol) {
                try {
                    results.push_back(Osmordred::calcAbrahamsFeatures(*mol));
                } catch (...) {
                    results.push_back(nanRow);
                }
            } else {
                results.push_back(nanRow);
            }
        }
        return results;
    }
    
    std::vector<std::future<std::vector<double>>> futures;
    futures.reserve(mols.size());
    
    for (const auto* mol : mols) {
        futures.emplace_back(std::async(std::launch::async, [mol, &nanRow]() {
            if (mol) {
                try {
                    return Osmordred::calcAbrahamsFeatures(*mol);
                } catch (...) {}
            }
            return nanRow;
        }));
    }
    
    for (auto& f : futures) {
        results.push_back(f.get());
    }
    
    return results;
}

}  // namespace SMARTS291
}  // namespace Descriptors
}  // namespace RDKit
