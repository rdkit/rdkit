#include <iostream>
#include <cmath>
#include <algorithm>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <CoordGen/CoordGen.h>

using namespace RDKit;

double getBondLength(const Conformer &conf, const Bond *bond) {
    auto pos1 = conf.getAtomPos(bond->getBeginAtomIdx());
    auto pos2 = conf.getAtomPos(bond->getEndAtomIdx());
    auto dx = pos1.x - pos2.x;
    auto dy = pos1.y - pos2.y;
    return std::sqrt(dx*dx + dy*dy);
}

void getBondStats(const ROMol &mol, double &median, double &maxLen) {
    if (mol.getNumConformers() == 0) {
        median = maxLen = 0.0;
        return;
    }

    const auto &conf = mol.getConformer();
    std::vector<double> lengths;

    for (const auto bond : mol.bonds()) {
        lengths.push_back(getBondLength(conf, bond));
    }

    std::sort(lengths.begin(), lengths.end());
    median = lengths[lengths.size() / 2];
    maxLen = *std::max_element(lengths.begin(), lengths.end());
}

int main() {
    // Test CXSMILES with elongated bond (21-carbon macrocycle)
    std::string cxsmiles = "C1CCCCCCCCCCCCCCCCCCCC1 |(-4.96999,3.68927,;-5.68003,2.4595,;-7.09998,2.45952,;-7.81002,1.22975,;-9.22997,1.22977,;-9.94001,0,;-9.23001,-1.22976,;-7.80997,-1.22974,;-7.09997,-2.45949,;-5.68002,-2.45952,;-4.26,-3.68927,;-2.84003,-2.45953,;-1.41999,-2.45951,;-0.710038,-1.22978,;0.71,-1.22976,;1.42004,0,;0.710038,1.22978,;-0.71,1.22976,;-1.42,2.45951,;-2.84004,2.45949,;-3.55004,3.68925,)|";

    std::cout << "Testing minimizeOnly with 21-carbon macrocycle\n\n";

    // Parse molecule with coordinates
    RWMol *mol = SmilesToMol(cxsmiles);
    if (!mol) {
        std::cerr << "Failed to parse CXSMILES\n";
        return 1;
    }

    // Get stats before minimization
    double median_before, max_before;
    getBondStats(*mol, median_before, max_before);

    std::cout << "BEFORE minimization:\n";
    std::cout << "  Median bond length: " << median_before << "\n";
    std::cout << "  Max bond length: " << max_before << "\n\n";

    // Set up coordgen parameters with minimizeOnly
    CoordGen::CoordGenParams params;
    params.minimizeOnly = true;
    params.minimizerPrecision = params.sketcherBestPrecision;

    std::cout << "Running coordgen minimization with minimizeOnly=true...\n";

    // Run minimization
    CoordGen::addCoords(*mol, &params);

    // Get stats after minimization
    double median_after, max_after;
    getBondStats(*mol, median_after, max_after);

    std::cout << "\nAFTER minimization:\n";
    std::cout << "  Median bond length: " << median_after << "\n";
    std::cout << "  Max bond length: " << max_after << "\n";
    std::cout << "  Improvement: " << (max_before - max_after) << "\n\n";

    // Output CXSMILES
    std::cout << "Output CXSMILES:\n";
    std::cout << MolToCXSmiles(*mol) << "\n";

    delete mol;
    return 0;
}
