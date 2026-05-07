#include <iostream>
#include <fstream>
#include <string>
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

    if (lengths.empty()) {
        median = maxLen = 0.0;
        return;
    }

    std::sort(lengths.begin(), lengths.end());
    median = lengths[lengths.size() / 2];
    maxLen = *std::max_element(lengths.begin(), lengths.end());
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.smi output.smi [threshold_factor]\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    double threshold_factor = 1.1;

    if (argc > 3) {
        threshold_factor = std::stod(argv[3]);
    }

    std::ifstream infile(input_file);
    std::ofstream outfile(output_file);

    if (!infile.is_open()) {
        std::cerr << "Failed to open input file: " << input_file << "\n";
        return 1;
    }

    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: " << output_file << "\n";
        return 1;
    }

    int line_num = 0;
    int minimized_count = 0;
    int total_count = 0;
    std::string line;

    // Set up coordgen parameters
    CoordGen::CoordGenParams params;
    params.minimizeOnly = true;
    params.minimizerPrecision = params.sketcherBestPrecision;

    while (std::getline(infile, line)) {
        line_num++;

        if (line.empty()) {
            continue;
        }

        total_count++;

        // Parse CXSMILES
        RWMol *mol = SmilesToMol(line);
        if (!mol) {
            std::cerr << "Line " << line_num << ": Failed to parse\n";
            outfile << line << "\n";
            continue;
        }

        // Get bond stats
        double median, maxLen;
        getBondStats(*mol, median, maxLen);

        bool needs_min = (median > 0 && maxLen > median * threshold_factor);

        if (needs_min) {
            std::cout << "Line " << line_num << ": median=" << median
                     << ", max=" << maxLen << " - MINIMIZING\n";

            try {
                // Run minimization
                CoordGen::addCoords(*mol, &params);

                // Get new stats
                double new_median, new_max;
                getBondStats(*mol, new_median, new_max);

                std::cout << "         After: median=" << new_median
                         << ", max=" << new_max << "\n";

                // Output minimized CXSMILES
                outfile << MolToCXSmiles(*mol) << "\n";
                minimized_count++;

            } catch (const std::exception &e) {
                std::cerr << "Line " << line_num << ": Minimization failed - "
                         << e.what() << "\n";
                outfile << line << "\n";
            }
        } else {
            std::cout << "Line " << line_num << ": median=" << median
                     << ", max=" << maxLen << " - OK\n";
            outfile << line << "\n";
        }

        delete mol;
    }

    std::cout << "\nProcessed " << total_count << " structures\n";
    std::cout << "Minimized " << minimized_count << " structures\n";
    std::cout << "Output written to: " << output_file << "\n";

    return 0;
}
