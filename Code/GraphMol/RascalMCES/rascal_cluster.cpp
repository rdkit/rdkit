//
// Created by David Cosgrove on 31/07/2023.
//
#include <cstdlib>
#include <iostream>

#include <RDGeneral/versions.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include "RascalMCES.h"
#include "RascalOptions.h"
#include "RascalResult.h"

bool read_file(std::string mol_file,
               std::vector<std::unique_ptr<RDKit::ROMol>> &mols) {
  if (mol_file.size() < 4) {
    std::cout << "Incorrect input filename " << mol_file << std::endl;
    return false;
  }
  std::string smi_suffix(".smi");
  std::string sdf_suffix(".sdf");
  std::unique_ptr<RDKit::MolSupplier> suppl;
  if (std::equal(begin(smi_suffix), end(smi_suffix),
                 end(mol_file) - smi_suffix.size())) {
    std::cout << "SMILES file : " << mol_file << std::endl;
    suppl.reset(
        new RDKit::SmilesMolSupplier(mol_file, " \t", 0, 1, false, true));
  } else if (std::equal(begin(sdf_suffix), end(sdf_suffix),
                        end(mol_file) - sdf_suffix.size())) {
    std::cout << "SDF file" << std::endl;
    suppl.reset(new RDKit::SDMolSupplier(mol_file));
  } else {
    std::cout << "Incorrect input filename " << mol_file << std::endl;
    return false;
  }
  while (!suppl->atEnd()) {
    mols.push_back(std::unique_ptr<RDKit::ROMol>(suppl->next()));
  }
  return true;
}

bool parse_args(int argc, char **argv, std::string &mol_file,
                RDKit::RascalMCES::RascalOptions &opts) {
  if (argc < 3) {
    std::cout
        << "ERROR : needs at least the command line options '--input-file'"
           " with a filename"
        << std::endl;
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg = std::string(argv[i]);
    if (arg == "--input-file") {
      if (i == argc - 1) {
        std::cout
            << "ERROR : --input-file requires 1 argument, the input file name."
            << std::endl;
        return false;
      } else {
        mol_file = std::string(argv[++i]);
      }
    } else if (arg == "--threshold") {
      if (i == argc - 1) {
        std::cout
            << "ERROR : --threshold requires 1 argument, the similarity threshold,"
               " a number between 0.0 and 1.0."
            << std::endl;
        return false;
      } else {
        opts.similarityThreshold = std::atof(argv[++i]);
      }
    } else if (arg == "--minimum-fragment-size") {
      if (i == argc - 1) {
        std::cout << "ERROR : --minimum-fragment-size requires 1 argument, an"
                     " integer greater than 1."
                  << std::endl;
        return false;
      } else {
        opts.minFragSize = std::atoi(argv[++i]);
      }
    } else if (arg == "--maximum-fragment-separation") {
      if (i == argc - 1) {
        std::cout
            << "ERROR : --maximum-fragment-separation requires 1 argument, an"
               " integer greater than 0."
            << std::endl;
        return false;
      } else {
        opts.maxFragSeparation = std::atoi(argv[++i]);
      }
    } else if (arg == "--timeout") {
      if (i == argc - 1) {
        std::cout << "ERROR : --timeout requires 1 argument, an"
                     " integer greater than 0."
                  << std::endl;
        return false;
      } else {
        opts.timeout = std::atoi(argv[++i]);
      }
    } else if (arg == "--complete-aromatic-rings") {
      opts.completeAromaticRings = true;
    } else if (arg == "--ring-matches-ring-only") {
      opts.ringMatchesRingOnly = true;
    } else if (arg == "--single-largest-frag") {
      opts.singleLargestFrag = true;
    } else if (arg == "--all-best-mces") {
      opts.allBestMCESs = true;
    } else if (arg == "--exact-chirality") {
      opts.exactChirality = true;
    } else {
      std::cout << "ERROR : unrecognised option " << arg << std::endl;
      return false;
    }
  }
  return true;
}

int main(int argc, char **argv) {
  std::cout << "Using RDKit version " << RDKit::rdkitVersion << std::endl;
  std::string mol_file;
  RDKit::RascalMCES::RascalOptions opts;
  if (!parse_args(argc, argv, mol_file, opts)) {
    exit(1);
  }

  std::vector<std::unique_ptr<RDKit::ROMol>> mols;
  if (!read_file(mol_file, mols)) {
    exit(1);
  }
  if (mols.size() < 2) {
    std::cout << "Need at least 2 molecules for MCES calculation." << std::endl;
  }
  std::cout << "Read " << mols.size() << " molecules." << std::endl;
  auto res = RDKit::RascalMCES::rascalMces(*mols[0], *mols[1], opts);
  std::cout << "Number of results : " << res.size() << std::endl;
#if 0
    for (const auto &r: res) {
        print_bond_matches(r, std::cout);
    }
#endif
}
