//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>

#include <boost/program_options.hpp>

using namespace std;
using namespace RDKit;
using namespace boost::program_options;
namespace options = boost::program_options;

// Example systems based on Brian's MultipleCores notebook for profiling
int main(int argc, char* argv[]) {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.debug");

  options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "dataset", options::value<std::string>()->default_value("rg-easy"),
      "built-in dataset rg-easy, rg-stereo");
  options::variables_map vm;
  options::store(options::parse_command_line(argc, argv, desc), vm);
  options::notify(vm);

  if (vm.count("help")) {
    cerr << desc << endl;
    return 0;
  }

  auto dataset = vm["dataset"].as<std::string>();
  std::string rdBase(getenv("RDBASE"));
  std::string file;
  std::vector<std::string> coreSmiles;
  if (dataset == "rg-easy") {
    cerr << "Using dataset rg-easy" << endl;
    file = rdBase + "/Docs/Notebooks/compounds.txt";
    coreSmiles = {"O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CS2",
                  "O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CC2",
                  "O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CO2",
                  "O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])C2",
                  "O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])C2",
                  "O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])S2",
                  "O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])S2",
                  "O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])O2",
                  "O=C1C([*:1])C([*:5])N1"};
  } else if (dataset == "rg-stereo") {
    cerr << "Using dataset rg-stereo" << endl;
    file = rdBase + "/Docs/Notebooks/compounds.txt";
    coreSmiles = {"O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CS2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CC2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CO2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])C2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])C2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])S2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])S2",
                  "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])O2",
                  "O=C1C([*:2])([*:1])C([*:6])([*:5])N1"};
  }

  vector<shared_ptr<ROMol>> molecules;
  {
    fstream fh;
    fh.open(file, ios::in);
    string line;
    getline(fh, line);

    while (getline(fh, line)) {
      int pos = line.find_last_of("\t");
      auto smiles = line.substr(pos + 1);
      shared_ptr<ROMol> mol(SmilesToMol(smiles));
      molecules.push_back(mol);
    }
  }
  cerr << "Read " << molecules.size() << endl;

  std::vector<ROMOL_SPTR> cores;
  for (const auto& s : coreSmiles) {
    cores.emplace_back(SmartsToMol(s));
  }

  RGroupDecompositionParameters parameters;
  parameters.scoreMethod = FingerprintVariance;
  parameters.matchingStrategy = GA;
  RGroupDecomposition decomposition(cores, parameters);

  int numberAdded(0);
  for (auto& molecule : molecules) {
    numberAdded = decomposition.add(*molecule);
  }
  cerr << "Added " << numberAdded << " compounds to decomposition" << endl;
  auto result = decomposition.processAndScore();
  cerr << "Results success " << result.success << " score " << result.score
       << endl;
}