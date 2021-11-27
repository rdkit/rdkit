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
#include <GraphMol/SmilesParse/SmilesWrite.h>
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
  desc.add_options()("help", "Help message")(
      "dataset", options::value<std::string>()->default_value("rg-easy"),
      "Built-in dataset rg-easy, rg-stereo")(
      "maximumOperations", options::value<int>()->default_value(-1),
      "Maximum number of GA operations")(
      "populationSize", options::value<int>()->default_value(-1),
      "GA population size")(
      "numberOperationsWithoutImprovement",
      options::value<int>()->default_value(-1),
      "number of operations without improvement before exiting")(
      "randomSeed", options::value<int>()->default_value(-1),
      "Random number seed (-1 for default, -2 for random)")(
      "matchingStrategy", options::value<std::string>()->default_value("GA"),
      "Matching strategy- GA or GreedyChunks")(
      "numberRuns", options::value<int>()->default_value(1),
      "Number of GA runs")("start", options::value<int>()->default_value(0),
                           "start")(
      "batch", options::value<int>()->default_value(10000000), "batch size");
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
  } else if (dataset == "beacon") {
    cerr << "Using Beacon dataset" << endl;
    file = rdBase + "/Docs/Notebooks/compounds2.txt";
    coreSmiles = {"N1([*:1])CCN([*:2])CC1", "C1(O[*:1])CCC(O[*:2])CC1",
                  "C1([*:1])CCC([*:2])CC1"};
  } else {
    cerr << "unknown dataset " << dataset;
    return 1;
  }

  vector<shared_ptr<ROMol>> molecules;
  {
    fstream fh;
    fh.open(file, ios::in);
    string line;
    getline(fh, line);

    int count = 0;
    int start = vm["start"].as<int>();
    int batch = vm["batch"].as<int>();
    while (getline(fh, line)) {
      count++;
      if (count < start) {
        continue;
      }
      if (count > start + batch) {
        break;
      }
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
  parameters.gaMaximumOperations = vm["maximumOperations"].as<int>();
  parameters.gaNumberOperationsWithoutImprovement =
      vm["numberOperationsWithoutImprovement"].as<int>();
  parameters.gaPopulationSize = vm["populationSize"].as<int>();
  parameters.gaRandomSeed = vm["randomSeed"].as<int>();
  parameters.gaNumberRuns = vm["numberRuns"].as<int>();
  auto strategyString = vm["matchingStrategy"].as<string>();
  if (strategyString == "GA") {
    parameters.matchingStrategy = GA;
  } else if (strategyString == "GreedyChunks") {
    parameters.matchingStrategy = GreedyChunks;
  } else {
    cerr << "Unknown matching strategy " << strategyString << endl;
    return 0;
  }
  RGroupDecomposition decomposition(cores, parameters);

  int numberAdded(0);
  for (auto& molecule : molecules) {
    auto added = decomposition.add(*molecule);
    if (added > -1) {
      auto smiles = MolToSmiles(*molecule);
      cerr << "\"" << smiles << "\", " << endl;
    }
    numberAdded = std::max(numberAdded, added);
  }
  cerr << "Added " << numberAdded << " compounds to decomposition" << endl;
  auto result = decomposition.processAndScore();
  cerr << "Results success " << result.success << " score " << result.score
       << endl;
}