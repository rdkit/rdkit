//
// Created by gareth on 11/4/20.
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

#include <boost/shared_ptr.hpp>

using namespace std;
using namespace RDKit;

// Example system based on Brian's MultipleCores notebook for profiling
int main() {
  RDLog::InitLogs();

  vector<shared_ptr<ROMol>> molecules;
  {
    fstream fh;
    std::string rdBase(getenv("RDBASE"));
    fh.open(rdBase + "/Docs/Notebooks/compounds.txt", ios::in);
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

  boost::shared_ptr<ROMol> cephem(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CS2"));
  boost::shared_ptr<ROMol> carbacephem(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CC2"));
  boost::shared_ptr<ROMol> oxacephem(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CO2"));
  boost::shared_ptr<ROMol> carbapenem(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])C2"));
  boost::shared_ptr<ROMol> carbapenam(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])C2"));
  boost::shared_ptr<ROMol> penem(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])S2"));
  boost::shared_ptr<ROMol> penam(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])S2"));
  boost::shared_ptr<ROMol> oxapenam(
      SmilesToMol("O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])O2"));
  boost::shared_ptr<ROMol> monobactam(SmilesToMol("O=C1C([*:1])C([*:5])N1"));

  vector<boost::shared_ptr<ROMol>> cores{cephem,     carbacephem, oxacephem,
                                         carbapenem, carbapenam,  penem,
                                         penam,      oxapenam,    monobactam};
  RGroupDecompositionParameters parameters;
  parameters.scoreMethod = FingerprintVariance;
  parameters.matchingStrategy = GA;
  RGroupDecomposition decomposition(cores, parameters);

  int numberAdded(0);
  for (auto &molecule : molecules) {
    numberAdded = decomposition.add(*molecule);
  }
  cerr << "Added " << numberAdded << " compounds to decomposition" << endl;
  auto result = decomposition.process();
  cerr << "Results success " << result.success << " score " << result.score;
}