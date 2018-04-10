//
// Writing molecules - example3.cpp

#include <iostream>
#include <string>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>

int main(int argc, char **argv) {
  std::string file_root = getenv("RDBASE");
  file_root += "/Docs/Book";
  std::string mol_file = file_root + "/data/chiral.mol";
  RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(mol_file));
  std::cout << RDKit::MolToSmiles(*mol, true) << std::endl;
  // 2nd parameter doIsomericSmiles defaults to true
  std::cout << RDKit::MolToSmiles(*mol, false) << std::endl;

  RDKit::ROMOL_SPTR mol1(RDKit::SmilesToMol("C1=CC=CN=C1"));
  std::cout << RDKit::MolToSmiles(*mol1) << std::endl;

  RDKit::ROMOL_SPTR mol2(RDKit::SmilesToMol("c1cccnc1"));
  std::cout << RDKit::MolToSmiles(*mol2) << std::endl;

  RDKit::ROMOL_SPTR mol3(RDKit::SmilesToMol("n1ccccc1"));
  std::cout << RDKit::MolToSmiles(*mol3) << std::endl;

  RDKit::RWMOL_SPTR mol4(new RDKit::RWMol(*mol));
  RDKit::MolOps::Kekulize(*mol4);
  std::cout << RDKit::MolToSmiles(*mol4) << std::endl;

  mol1.reset(RDKit::SmilesToMol("C1CCC1"));
  std::cout << RDKit::MolToMolBlock(*mol1) << std::endl;

  mol1->setProp("_Name", "cyclobutane");
  std::cout << RDKit::MolToMolBlock(*mol1) << std::endl;
}
