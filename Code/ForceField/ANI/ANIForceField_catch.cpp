#include "RDGeneral/test.h"
#include "catch.hpp"

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <ForceField/ANI/AtomicContrib.h>
#include <ForceField/ForceField.h>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <Eigen/Dense>

using namespace Eigen;

TEST_CASE("ANI-1ccx NN Forward Pass", "[ANI Force Field]") {
  SECTION("CH4") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-40.0553)) < 0.01);
  }
  SECTION("SO2") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-150.1208)) < 0.01);
  }
  SECTION("NH3") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/NH3.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-56.1652)) < 0.01);
  }
  SECTION("Ethanol") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethanol.sdf";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-154.8904)) < 0.01);
  }
}

TEST_CASE("ANI-1x NN Forward Pass", "[ANI Force Field]") {
  SECTION("CH4") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-40.0517)) < 0.01);
  }
  SECTION("SO2") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-150.0904)) < 0.01);
  }
  SECTION("NH3") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/NH3.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-56.1592)) < 0.01);
  }
  SECTION("Ethanol") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethanol.sdf";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);
    auto species_vec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    ForceFields::ForceField field;
    auto &ps = field.positions();
    double *pos;
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      ps.push_back(&atom);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;

      ForceFields::ANI::ANIAtomContrib *ac;
      ac = new ForceFields::ANI::ANIAtomContrib(
          &field, species_vec(i), i, species_vec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    double total = 0;
    for (auto i : field.contribs()) {
      total += i->getEnergy(pos);
    }
    CHECK(std::fabs(total - (-154.9892)) < 0.01);
  }
}