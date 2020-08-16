//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <ForceField/ANI/AtomicContrib.h>
#include <ForceField/ForceField.h>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/DistanceConstraint.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <ForceField/MMFF/PositionConstraint.h>
#include <Eigen/Dense>
#include <chrono>

using namespace Eigen;
using namespace Catch::literals;

TEST_CASE("ANI-1ccx NN Forward Pass", "[ANI Force Field]") {
  std::string pathName = getenv("RDBASE");
  std::string dirPath = pathName + "/Code/GraphMol/Descriptors/test_data/";
  SECTION("CH4") {
    std::string molFile = dirPath + "CH4.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-40.0553_a).margin(0.05));
  }
  SECTION("NH3") {
    std::string molFile = dirPath + "NH3.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-56.1652_a).margin(0.05));
  }
  SECTION("Ethanol") {
    std::string molFile = dirPath + "ethanol.sdf";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1ccx");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-154.8904_a).margin(0.05));
  }
}

TEST_CASE("ANI-1x NN Forward Pass", "[ANI Force Field]") {
  std::string pathName = getenv("RDBASE");
  std::string dirPath = pathName + "/Code/GraphMol/Descriptors/test_data/";
  SECTION("CH4") {
    std::string molFile = dirPath + "CH4.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-40.0517_a).margin(0.05));
  }
  SECTION("NH3") {
    std::string molFile = dirPath + "NH3.mol";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-56.1592_a).margin(0.05));
  }
  SECTION("Ethanol") {
    std::string molFile = dirPath + "ethanol.sdf";

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;
    auto conf = mol->getConformer(confId);
    auto numAtoms = mol->getNumAtoms();
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
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

      auto ac = new ForceFields::ANI::ANIAtomContrib(
          &field, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      field.contribs().push_back(ForceFields::ContribPtr(ac));
    }
    field.initialize();
    CHECK(field.calcEnergy(pos) == (-154.9892_a).margin(0.05));
  }
}
TEST_CASE("QM7 Test Cases") {
  SECTION("3C 1N and 5H Molecule") {
    int species[] = {6, 6, 6, 6, 7, 1, 1, 1, 1, 1};
    unsigned int numAtoms = 10;
    auto speciesVec =
        RDKit::Descriptors::ANI::GenerateSpeciesVector(species, numAtoms);
    ForceFields::ForceField fieldANI1x, fieldANI1ccx;

    double pos[] = {
        1.8633455038070679, 0.02800574153661728,     -0.012358808889985085,
        4.697613716125488,  -0.00013228082389105111, -0.005952637176960707,
        5.985896587371826,  -2.188642978668213,      0.009750986471772194,
        6.029568195343018,  2.347965717315674,       -0.011772993952035904,
        7.106976509094238,  4.2567973136901855,      -0.015174500644207,
        1.120683193206787,  1.9597971439361572,      0.04142279550433159,
        1.1187934875488281, -0.9710546731948853,     1.640320062637329,
        1.1301696300506592, -0.8761904239654541,     -1.7236380577087402,
        5.016429424285889,  -3.997413396835327,      0.014909938909113407,
        8.040104866027832,  -2.2334485054016113,     0.015911493450403214};
    for (unsigned int i = 0; i < numAtoms; i++) {
      ForceFields::ANI::ANIAtomContrib *ac1;
      ac1 = new ForceFields::ANI::ANIAtomContrib(
          &fieldANI1x, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      ForceFields::ANI::ANIAtomContrib *ac2;
      ac2 = new ForceFields::ANI::ANIAtomContrib(&fieldANI1ccx, speciesVec(i),
                                                 i, speciesVec, numAtoms, 4, 8,
                                                 "ANI-1ccx");
      fieldANI1x.contribs().push_back(ForceFields::ContribPtr(ac1));
      fieldANI1ccx.contribs().push_back(ForceFields::ContribPtr(ac2));
      auto position =
          RDGeom::Point3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      fieldANI1ccx.positions().push_back(&position);
      fieldANI1x.positions().push_back(&position);
    }
    fieldANI1x.initialize();
    fieldANI1ccx.initialize();
    CHECK(fieldANI1x.calcEnergy(pos) == (-209.0046_a).margin(0.05));
    CHECK(fieldANI1ccx.calcEnergy(pos) == (-208.9733_a).margin(0.05));
  }
  SECTION("1N 3C 1O and 9H") {
    int species[] = {7, 6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    unsigned int numAtoms = 14;
    auto speciesVec =
        RDKit::Descriptors::ANI::GenerateSpeciesVector(species, numAtoms);
    ForceFields::ForceField fieldANI1x, fieldANI1ccx;

    double pos[] = {
        2.2074835300445557, -0.3340090811252594,  0.10053343325853348,
        4.9639892578125,    -0.22546322643756866, 0.1855711042881012,
        5.859700679779053,  2.5101420879364014,   -0.05463198199868202,
        8.73183822631836,   2.706106662750244,    -0.03679296746850014,
        9.749720573425293,  1.7906666994094849,   2.2751736640930176,
        1.6450822353363037, -2.170199394226074,   0.2544138431549072,
        1.5863306522369385, 0.2723284363746643,   -1.618361473083496,
        5.743180274963379,  -1.3758339881896973,  -1.3486785888671875,
        5.6071953773498535, -1.0400863885879517,  1.9756519794464111,
        5.088220119476318,  3.6159720420837402,   1.5213806629180908,
        5.141377925872803,  3.341905117034912,    -1.8118693828582764,
        9.32322883605957,   4.67402982711792,     -0.26873794198036194,
        9.555381774902344,  1.6062861680984497,   -1.5811527967453003,
        9.00974178314209,   2.783963441848755,    3.6356818675994873};
    for (unsigned int i = 0; i < numAtoms; i++) {
      ForceFields::ANI::ANIAtomContrib *ac1;
      ac1 = new ForceFields::ANI::ANIAtomContrib(
          &fieldANI1x, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      ForceFields::ANI::ANIAtomContrib *ac2;
      ac2 = new ForceFields::ANI::ANIAtomContrib(&fieldANI1ccx, speciesVec(i),
                                                 i, speciesVec, numAtoms, 4, 8,
                                                 "ANI-1ccx");
      fieldANI1x.contribs().push_back(ForceFields::ContribPtr(ac1));
      fieldANI1ccx.contribs().push_back(ForceFields::ContribPtr(ac2));
      auto position =
          RDGeom::Point3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      fieldANI1ccx.positions().push_back(&position);
      fieldANI1x.positions().push_back(&position);
    }
    fieldANI1x.initialize();
    fieldANI1ccx.initialize();
    CHECK(fieldANI1x.calcEnergy(pos) == (-248.1921_a).margin(0.05));
    CHECK(fieldANI1ccx.calcEnergy(pos) == (-248.2279_a).margin(0.05));
  }
  SECTION("1N 5C 9H") {
    int species[] = {7, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    unsigned int numAtoms = 15;
    auto speciesVec =
        RDKit::Descriptors::ANI::GenerateSpeciesVector(species, numAtoms);
    ForceFields::ForceField fieldANI1x, fieldANI1ccx;

    double pos[] = {
        2.1612987518310547, -0.6327747702598572,  0.43304964900016785,
        4.886623859405518,  -0.5855505466461182,  0.8615450263023376,
        6.046708106994629,  1.7641916275024414,   -0.3568369746208191,
        8.887079238891602,  1.9741591215133667,   0.131562739610672,
        10.33105754852295,  0.006821911316365004, -1.1970659494400024,
        11.513798713684082, -1.6008626222610474,  -2.2824490070343018,
        1.4307305812835693, -2.2506449222564697,  1.1816647052764893,
        1.8104143142700195, -0.7166597247123718,  -1.4598512649536133,
        5.720730304718018,  -2.310190200805664,   0.07972754538059235,
        5.243819713592529,  -0.5787853002548218,  2.9002761840820312,
        5.127204895019531,  3.4552319049835205,   0.416646808385849,
        5.6790995597839355, 1.7595807313919067,   -2.398723840713501,
        9.267745971679688,  1.8809388875961304,   2.165493965148926,
        9.537825584411621,  3.8325536251068115,   -0.5140811204910278,
        12.562162399291992, -3.02358078956604,    -3.24534010887146};
    for (unsigned int i = 0; i < numAtoms; i++) {
      ForceFields::ANI::ANIAtomContrib *ac1;
      ac1 = new ForceFields::ANI::ANIAtomContrib(
          &fieldANI1x, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      ForceFields::ANI::ANIAtomContrib *ac2;
      ac2 = new ForceFields::ANI::ANIAtomContrib(&fieldANI1ccx, speciesVec(i),
                                                 i, speciesVec, numAtoms, 4, 8,
                                                 "ANI-1ccx");
      fieldANI1x.contribs().push_back(ForceFields::ContribPtr(ac1));
      fieldANI1ccx.contribs().push_back(ForceFields::ContribPtr(ac2));
      auto position =
          RDGeom::Point3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      fieldANI1ccx.positions().push_back(&position);
      fieldANI1x.positions().push_back(&position);
    }
    fieldANI1x.initialize();
    fieldANI1ccx.initialize();
    CHECK(fieldANI1x.calcEnergy(pos) == (-249.0084_a).margin(0.05));
    CHECK(fieldANI1ccx.calcEnergy(pos) == (-249.0295_a).margin(0.05));
  }
  SECTION("5C 1O 8H") {
    int species[] = {6, 6, 6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1};
    unsigned int numAtoms = 14;
    auto speciesVec =
        RDKit::Descriptors::ANI::GenerateSpeciesVector(species, numAtoms);
    ForceFields::ForceField fieldANI1x, fieldANI1ccx;

    double pos[] = {
        1.9601751565933228, -0.061945222318172455, -0.08331802487373352,
        4.777057647705078,  -0.06096256524324417,  -0.04117713123559952,
        6.238533973693848,  -1.251603364944458,    -2.1706528663635254,
        6.237986087799072,  1.5869730710983276,    -1.860718846321106,
        8.574311256408691,  2.8852715492248535,    -1.034757375717163,
        10.139514923095703, 3.6033108234405518,    -2.578569173812866,
        1.2324415445327759, 1.4809783697128296,    1.0859311819076538,
        1.2382619380950928, -1.844051480293274,    0.6785250902175903,
        1.2088199853897095, 0.1712280809879303,    -1.997175931930542,
        5.501691818237305,  -0.3164535462856293,   1.858847975730896,
        5.213395118713379,  -2.030416250228882,    -3.7652416229248047,
        7.952043056488037,  -2.279557704925537,    -1.712791085243225,
        5.169175624847412,  2.650114059448242,     -3.2528233528137207,
        8.815553665161133,  3.1873443126678467,    1.0119673013687134,
    };
    for (unsigned int i = 0; i < numAtoms; i++) {
      ForceFields::ANI::ANIAtomContrib *ac1;
      ac1 = new ForceFields::ANI::ANIAtomContrib(
          &fieldANI1x, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      ForceFields::ANI::ANIAtomContrib *ac2;
      ac2 = new ForceFields::ANI::ANIAtomContrib(&fieldANI1ccx, speciesVec(i),
                                                 i, speciesVec, numAtoms, 4, 8,
                                                 "ANI-1ccx");
      fieldANI1x.contribs().push_back(ForceFields::ContribPtr(ac1));
      fieldANI1ccx.contribs().push_back(ForceFields::ContribPtr(ac2));
      auto position =
          RDGeom::Point3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      fieldANI1ccx.positions().push_back(&position);
      fieldANI1x.positions().push_back(&position);
    }
    fieldANI1x.initialize();
    fieldANI1ccx.initialize();
    CHECK(fieldANI1x.calcEnergy(pos) == (-269.0078_a).margin(0.05));
    CHECK(fieldANI1ccx.calcEnergy(pos) == (-269.0295_a).margin(0.05));
  }
  SECTION("6C 1N 13H") {
    int species[] = {6, 6, 6, 7, 6, 6, 6, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    unsigned int numAtoms = 20;
    auto speciesVec =
        RDKit::Descriptors::ANI::GenerateSpeciesVector(species, numAtoms);
    ForceFields::ForceField fieldANI1x, fieldANI1ccx;

    double pos[] = {
        1.6804956197738647, -0.2058289647102356,  0.3033577501773834,
        4.54233455657959,   -0.02983877621591091, -0.05497213453054428,
        5.335660457611084,  -0.09013993293046951, -2.8519558906555176,
        4.461907863616943,  2.1463699340820312,   -4.214334964752197,
        5.836381435394287,  2.0427372455596924,   1.597517728805542,
        7.893631935119629,  0.09807678312063217,  2.2748334407806396,
        6.003244400024414,  -2.0027129650115967,  1.5839306116104126,
        1.1721593141555786, -0.3195337951183319,  2.3073368072509766,
        0.9210147261619568, -1.8902552127838135,  -0.6298835277557373,
        0.7200612425804138, 1.4526702165603638,   -0.47447243332862854,
        4.586913585662842,  -1.7891360521316528,  -3.771118640899658,
        7.397408485412598,  -0.17650042474269867, -3.0236940383911133,
        2.5398108959198,    2.1276614665985107,   -4.326868057250977,
        5.0936431884765625, 2.0650928020477295,   -6.032308101654053,
        4.723445892333984,  2.6494150161743164,   3.2372143268585205,
        6.513205528259277,  3.728240728378296,    0.6133106350898743,
        9.554058074951172,  0.17268317937850952,  1.0386879444122314,
        8.529033660888672,  0.1181078851222992,   4.241074562072754,
        6.8076629638671875, -3.61893892288208,    0.5778782367706299,
        4.953671455383301,  -2.7142891883850098,  3.2221720218658447,

    };
    for (unsigned int i = 0; i < numAtoms; i++) {
      ForceFields::ANI::ANIAtomContrib *ac1;
      ac1 = new ForceFields::ANI::ANIAtomContrib(
          &fieldANI1x, speciesVec(i), i, speciesVec, numAtoms, 4, 8, "ANI-1x");
      ForceFields::ANI::ANIAtomContrib *ac2;
      ac2 = new ForceFields::ANI::ANIAtomContrib(&fieldANI1ccx, speciesVec(i),
                                                 i, speciesVec, numAtoms, 4, 8,
                                                 "ANI-1ccx");
      fieldANI1x.contribs().push_back(ForceFields::ContribPtr(ac1));
      fieldANI1ccx.contribs().push_back(ForceFields::ContribPtr(ac2));
      auto position =
          RDGeom::Point3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      fieldANI1ccx.positions().push_back(&position);
      fieldANI1x.positions().push_back(&position);
    }
    fieldANI1x.initialize();
    fieldANI1ccx.initialize();
    CHECK(fieldANI1x.calcEnergy(pos) == (-289.1283_a).margin(0.05));
    CHECK(fieldANI1ccx.calcEnergy(pos) == (-289.1044_a).margin(0.05));
  }
}
