//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/point.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/MolInteractionFields/MIFDescriptors.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace RDGeom;
using namespace RDKit;
using namespace RDMIF;

namespace RDMIF {

class testfunctor {
 public:
  testfunctor() {}
  double operator()(const double &, const double &, const double &,
                    double) const {
    return 1.0;
  }
};
}  // namespace RDMIF

TEST_CASE("constructGrid") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";
  //	//Generate Molecule with 3D Coordinates and Partial Charges
  //	RWMol mol=*SmilesToMol("Cl");
  //	MolOps::addHs(*mol);
  //	DGeomHelpers::EmbedMolecule(*mol);

  //	MolToMolFile(*mol, path + "HCl.mol");

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  REQUIRE(mol);
  auto grd = *constructGrid(*mol, 0, 5.0, 0.5);

  Point3D bond =
      mol->getConformer().getAtomPos(1) - mol->getConformer().getAtomPos(0);
  CHECK(feq(grd.getSpacing(), 0.5));
  CHECK(feq(grd.getNumX(), std::floor((fabs(bond.x) + 10.0) / 0.5 + 0.5)));
  CHECK(feq(grd.getNumY(), std::floor((fabs(bond.y) + 10.0) / 0.5 + 0.5)));
  CHECK(feq(grd.getNumZ(), std::floor((fabs(bond.z) + 10.0) / 0.5 + 0.5)));
}

TEST_CASE("CalculateDescriptors") {
  RealValueVect data(125, 0.0);
  Point3D o(0.0, 0.0, 0.0);
  UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 1.0, &o, &data);

  calculateDescriptors(grd, testfunctor());
  CHECK(grd.getOccupancyVect()->getTotalVal() == grd.getSize());
}

TEST_CASE("CubeFiles") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  REQUIRE(mol);
  RealValueVect data(125, 0.0);
  Point3D o(0.0, 0.0, 0.0);
  UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 1.0, &o, &data);
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, double(i / 10.0));
  }
  writeToCubeFile(grd, path + "test3.cube", mol.get());
  UniformRealValueGrid3D grd2;
  auto mol2 = readFromCubeFile(grd2, path + "test3.cube");

  CHECK(grd.getSize() == grd2.getSize());
  for (unsigned int i = 0; i < grd2.getSize(); i++) {
    CHECK(feq(grd2.getVal(i), double(i / 10.0)));
  }
  CHECK(grd.getNumX() == grd2.getNumX());
  CHECK(grd.getNumY() == grd2.getNumY());
  CHECK(grd.getNumZ() == grd2.getNumZ());
  CHECK(feq(grd.getOffset().x, grd2.getOffset().x));
  CHECK(feq(grd.getSpacing(), grd2.getSpacing()));
  CHECK(grd.compareVectors(grd2));
  CHECK(grd.compareParams(grd2));
  CHECK(grd.compareGrids(grd2));

  CHECK(mol->getNumAtoms() == mol2->getNumAtoms());
  for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
    CHECK(mol->getAtomWithIdx(i)->getAtomicNum() ==
          mol2->getAtomWithIdx(i)->getAtomicNum());
    CHECK(feq(mol->getConformer().getAtomPos(i).x,
              mol2->getConformer().getAtomPos(i).x));
    CHECK(feq(mol->getConformer().getAtomPos(i).y,
              mol2->getConformer().getAtomPos(i).y));
    CHECK(feq(mol->getConformer().getAtomPos(i).z,
              mol2->getConformer().getAtomPos(i).z));
  }
}

TEST_CASE("Coulomb") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  REQUIRE(mol);

  computeGasteigerCharges(*mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol->getConformer(0);
  for (auto i = 0u; i < mol->getNumAtoms(); ++i) {
    charges.push_back(
        mol->getAtomWithIdx(i)->getProp<double>("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  auto grd = *constructGrid(*mol);
  auto grd2 = *constructGrid(*mol);

  Coulomb coul(*mol);

  calculateDescriptors<Coulomb>(grd, coul);
  calculateDescriptors<Coulomb>(grd2, Coulomb(charges, pos));

  CHECK(grd.compareGrids(grd2));
  CHECK(feq(coul(0.0, 0.0, 0.0, 1000), 0.0));
  CHECK(coul(2.0, 0.0, 0.0, 1000) < 0);
  CHECK(coul(-2.0, 0.0, 0.0, 1000) > 0);
  CHECK(feq(coul(0.0, 0.0, 0.0, 0.1), 0.0));

  calculateDescriptors<Coulomb>(grd, Coulomb(*mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK(grd.getVal(i) <= 0.0);
  }

  Coulomb coul1(*mol, 0, -1.0, false, "_GasteigerCharge", 0.0, 0.01);
  CHECK(coul1(-2.0, 0.0, 0.0, 1000) < 0);
  CHECK(coul1(2.0, 0.0, 0.0, 1000) > 0);

  Coulomb coul2 = Coulomb(*mol, 0, -.5, false, "_GasteigerCharge", 0.0, 0.01);
  CHECK(coul1(-2.0, 0.0, 0.0, 1000) < coul2(-2.0, 0.0, 0.0, 1000));

  Coulomb coul3(*mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0);
  CHECK(coul3(0.0, 0.0, 0.0, 1000) > coul3(0.1, 0.0, 0.0, 1000));
  CHECK(coul3(0.66, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000));
  CHECK(coul3(0.70, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000));
  CHECK(feq(coul3(0.0, 0.0, 0.0, 0.1), 0.0));
}

TEST_CASE("CoulombDielectric") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  REQUIRE(mol);

  computeGasteigerCharges(*mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol->getConformer(0);
  for (auto i = 0u; i < mol->getNumAtoms(); ++i) {
    charges.push_back(
        mol->getAtomWithIdx(i)->getProp<double>("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  auto grd = *constructGrid(*mol);
  auto grd2 = *constructGrid(*mol);

  CoulombDielectric couldiele(*mol);

  calculateDescriptors<CoulombDielectric>(grd, couldiele);
  calculateDescriptors<CoulombDielectric>(grd2,
                                          CoulombDielectric(charges, pos));

  CHECK(grd.compareGrids(grd2));
  CHECK(feq(couldiele(0.0, 0.0, 0.0, 1000), 0.0));
  CHECK(couldiele(2.0, 0.0, 0.0, 1000) < 0);
  CHECK(couldiele(-2.0, 0.0, 0.0, 1000) > 0);

  calculateDescriptors<CoulombDielectric>(
      grd, CoulombDielectric(*mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK(grd.getVal(i) <= 0.0);
  }

  CHECK(feq(couldiele(0.0, 0.0, 0.0, 1000), 0.0));
  CHECK(couldiele(2.0, 0.0, 0.0, 1000) < 0);
  CHECK(couldiele(-2.0, 0.0, 0.0, 1000) > 0);

  calculateDescriptors<CoulombDielectric>(
      grd, CoulombDielectric(*mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK(grd.getVal(i) <= 0.0);
  }

  CoulombDielectric couldiele1(*mol, 0, -1.0);
  CHECK(couldiele1(-2.0, 0.0, 0.0, 1000) < 0);
  CHECK(couldiele1(2.0, 0.0, 0.0, 1000) > 0);

  CoulombDielectric couldiele2(*mol, 0, -.5);

  CHECK(couldiele1(-2.0, 0.0, 0.0, 1000) < couldiele2(-2.0, 0.0, 0.0, 1000));

  CoulombDielectric couldiele3(*mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0);
  CHECK(couldiele3(0.0, 0.0, 0.0, 1000) > couldiele3(0.1, 0.0, 0.0, 1000));
  CHECK(couldiele3(0.66, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000));
  CHECK(couldiele3(0.70, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000));

  mol = v2::FileParsers::MolFromMolFile(path + "glucose.mol", fopts);
  REQUIRE(mol);
  computeGasteigerCharges(*mol);

  CoulombDielectric couldiele4(*mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 80.0, 4.0);
  CoulombDielectric couldiele5(*mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 200.0, 4.0);
  CoulombDielectric couldiele6(*mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 80.0, 10.0);
  CHECK(couldiele5(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000));

  CHECK(couldiele6(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000));
}

TEST_CASE("VdWaals") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  REQUIRE(mol);
  try {
    MMFFVdWaals vdw(*mol, 0, 6, false, 1.0);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.what() << "\n";
  }

  mol = v2::FileParsers::MolFromMolFile(path + "HCN.mol", fopts);
  REQUIRE(mol);
  {
    MMFFVdWaals vdw(*mol, 0, 6, false, 1.0);

    CHECK(vdw(-5.0, 0, 0, 1000) < 0);
    CHECK(vdw(-1.68, 0, 0, 1000) > vdw(-5.0, 0, 0, 1000));
    CHECK(vdw(-5.0, 0, 0, 1000) < vdw(-10.0, 0, 0, 1000));
  }

  auto mol2 = v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);
  REQUIRE(mol2);
  MMFFVdWaals vdw(*mol2, 0, 6, false, 1.0);
  MMFFVdWaals vdw2(*mol2, 0, 6, true, 1.0);
  CHECK(fabs(vdw2(-3.0, 0, 0, 1000) - vdw(-3.0, 0, 0, 1000)) > 0.0001);

  UFFVdWaals vdw3(*mol, 0, "O_3", 1.0);
  CHECK(vdw3(-5.0, 0, 0, 1000) < 0);
  CHECK(vdw3(-1.68, 0, 0, 1000) > vdw3(-5.0, 0, 0, 1000));
  CHECK(vdw3(-5.0, 0, 0, 1000) < vdw3(-10.0, 0, 0, 1000));

  std::vector<std::string> names{
      "acetone",       "aceticacid",    "phenol",       "phenolate",
      "serine",        "threonine",     "ethanol",      "diethylether",
      "h2o",           "ammonia",       "ethylamine",   "imine",
      "acetonitrile",  "histidine",     "phenylamine",  "methylammonium",
      "fluoromethane", "chloromethane", "bromomethane", "glycine",
      "glyphe",        "glysergly",     "glythrgly",    "glucose"};
  for (const auto &name : names) {
    mol = v2::FileParsers::MolFromMolFile(path + name + ".mol", fopts);
    REQUIRE(mol);
    MMFFVdWaals mmffVdw(*mol);
    CHECK(mmffVdw(0.0, 0.0, 0.0, 1000));
    UFFVdWaals uffVdw(*mol);
    CHECK(uffVdw(0.0, 0.0, 0.0, 1000));
  }
}

TEST_CASE("HBond") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  // Generate Molecule with 3D Coordinates and Partial Charges
  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol =
      v2::FileParsers::MolFromMolFile(path + "ethane.mol", fopts);  // Ethane
  REQUIRE(mol);
  auto grd = *constructGrid(*mol, 0, 5.0, 1);
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, 1.0);
  }

  UniformRealValueGrid3D grd1(grd);

  HBond hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  TEST_ASSERT(!grd.compareGrids(grd1));
  const RealValueVect *vect = grd.getOccupancyVect();
  const RealValueVect *vect1 = grd1.getOccupancyVect();
  CHECK((unsigned int)fabs(((*vect1 - *vect).getTotalVal())) == grd.getSize());

  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);

  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, 1.0);
  }
  calculateDescriptors<HBond>(grd, hbonddes);
  TEST_ASSERT(!grd.compareGrids(grd1));
  CHECK((unsigned int)fabs(((*vect1 - *vect).getTotalVal())));

  mol = v2::FileParsers::MolFromMolFile(path + "aceticacid.mol",
                                        fopts);  // Acetic Acid
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 2);
  calculateDescriptors<HBond>(grd, hbonddes);
  HBond hbonddes1(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes1.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes1);
  HBond hbonddes2(*mol, 0, "O", false, 0.001);
  CHECK(hbonddes1(4.0, 0.0, 1.0, 1000) > hbonddes2(4.0, 0.0, 1.0, 1000));
  HBond hbonddes3(*mol, 0, "NH", true, 0.001);
  CHECK(hbonddes3.getNumInteractions() == 2);
  CHECK(hbonddes(2.0, 2.0, 2.0, 1000) < hbonddes3(2.0, 2.0, 2.0, 1000));
  HBond hbonddes4(*mol, 0, "N", true, 0.001);
  CHECK(hbonddes4.getNumInteractions() == 1);
  CHECK(hbonddes1(3.0, 0.0, 0.0, 1000) < hbonddes4(3.0, 0.0, 0.0, 1000));

  mol =
      v2::FileParsers::MolFromMolFile(path + "acetone.mol", fopts);  // Acetone
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "diethylether.mol",
                                        fopts);  // Et2O
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);  // H2O
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 2);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol =
      v2::FileParsers::MolFromMolFile(path + "ammonia.mol", fopts);  // ammonia
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 3);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "imine.mol", fopts);  // imine
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "methylammonium.mol",
                                        fopts);  // methylammonium
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 3);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "chloromethane.mol",
                                        fopts);  // Chloromethane
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "phosphonate.mol",
                                        fopts);  // Phosphonate
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 3);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "phosphatediester.mol",
                                        fopts);  // Phosphatediester
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 4);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "hydrogenphosphatediester.mol",
                                        fopts);  // Hydrogenphosphatediester
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 4);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "mustardgas.mol",
                                        fopts);  // mustard gas
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 2);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "alicin.mol", fopts);  // Alicin
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 1);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 0);
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = v2::FileParsers::MolFromMolFile(path + "sulfanilamide.mol",
                                        fopts);  // Sulfanilamide
  REQUIRE(mol);
  grd = *constructGrid(*mol, 0, 5.0, 1);
  hbonddes = HBond(*mol, 0, "OH", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 3);
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(*mol, 0, "O", true, 0.001);
  CHECK(hbonddes.getNumInteractions() == 4);
  calculateDescriptors<HBond>(grd, hbonddes);
}

TEST_CASE("Hydrophilic") {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);
  REQUIRE(mol);

  Hydrophilic hydro(*mol);
  HBond hbondOH(*mol, 0, "OH");
  HBond hbondO(*mol, 0, "O");

  double hyd = hydro(0.0, 0.0, 0.0, 1000), hOH = hbondOH(0.0, 0.0, 0.0, 1000),
         hO = hbondO(0.0, 0.0, 0.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(1.0, 1.5, 2.0, 1000);
  hOH = hbondOH(1.0, 1.5, 2.0, 1000);
  hO = hbondO(1.0, 1.5, 2.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(2.0, 1.5, -3.0, 1000);
  hOH = hbondOH(2.0, 1.5, -3.0, 1000);
  hO = hbondO(2.0, 1.5, -3.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(-2.5, 0.5, 3.0, 1000);
  hOH = hbondOH(-2.5, 0.5, 3.0, 1000);
  hO = hbondO(-2.5, 0.5, 3.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(10.0, 1.5, 1.0, 1000);
  hOH = hbondOH(10.0, 1.5, 1.0, 1000);
  hO = hbondO(10.0, 1.5, 1.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(6.0, -5.0, 0.0, 1000);
  hOH = hbondOH(6.0, -5.0, 0.0, 1000);
  hO = hbondO(6.0, -5.0, 0.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(-3.0, -3.0, 7.0, 1000);
  hOH = hbondOH(-3.0, -3.0, 7.0, 1000);
  hO = hbondO(-3.0, -3.0, 7.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(1.0, 0.0, 0.0, 1000);
  hOH = hbondOH(1.0, 0.0, 0.0, 1000);
  hO = hbondO(1.0, 0.0, 0.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(0.0, 2.0, 2.0, 1000);
  hOH = hbondOH(0.0, 2.0, 2.0, 1000);
  hO = hbondO(0.0, 2.0, 2.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(2.0, -2.0, 0.0, 1000);
  hOH = hbondOH(2.0, -2.0, 0.0, 1000);
  hO = hbondO(2.0, -2.0, 0.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));

  hyd = hydro(2.0, -2.0, -3.0, 1000);
  hOH = hbondOH(2.0, -2.0, -3.0, 1000);
  hO = hbondO(2.0, -2.0, -3.0, 1000);
  CHECK(feq(std::min(hOH, hO), hyd));
}
