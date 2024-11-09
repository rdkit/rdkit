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
#include <stdlib.h>

using namespace RDGeom;
using namespace RDKit;
using namespace RDMIF;

namespace RDMIF {

class testfunctor {
 public:
  testfunctor() {}
  double operator()(const double &, const double &, const double &, double) {
    return 1.0;
  }
};
}  // namespace RDMIF

void test1ConstructGrid() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";
  //	//Generate Molecule with 3D Coordinates and Partial Charges
  //	RWMol mol=*SmilesToMol("Cl");
  //	MolOps::addHs(mol);
  //	DGeomHelpers::EmbedMolecule(mol);

  //	MolToMolFile(mol, path + "HCl.mol");

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = *v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  auto grd = *constructGrid(mol, 0, 5.0, 0.5);

  Point3D bond =
      mol.getConformer().getAtomPos(1) - mol.getConformer().getAtomPos(0);
  CHECK_INVARIANT(feq(grd.getSpacing(), 0.5),
                  "Spacing of grid is not correct.");
  CHECK_INVARIANT(
      feq(grd.getNumX(), std::floor((fabs(bond.x) + 10.0) / 0.5 + 0.5)),
      "X length of grid is incorrect.");
  CHECK_INVARIANT(
      feq(grd.getNumY(), std::floor((fabs(bond.y) + 10.0) / 0.5 + 0.5)),
      "Y length of grid is incorrect.");
  CHECK_INVARIANT(
      feq(grd.getNumZ(), std::floor((fabs(bond.z) + 10.0) / 0.5 + 0.5)),
      "Z length of grid is incorrect.");
}

void test2CalculateDescriptors() {
  auto *data = new RealValueVect(0.0, 125);
  Point3D o(0.0, 0.0, 0.0);
  UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 1.0, &o, data);

  calculateDescriptors(grd, testfunctor());
  CHECK_INVARIANT(feq(data->getTotalVal(), grd.getSize()),
                  "Descriptor Calculation does not work correctly.");
}

void test3CubeFiles() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  RWMol mol = *v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  auto *data = new RealValueVect(0.0, 125);
  Point3D o(0.0, 0.0, 0.0);
  UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 1.0, &o, data);
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, double(i / 10.0));
  }
  writeToCubeFile(grd, mol, path + "test3.cube");
  UniformRealValueGrid3D grd2;
  auto mol2 = *readFromCubeFile(grd2, path + "test3.cube");

  CHECK_INVARIANT(grd.getSize() == grd2.getSize(),
                  "I/O: grid sizes are not the same.");
  for (unsigned int i = 0; i < grd2.getSize(); i++) {
    CHECK_INVARIANT(feq(grd2.getVal(i), double(i / 10.0)),
                    "I/O: values in grid are not correct.");
  }
  CHECK_INVARIANT(grd.getNumX() == grd2.getNumX(),
                  "I/O: grids are not the same.");
  CHECK_INVARIANT(grd.getNumY() == grd2.getNumY(),
                  "I/O: grids are not the same.");
  CHECK_INVARIANT(grd.getNumZ() == grd2.getNumZ(),
                  "I/O: grids are not the same.");
  CHECK_INVARIANT(feq(grd.getOffset().x, grd2.getOffset().x),
                  "I/O: grids are not the same.");
  CHECK_INVARIANT(feq(grd.getSpacing(), grd2.getSpacing()),
                  "I/O: grids are not the same.");
  CHECK_INVARIANT(grd.compareVectors(grd2), "I/O: grids are not the same.");
  CHECK_INVARIANT(grd.compareParams(grd2), "I/O: grids are not the same.");
  CHECK_INVARIANT(grd.compareGrids(grd2), "I/O: grids are not the same.");

  CHECK_INVARIANT(mol.getNumAtoms() == mol2.getNumAtoms(),
                  "I/O: number of atoms are not the same.");
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    CHECK_INVARIANT(mol.getAtomWithIdx(i)->getAtomicNum() ==
                        mol2.getAtomWithIdx(i)->getAtomicNum(),
                    "I/O: atoms are not the same");
    CHECK_INVARIANT(feq(mol.getConformer().getAtomPos(i).x,
                        mol2.getConformer().getAtomPos(i).x),
                    "I/O: atom positions are not the same");
    CHECK_INVARIANT(feq(mol.getConformer().getAtomPos(i).y,
                        mol2.getConformer().getAtomPos(i).y),
                    "I/O: atom positions are not the same");
    CHECK_INVARIANT(feq(mol.getConformer().getAtomPos(i).z,
                        mol2.getConformer().getAtomPos(i).z),
                    "I/O: atom positions are not the same");
  }
}

void test4Coulomb() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = *v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);

  computeGasteigerCharges(mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol.getConformer(0);
  for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
    charges.push_back(
        mol.getAtomWithIdx(i)->getProp<double>("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  auto grd = *constructGrid(mol);
  auto grd2 = *constructGrid(mol);

  Coulomb coul(mol);

  calculateDescriptors<Coulomb>(grd, coul);
  calculateDescriptors<Coulomb>(grd2, Coulomb(charges, pos));

  CHECK_INVARIANT(
      grd.compareGrids(grd2),
      "Coulomb: Different constructors do not yield the same descriptor.");
  CHECK_INVARIANT(feq(coul(0.0, 0.0, 0.0, 1000), 0.0),
                  "Coulomb: Potential between atoms wrong.(should be 0)");
  CHECK_INVARIANT(coul(2.0, 0.0, 0.0, 1000) < 0,
                  "Coulomb: Potential between positive charges not positive.");
  CHECK_INVARIANT(
      coul(-2.0, 0.0, 0.0, 1000) > 0,
      "Coulomb: Potential between positive and negative charges not negative.");
  CHECK_INVARIANT(feq(coul(0.0, 0.0, 0.0, 0.1), 0.0),
                  "Coulomb: Small threshold dist does not give 0.");

  calculateDescriptors<Coulomb>(grd, Coulomb(mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK_INVARIANT(grd.getVal(i) <= 0.0,
                    "Coulomb: Absolute value field not negative");
  }

  Coulomb coul1(mol, 0, -1.0, false, "_GasteigerCharge", 0.0, 0.01);
  CHECK_INVARIANT(coul1(-2.0, 0.0, 0.0, 1000) < 0,
                  "Coulomb: Potential between negative charges not positive.");
  CHECK_INVARIANT(
      coul1(2.0, 0.0, 0.0, 1000) > 0,
      "Coulomb: Potential between positive and negative charges not negative.");

  Coulomb coul2 = Coulomb(mol, 0, -.5, false, "_GasteigerCharge", 0.0, 0.01);
  CHECK_INVARIANT(
      coul1(-2.0, 0.0, 0.0, 1000) < coul2(-2.0, 0.0, 0.0, 1000),
      "Coulomb: Higher probecharge does not result in stronger forces.");

  Coulomb coul3(mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0);
  CHECK_INVARIANT(coul3(0.0, 0.0, 0.0, 1000) > coul3(0.1, 0.0, 0.0, 1000),
                  "Coulomb: Softcore interaction wrong.");
  CHECK_INVARIANT(coul3(0.66, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000),
                  "Coulomb: Softcore interaction wrong.");
  CHECK_INVARIANT(coul3(0.70, 0.0, 0.0, 1000) > coul3(0.68, 0.0, 0.0, 1000),
                  "Coulomb: Softcore interaction wrong.");
  CHECK_INVARIANT(feq(coul3(0.0, 0.0, 0.0, 0.1), 0.0),
                  "Coulomb: Small threshold dist does not give 0.");
}

void test5CoulombDielectric() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  RWMol mol = *v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);

  computeGasteigerCharges(mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol.getConformer(0);
  for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
    charges.push_back(
        mol.getAtomWithIdx(i)->getProp<double>("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  auto grd = *constructGrid(mol);
  auto grd2 = *constructGrid(mol);

  CoulombDielectric couldiele(mol);

  calculateDescriptors<CoulombDielectric>(grd, couldiele);
  calculateDescriptors<CoulombDielectric>(grd2,
                                          CoulombDielectric(charges, pos));

  CHECK_INVARIANT(
      grd.compareGrids(grd2),
      "CoulombDielectric: Different constructors do not yield the same descriptor.");
  CHECK_INVARIANT(
      feq(couldiele(0.0, 0.0, 0.0, 1000), 0.0),
      "CoulombDielectric: Potential between atoms wrong.(should be 0)");
  CHECK_INVARIANT(
      couldiele(2.0, 0.0, 0.0, 1000) < 0,
      "CoulombDielectric: Potential between positive charges not positive.");
  CHECK_INVARIANT(
      couldiele(-2.0, 0.0, 0.0, 1000) > 0,
      "CoulombDielectric: Potential between positive and negative charges not negative.");

  calculateDescriptors<CoulombDielectric>(grd,
                                          CoulombDielectric(mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK_INVARIANT(grd.getVal(i) <= 0.0,
                    "CoulombDielectric: Absolute value field not negative");
  }

  CHECK_INVARIANT(
      feq(couldiele(0.0, 0.0, 0.0, 1000), 0.0),
      "CoulombDielectric: Potential between atoms wrong.(should be 0)");
  CHECK_INVARIANT(
      couldiele(2.0, 0.0, 0.0, 1000) < 0,
      "CoulombDielectric: Potential between positive charges not positive.");
  CHECK_INVARIANT(
      couldiele(-2.0, 0.0, 0.0, 1000) > 0,
      "CoulombDielectric: Potential between positive and negative charges not negative.");

  calculateDescriptors<CoulombDielectric>(grd,
                                          CoulombDielectric(mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK_INVARIANT(grd.getVal(i) <= 0.0,
                    "CoulombDielectric: Absolute value field not negative");
  }

  CoulombDielectric couldiele1(mol, 0, -1.0);
  CHECK_INVARIANT(
      couldiele1(-2.0, 0.0, 0.0, 1000) < 0,
      "CoulombDielectric: Potential between negative charges not positive.");
  CHECK_INVARIANT(
      couldiele1(2.0, 0.0, 0.0, 1000) > 0,
      "CoulombDielectric: Potential between positive and negative charges not negative.");

  CoulombDielectric couldiele2(mol, 0, -.5);

  CHECK_INVARIANT(
      couldiele1(-2.0, 0.0, 0.0, 1000) < couldiele2(-2.0, 0.0, 0.0, 1000),
      "CoulombDielectric: Higher probecharge does not result in stronger forces.");

  CoulombDielectric couldiele3(mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0);
  CHECK_INVARIANT(
      couldiele3(0.0, 0.0, 0.0, 1000) > couldiele3(0.1, 0.0, 0.0, 1000),
      "CoulombDielectric: Softcore interaction wrong.");
  CHECK_INVARIANT(
      couldiele3(0.66, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000),
      "CoulombDielectric: Softcore interaction wrong.");
  CHECK_INVARIANT(
      couldiele3(0.70, 0.0, 0.0, 1000) > couldiele3(0.68, 0.0, 0.0, 1000),
      "CoulombDielectric: Softcore interaction wrong.");

  mol = *v2::FileParsers::MolFromMolFile(path + "glucose.mol", fopts);
  computeGasteigerCharges(mol);

  CoulombDielectric couldiele4(mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 80.0, 4.0);
  CoulombDielectric couldiele5(mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 200.0, 4.0);
  CoulombDielectric couldiele6(mol, 0, 1.0, false, "_GasteigerCharge", 0.01,
                               1.0, 80.0, 10.0);
  CHECK_INVARIANT(
      couldiele5(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000),
      "CoulombDielectric: solvent permittivity scaling wrong.");

  CHECK_INVARIANT(
      couldiele6(-1.0, 0.0, 0.0, 1000) < couldiele4(-1.0, 0.0, 0.0, 1000),
      "CoulombDielectric: solute permittivity scaling wrong.");
}

void test6VdWaals() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = *v2::FileParsers::MolFromMolFile(path + "HCl.mol", fopts);
  try {
    VdWaals vdw = constructVdWaalsMMFF(mol, 0, 6, false, 1.0);
  } catch (ValueErrorException &dexp) {
    BOOST_LOG(rdInfoLog) << "Expected failure: " << dexp.what() << "\n";
  }

  mol = *v2::FileParsers::MolFromMolFile(path + "HCN.mol", fopts);
  VdWaals vdw = constructVdWaalsMMFF(mol, 0, 6, false, 1.0);

  CHECK_INVARIANT(vdw(-5.0, 0, 0, 1000) < 0,
                  "VdWMMFF: Potential not negative in favorable region.");
  CHECK_INVARIANT(
      vdw(-1.68, 0, 0, 1000) > vdw(-5.0, 0, 0, 1000),
      "VdWMMFF: Potential next to core not higher than further apart.");
  CHECK_INVARIANT(
      vdw(-5.0, 0, 0, 1000) < vdw(-10.0, 0, 0, 1000),
      "VdWMMFF: Potential very far apart not higher than in favorable distance to core.");

  auto mol2 = *v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);
  vdw = constructVdWaalsMMFF(mol2, 0, 6, false, 1.0);
  VdWaals vdw2 = constructVdWaalsMMFF(mol2, 0, 6, true, 1.0);
  CHECK_INVARIANT(fabs(vdw2(-3.0, 0, 0, 1000) - vdw(-3.0, 0, 0, 1000)) > 0.0001,
                  "VdWMMFF: No scaling of interactions.");

  VdWaals vdw3 = constructVdWaalsUFF(mol, 0, "O_3", 1.0);
  CHECK_INVARIANT(vdw3(-5.0, 0, 0, 1000) < 0,
                  "VdWMMFF: Potential not negative in favorable region.");
  CHECK_INVARIANT(
      vdw3(-1.68, 0, 0, 1000) > vdw3(-5.0, 0, 0, 1000),
      "VdWUFF: Potential next to core not higher than further apart.");
  CHECK_INVARIANT(
      vdw3(-5.0, 0, 0, 1000) < vdw3(-10.0, 0, 0, 1000),
      "VdWUFF: Potential very far apart not higher than in favorable distance to core.");

  std::string names[] = {
      "acetone",       "aceticacid",    "phenol",       "phenolate",
      "serine",        "threonine",     "ethanol",      "diethylether",
      "h2o",           "ammonia",       "ethylamine",   "imine",
      "acetonitrile",  "histidine",     "phenylamine",  "methylammonium",
      "fluoromethane", "chloromethane", "bromomethane", "glycine",
      "glyphe",        "glysergly",     "glythrgly",    "glucose"};
  for (const auto & name : names) {
    mol = *v2::FileParsers::MolFromMolFile(path + name + ".mol", fopts);
    vdw = constructVdWaalsMMFF(mol);
    CHECK_INVARIANT(vdw(0.0, 0.0, 0.0, 1000),
                    "VdWMMFF: crashed with " + name);
    vdw = constructVdWaalsUFF(mol);
    CHECK_INVARIANT(vdw(0.0, 0.0, 0.0, 1000),
                    "VdWUFF: crashed with " + name);
  }
}

void test7HBond() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  // Generate Molecule with 3D Coordinates and Partial Charges
  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol =
      *v2::FileParsers::MolFromMolFile(path + "ethane.mol", fopts);  // Ethane
  auto grd = *constructGrid(mol, 0, 5.0, 1);
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, 1.0);
  }

  UniformRealValueGrid3D grd1(grd);

  HBond hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  TEST_ASSERT(!grd.compareGrids(grd1));
  const RealValueVect *vect = grd.getOccupancyVect();
  const RealValueVect *vect1 = grd1.getOccupancyVect();
  CHECK_INVARIANT(
      (unsigned int)fabs(((*vect1 - *vect).getTotalVal())) == grd.getSize(),
      "");

  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");

  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, 1.0);
  }
  calculateDescriptors<HBond>(grd, hbonddes);
  TEST_ASSERT(!grd.compareGrids(grd1));
  CHECK_INVARIANT((unsigned int)fabs(((*vect1 - *vect).getTotalVal())), "");

  mol = *v2::FileParsers::MolFromMolFile(path + "aceticacid.mol",
                                         fopts);  // Acetic Acid
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  HBond hbonddes1(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes1.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes1);
  HBond hbonddes2(mol, 0, "O", false, 0.001);
  CHECK_INVARIANT(
      hbonddes1(4.0, 0.0, 1.0, 1000) > hbonddes2(4.0, 0.0, 1.0, 1000),
      "HBond: Flexible bonds do not yield more negative potential.");
  HBond hbonddes3(mol, 0, "NH", true, 0.001);
  CHECK_INVARIANT(hbonddes3.getNumInteractions() == 2, "");
  CHECK_INVARIANT(
      hbonddes(2.0, 2.0, 2.0, 1000) < hbonddes3(2.0, 2.0, 2.0, 1000),
      "HBond: N probe stronger interaction than O probe");
  HBond hbonddes4(mol, 0, "N", true, 0.001);
  CHECK_INVARIANT(hbonddes4.getNumInteractions() == 1, "");
  CHECK_INVARIANT(
      hbonddes1(3.0, 0.0, 0.0, 1000) < hbonddes4(3.0, 0.0, 0.0, 1000),
      "HBond: N probe stronger interaction than O probe");

  mol =
      *v2::FileParsers::MolFromMolFile(path + "acetone.mol", fopts);  // Acetone
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "diethylether.mol",
                                         fopts);  // Et2O
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);  // H2O
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol =
      *v2::FileParsers::MolFromMolFile(path + "ammonia.mol", fopts);  // ammonia
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "imine.mol", fopts);  // imine
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "methylammonium.mol",
                                         fopts);  // methylammonium
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "chloromethane.mol",
                                         fopts);  // Chloromethane
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "phosphonate.mol",
                                         fopts);  // Phosphonate
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "phosphatediester.mol",
                                         fopts);  // Phosphatediester
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 4, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "hydrogenphosphatediester.mol",
                                         fopts);  // Hydrogenphosphatediester
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 4, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "mustardgas.mol",
                                         fopts);  // mustard gas
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "alicin.mol", fopts);  // Alicin
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *v2::FileParsers::MolFromMolFile(path + "sulfanilamide.mol",
                                         fopts);  // Sulfanilamide
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 4, "");
  calculateDescriptors<HBond>(grd, hbonddes);
}

void test8Hydrophilic() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MolInteractionFields/test_data/";

  auto fopts = v2::FileParsers::MolFileParserParams();
  fopts.removeHs = false;
  auto mol = *v2::FileParsers::MolFromMolFile(path + "h2o.mol", fopts);

  Hydrophilic hydro(mol);
  HBond hbondOH(mol, 0, "OH");
  HBond hbondO(mol, 0, "O");

  double hyd = hydro(0.0, 0.0, 0.0, 1000), hOH = hbondOH(0.0, 0.0, 0.0, 1000),
         hO = hbondO(0.0, 0.0, 0.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(1.0, 1.5, 2.0, 1000);
  hOH = hbondOH(1.0, 1.5, 2.0, 1000);
  hO = hbondO(1.0, 1.5, 2.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(2.0, 1.5, -3.0, 1000);
  hOH = hbondOH(2.0, 1.5, -3.0, 1000);
  hO = hbondO(2.0, 1.5, -3.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(-2.5, 0.5, 3.0, 1000);
  hOH = hbondOH(-2.5, 0.5, 3.0, 1000);
  hO = hbondO(-2.5, 0.5, 3.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(10.0, 1.5, 1.0, 1000);
  hOH = hbondOH(10.0, 1.5, 1.0, 1000);
  hO = hbondO(10.0, 1.5, 1.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(6.0, -5.0, 0.0, 1000);
  hOH = hbondOH(6.0, -5.0, 0.0, 1000);
  hO = hbondO(6.0, -5.0, 0.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(-3.0, -3.0, 7.0, 1000);
  hOH = hbondOH(-3.0, -3.0, 7.0, 1000);
  hO = hbondO(-3.0, -3.0, 7.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(1.0, 0.0, 0.0, 1000);
  hOH = hbondOH(1.0, 0.0, 0.0, 1000);
  hO = hbondO(1.0, 0.0, 0.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(0.0, 2.0, 2.0, 1000);
  hOH = hbondOH(0.0, 2.0, 2.0, 1000);
  hO = hbondO(0.0, 2.0, 2.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(2.0, -2.0, 0.0, 1000);
  hOH = hbondOH(2.0, -2.0, 0.0, 1000);
  hO = hbondO(2.0, -2.0, 0.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");

  hyd = hydro(2.0, -2.0, -3.0, 1000);
  hOH = hbondOH(2.0, -2.0, -3.0, 1000);
  hO = hbondO(2.0, -2.0, -3.0, 1000);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd),
                  "Hydrophilic: Not working correctly.");
}

int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing MIF\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1: Constructing Grid \n\n";
  test1ConstructGrid();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2: Calculating Descriptors \n\n";
  test2CalculateDescriptors();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test3: I/O with .cube Files \n\n";
  test3CubeFiles();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test4: Coulomb (vacuum) descriptor \n\n";
  test4Coulomb();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test5: Coulomb (dielectric) descriptor \n\n";
  test5CoulombDielectric();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test6: van der Waals descriptor \n\n";
  test6VdWaals();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test7: Hydrogen Bond Descriptors \n\n";
  test7HBond();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test8: Hydrophilic Descriptors \n\n";
  test8Hydrophilic();

  return 0;
}
