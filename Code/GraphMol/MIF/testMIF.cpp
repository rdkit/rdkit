//  testMIF.cpp
//  Created on: Apr 4, 2014
//  Author: hahnda6
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/point.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/MIF/MIFDescriptors.h>
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
    testfunctor() {
    }
    double operator()(const Point3D &pt) {
      return 1.0;
    }
    double operator()(const double &x, const double &y, const double &z) {
      return 1.0;
    }
  };
}

void test1ConstructGrid() {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";
  //	//Generate Molecule with 3D Coordinates and Partial Charges
  //	RWMol mol=*SmilesToMol("Cl");
  //	MolOps::addHs(mol);
  //	DGeomHelpers::EmbedMolecule(mol);

  //	MolToMolFile(mol, path + "HCl.mol");

  RWMol mol = *MolFileToMol(path + "HCl.mol", true, false);
  UniformRealValueGrid3D grd = *constructGrid(mol, 0, 5.0, 0.5);

  Point3D bond = mol.getConformer().getAtomPos(1) - mol.getConformer().getAtomPos(0);
  CHECK_INVARIANT(feq(grd.getSpacing (), 0.5),
                  "Spacing of grid is not correct.");
  CHECK_INVARIANT(feq(grd.getNumX (), std::floor ((fabs (bond.x) + 10.0) / 0.5 + 0.5)),
                  "X length of grid is incorrect.");
  CHECK_INVARIANT(feq(grd.getNumY (), std::floor ((fabs (bond.y) + 10.0) / 0.5 + 0.5)),
                  "Y length of grid is incorrect.");
  CHECK_INVARIANT(feq(grd.getNumZ (), std::floor ((fabs (bond.z) + 10.0) / 0.5 + 0.5)),
                  "Z length of grid is incorrect.");
}

void test2CalculateDescriptors () {
  RealValueVect *data = new RealValueVect (0.0, 125);
  UniformRealValueGrid3D grd (5.0, 5.0, 5.0, 1.0, new Point3D (0.0, 0.0, 0.0),
                              data);

  calculateDescriptors (grd, testfunctor());
  CHECK_INVARIANT(feq (data->getTotalVal (), grd.getSize ()),
                  "Descriptor Calculation does not work correctly.");
}

void test3CubeFiles () {
  std::string path = getenv ("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";

  RWMol mol = *MolFileToMol(path + "HCl.mol", true, false);
  RealValueVect *data = new RealValueVect (0.0, 125);
  UniformRealValueGrid3D grd (5.0, 5.0, 5.0, 1.0, new Point3D (0.0, 0.0, 0.0),
                              data);
  for (unsigned int i = 0; i < grd.getSize (); i++) {
    grd.setVal (i, double (i / 10.0));
  }
  writeToCubeFile (grd, mol, path + "testCube.cube");
  UniformRealValueGrid3D grd2;
  RWMol mol2 = *readFromCubeFile (grd2, path + "testCube.cube");

  CHECK_INVARIANT(grd.getSize () == grd2.getSize (),
                  "I/O: grid sizes are not the same.");
  for (unsigned int i = 0; i < grd2.getSize (); i++) {
    CHECK_INVARIANT(feq (grd2.getVal (i), double (i / 10.0)),
                    "I/O: values in grid are not correct.");
  }
  CHECK_INVARIANT(grd.compareGrids (grd2), "I/O: grids are not the same.");

  CHECK_INVARIANT(mol.getNumAtoms () == mol2.getNumAtoms (),
                  "I/O: number of atoms are not the same.");
  for (unsigned int i = 0; i < mol.getNumAtoms (); i++) {
    CHECK_INVARIANT(mol.getAtomWithIdx (i)->getAtomicNum ()
                    == mol2.getAtomWithIdx (i)->getAtomicNum (),
                    "I/O: atoms are not the same");
    CHECK_INVARIANT(feq (mol.getConformer ().getAtomPos (i).x,
                         mol2.getConformer ().getAtomPos (i).x),
                    "I/O: atom positions are not the same");
    CHECK_INVARIANT(feq (mol.getConformer ().getAtomPos (i).y,
                         mol2.getConformer ().getAtomPos (i).y),
                    "I/O: atom positions are not the same");
    CHECK_INVARIANT(feq (mol.getConformer ().getAtomPos (i).z,
                         mol2.getConformer ().getAtomPos (i).z),
                    "I/O: atom positions are not the same");
  }
}

void test4Coulomb () {
  std::string path = getenv("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";

  RWMol mol = *MolFileToMol(path + "HCl.mol", true, false);

  computeGasteigerCharges(mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol.getConformer(0);
  for (int i = 0; i < mol.getNumAtoms(); ++i) {
    charges.push_back(mol.getAtomWithIdx (i)->getProp<double> ("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  UniformRealValueGrid3D grd = *constructGrid(mol);
  UniformRealValueGrid3D grd2 = *constructGrid(mol);

  Coulomb coul(mol);

  calculateDescriptors<Coulomb>(grd, coul);
  calculateDescriptors<Coulomb>(grd2, Coulomb (charges, pos));

  CHECK_INVARIANT(grd.compareGrids (grd2),
                  "Coulomb: Different constructors do not yield the same descriptor.");
  CHECK_INVARIANT(feq (coul (Point3D (0.0, 0.0, 0.0)), 0.0),
                  "Coulomb: Potential between atoms wrong.(should be 0)");
  CHECK_INVARIANT(coul (Point3D (2.0, 0.0, 0.0)) < 0,
                  "Coulomb: Potential between positive charges not positive.");
  CHECK_INVARIANT(coul (Point3D (-2.0, 0.0, 0.0)) > 0,
                  "Coulomb: Potential between positive and negative charges not negative.");

  calculateDescriptors<Coulomb>(grd, Coulomb(mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK_INVARIANT(grd.getVal (i) <= 0.0, "Coulomb: Absolute value field not negative");
  }

  Coulomb coul1(mol, 0, -1.0);
  CHECK_INVARIANT(coul1(Point3D(-2.0, 0.0, 0.0)) < 0, "Coulomb: Potential between negative charges not positive.");
  CHECK_INVARIANT(coul1 (Point3D (2.0, 0.0, 0.0)) > 0, "Coulomb: Potential between positive and negative charges not negative.");

  Coulomb coul2 = Coulomb(mol, 0, -.5);
  CHECK_INVARIANT(coul1 (Point3D (-2.0, 0.0, 0.0)) < coul2 (Point3D (-2.0, 0.0, 0.0)),
                  "Coulomb: Higher probecharge does not result in stronger forces.");

  Coulomb coul3(mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0);
  CHECK_INVARIANT(coul3(Point3D(0.0, 0.0, 0.0)) > coul3(Point3D(0.1, 0.0, 0.0)),
                  "Coulomb: Softcore interaction wrong.");
  CHECK_INVARIANT(coul3(Point3D(0.66, 0.0, 0.0)) > coul3(Point3D(0.68, 0.0, 0.0)),
                  "Coulomb: Softcore interaction wrong.");
  CHECK_INVARIANT(coul3(Point3D(0.70, 0.0, 0.0)) > coul3(Point3D(0.68, 0.0, 0.0)),
                  "Coulomb: Softcore interaction wrong.");
}

void test5CoulombDielectric () {
  std::string path = getenv ("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";

  RWMol mol = *MolFileToMol(path + "HCl.mol", true, false);

  computeGasteigerCharges(mol);

  std::vector<double> charges;
  std::vector<Point3D> pos;
  Conformer conf = mol.getConformer(0);
  for (int i = 0; i < mol.getNumAtoms(); ++i) {
    charges.push_back(mol.getAtomWithIdx(i)->getProp<double>("_GasteigerCharge"));
    pos.push_back(conf.getAtomPos(i));
  }

  UniformRealValueGrid3D grd = *constructGrid(mol);
  UniformRealValueGrid3D grd2 = *constructGrid(mol);

  CoulombDielectric couldiele(mol);

  calculateDescriptors<CoulombDielectric>(grd, couldiele);
  calculateDescriptors<CoulombDielectric>(grd2, CoulombDielectric (charges, pos));

  CHECK_INVARIANT(grd.compareGrids(grd2), "CoulombDielectric: Different constructors do not yield the same descriptor.");
  CHECK_INVARIANT(feq (couldiele (Point3D (0.0, 0.0, 0.0)), 0.0), "CoulombDielectric: Potential between atoms wrong.(should be 0)");
  CHECK_INVARIANT(couldiele (Point3D (2.0, 0.0, 0.0)) < 0, "CoulombDielectric: Potential between positive charges not positive.");
  CHECK_INVARIANT(couldiele (Point3D (-2.0, 0.0, 0.0)) > 0, "CoulombDielectric: Potential between positive and negative charges not negative.");

  calculateDescriptors<CoulombDielectric>(grd, CoulombDielectric (mol, 0, 1.0, true));
  for (unsigned int i = 0; i < grd.getSize(); i++) {
    CHECK_INVARIANT(grd.getVal (i) <= 0.0, "CoulombDielectric: Absolute value field not negative");
  }

  CoulombDielectric couldiele1(mol, 0, -1.0);
  CHECK_INVARIANT(couldiele1 (Point3D (-2.0, 0.0, 0.0)) < 0, "CoulombDielectric: Potential between negative charges not positive.");
  CHECK_INVARIANT(couldiele1 (Point3D (2.0, 0.0, 0.0)) > 0, "CoulombDielectric: Potential between positive and negative charges not negative.");

  CoulombDielectric couldiele2 (mol, 0, -.5);

  CHECK_INVARIANT(couldiele1 (Point3D (-2.0, 0.0, 0.0)) < couldiele2 (Point3D (-2.0, 0.0, 0.0)),
                  "CoulombDielectric: Higher probecharge does not result in stronger forces.");

  CoulombDielectric couldiele3(mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0);
  CHECK_INVARIANT(couldiele3(Point3D (0.0, 0.0, 0.0)) > couldiele3 (Point3D (0.1, 0.0, 0.0)),
                  "CoulombDielectric: Softcore interaction wrong.");
  CHECK_INVARIANT(couldiele3(Point3D (0.66, 0.0, 0.0)) > couldiele3 (Point3D (0.68, 0.0, 0.0)),
                  "CoulombDielectric: Softcore interaction wrong.");
  CHECK_INVARIANT(couldiele3(Point3D (0.70, 0.0, 0.0)) > couldiele3 (Point3D (0.68, 0.0, 0.0)),
                  "CoulombDielectric: Softcore interaction wrong.");



  mol = *MolFileToMol(path + "glucose.mol", true, false);
  computeGasteigerCharges(mol);

  CoulombDielectric couldiele4 (mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0,  80.0,  4.0);
  CoulombDielectric couldiele5 (mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0, 200.0,  4.0);
  CoulombDielectric couldiele6 (mol, 0, 1.0, false, "_GasteigerCharge", 0.01, 1.0,  80.0, 10.0);
  CHECK_INVARIANT(couldiele5 (Point3D (-1.0, 0.0, 0.0)) < couldiele4 (Point3D (-1.0, 0.0, 0.0)),
                  "CoulombDielectric: solvent permittivity scaling wrong.");


  CHECK_INVARIANT(couldiele6 (Point3D (-1.0, 0.0, 0.0)) < couldiele4 (Point3D (-1.0, 0.0, 0.0)),
                  "CoulombDielectric: solute permittivity scaling wrong.");
}

void test6VdWaals () {
  std::string path = getenv ("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";

  RWMol mol = *MolFileToMol(path + "HCl.mol", true, false);

  VdWaals vdw = constructVdWaalsMMFF(mol, 0, 6, false, 1.0);

  CHECK_INVARIANT(vdw(Point3D(-5.0, 0, 0)) < 0, "VdWMMFF: Potential not negative in favorable region.");
  CHECK_INVARIANT(vdw(Point3D(-1.68, 0, 0)) > vdw(Point3D(-5.0, 0, 0)), "VdWMMFF: Potential next to core not higher than further apart.");
  CHECK_INVARIANT(vdw(Point3D(-5.0, 0, 0)) < vdw(Point3D(-10.0, 0, 0)), "VdWMMFF: Potential very far apart not higher than in favorable distance to core.");

  RWMol mol2 = *MolFileToMol(path + "h2o.mol", true, false);
  vdw  = constructVdWaalsMMFF(mol2, 0, 6, false, 1.0);
  VdWaals vdw2 = constructVdWaalsMMFF(mol2, 0, 6, true, 1.0);
  CHECK_INVARIANT(fabs(vdw2(Point3D(-3.0, 0, 0)) - vdw(Point3D(-3.0, 0, 0))) > 0.0001,
                  "VdWMMFF: No scaling of interactions.");


  VdWaals vdw3 = constructVdWaalsUFF(mol, 0, "O_3", 1.0);
  CHECK_INVARIANT(vdw3(Point3D(-5.0, 0, 0)) < 0, "VdWMMFF: Potential not negative in favorable region.");
  CHECK_INVARIANT(vdw3(Point3D(-1.68, 0, 0)) > vdw3(Point3D(-5.0, 0, 0)),
                  "VdWUFF: Potential next to core not higher than further apart.");
  CHECK_INVARIANT(vdw3(Point3D(-5.0, 0, 0)) < vdw3(Point3D(-10.0, 0, 0)),
                  "VdWUFF: Potential very far apart not higher than in favorable distance to core.");

  std::string names[] = { "acetone", "aceticacid", "phenol", "phenolate",
      "serine", "threonine", "ethanol", "diethylether", "h2o", "ammonia",
      "ethylamine", "imine", "acetonitrile", "histidine", "phenylamine",
      "methylammonium", "fluoromethane", "chloromethane", "bromomethane",
      "glycine", "glyphe", "glysergly", "glythrgly", "glucose" };
  for ( unsigned int i = 0; i < 24; i++ ){
    mol = *MolFileToMol(path + names[i] + ".mol", true, false);
    vdw = constructVdWaalsMMFF(mol);
    CHECK_INVARIANT(vdw(Point3D(0.0,0.0,0.0)), "VdWMMFF: crashed with " + names[i]);
    vdw = constructVdWaalsUFF(mol);
    CHECK_INVARIANT(vdw(Point3D(0.0,0.0,0.0)), "VdWUFF: crashed with " + names[i]);
  }
}

void test7HBond() {
  std::string path = getenv ("RDBASE");
  path += "/Code/GraphMol/MIF/test_data/";

  UniformRealValueGrid3D grd;
  RWMol mol;

  //Generate Molecule with 3D Coordinates and Partial Charges
  mol = *MolFileToMol(path + "ethane.mol", true, false);              //Ethane
  grd = *constructGrid(mol, 0, 5.0, 1);
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
  CHECK_INVARIANT((unsigned int )fabs(((*vect1 - *vect).getTotalVal())) == grd.getSize(),"");

  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");

  for (unsigned int i = 0; i < grd.getSize(); i++) {
    grd.setVal(i, 1.0);
  }
  calculateDescriptors<HBond>(grd, hbonddes);
  TEST_ASSERT(!grd.compareGrids(grd1));
  CHECK_INVARIANT((unsigned int)fabs(((*vect1 - *vect).getTotalVal())), "");

  mol = *MolFileToMol(path + "aceticacid.mol", true, false);          //Acetic Acid
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  HBond hbonddes1(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes1.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes1);
  HBond hbonddes2(mol, 0, "O", false, 0.001);
  CHECK_INVARIANT(hbonddes1(Point3D(4.0, 0.0, 1.0)) > hbonddes2(Point3D(4.0, 0.0, 1.0)),
                  "HBond: Flexible bonds do not yield more negative potential.");
  HBond hbonddes3(mol, 0, "NH", true, 0.001);
  CHECK_INVARIANT(hbonddes3.getNumInteractions() == 2, "");
  CHECK_INVARIANT(hbonddes(Point3D(2.0, 2.0, 1.0)) < hbonddes3(Point3D(2.0, 2.0, 1.0)), "HBond: N probe stronger interaction than O probe");
  HBond hbonddes4(mol, 0, "N", true, 0.001);
  CHECK_INVARIANT(hbonddes4.getNumInteractions() == 1, "");
  CHECK_INVARIANT(hbonddes1(Point3D(3.0, 0.0,0.0)) < hbonddes4(Point3D(3.0,0.0,0.0)), "HBond: N probe stronger interaction than O probe");

  mol = *MolFileToMol(path + "acetone.mol", true, false);     //Acetone
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "diethylether.mol", true, false);                //Et2O
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "h2o.mol", true, false);         //H2O
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "ammonia.mol", true, false);             //ammonia
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);


  mol = *MolFileToMol(path + "imine.mol", true, false);               //imine
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "methylammonium.mol", true, false);//methylammonium
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "chloromethane.mol", true, false);       //Chloromethane
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "phosphonate.mol", true, false);         //Phosphonate
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 3, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "phosphatediester.mol", true, false);//Phosphatediester
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 4, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "hydrogenphosphatediester.mol", true, false);//Hydrogenphosphatediester
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 4, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "mustardgas.mol", true, false);          //mustard gas
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 2, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "alicin.mol", true, false);              //Alicin
  grd = *constructGrid(mol, 0, 5.0, 1);
  hbonddes = HBond(mol, 0, "OH", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 1, "");
  calculateDescriptors<HBond>(grd, hbonddes);
  hbonddes = HBond(mol, 0, "O", true, 0.001);
  CHECK_INVARIANT(hbonddes.getNumInteractions() == 0, "");
  calculateDescriptors<HBond>(grd, hbonddes);

  mol = *MolFileToMol(path + "sulfanilamide.mol", true, false);       //Sulfanilamide
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
  path += "/Code/GraphMol/MIF/test_data/";

  RWMol mol = *MolFileToMol(path + "h2o.mol", true, false);

  Hydrophilic hydro(mol);
  HBond hbondOH(mol, 0, "OH");
  HBond hbondO(mol, 0, "O");

  Point3D pt(0.0,0.0,0.0);
  double hyd = hydro(pt),
      hOH = hbondOH(pt),
      hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(1.0, 1.5, 2.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(2.0, 1.5, -3.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(-2.5, 0.5, 3.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(10.0, 1.5, 1.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(6.0, -5.0, 0.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(-3.0, -3.0, 7.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(1.0, 0.0, 0.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(0.0, 2.0, 2.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(2.0, -2.0, 0.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");

  pt = Point3D(2.0, -2.0, -3.0);
  hyd = hydro(pt);
  hOH = hbondOH(pt);
  hO  = hbondO(pt);
  CHECK_INVARIANT(feq(std::min(hOH, hO), hyd), "Hydrophilic: Not working correctly.");
}

int main () {
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
