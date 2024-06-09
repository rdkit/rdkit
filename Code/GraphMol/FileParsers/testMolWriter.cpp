//
//   Copyright (C) 2002-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "MolSupplier.h"
#include "MolWriters.h"
#include "FileParsers.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

void testSmilesWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    names.push_back(mname);
    props.push_back(pval);
    writer->write(*mol);
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  writer->close();
  delete writer;
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname);
  int i = 0;
  mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    CHECK_INVARIANT(mname == names[i], "");
    CHECK_INVARIANT(pval == props[i], "");
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  delete nSup;
}

void testSmilesWriter2() {
  {
    std::stringstream ss;
    bool takeOwnership = false, includeHeader = false, isomericSmiles = false;
    SmilesWriter *writer = new SmilesWriter(&ss, " ", "Name", takeOwnership,
                                            includeHeader, isomericSmiles);
    RWMol *mol;

    mol = SmilesToMol("c1ccccc1");
    // MolOps::Kekulize(*mol);
    writer->write(*mol);
    delete mol;

    mol = SmilesToMol("F[C@H](Cl)Br");
    writer->write(*mol);
    delete mol;
    writer->close();
    TEST_ASSERT(ss.str() == "c1ccccc1 0\nFC(Cl)Br 1\n");
    delete writer;
  }
  {
    std::stringstream ss;
    bool takeOwnership = false, includeHeader = false, isomericSmiles = true;
    SmilesWriter *writer = new SmilesWriter(&ss, " ", "Name", takeOwnership,
                                            includeHeader, isomericSmiles);
    RWMol *mol;

    mol = SmilesToMol("c1ccccc1");
    MolOps::Kekulize(*mol);
    writer->write(*mol);
    delete mol;

    mol = SmilesToMol("F[C@H](Cl)Br");
    writer->write(*mol);
    delete mol;
    writer->close();
    TEST_ASSERT(ss.str() == "C1=CC=CC=C1 0\nF[C@H](Cl)Br 1\n");
    delete writer;
  }
}

void testSmilesWriterNoNames() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles_molwriter.csv";

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ", "");
  writer->setProps(propNames);

  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("Column_2", pval);
    mol->setProp(common_properties::_Name, "bogus");
    props.push_back(pval);
    writer->write(*mol);
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  writer->flush();
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname, " ", 0, -1);
  int i = 0;
  mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    delete mol;
    TEST_ASSERT(mname != "bogus");
    TEST_ASSERT(pval == props[i]);
    i++;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  TEST_ASSERT(writer->numMols() == nSup->length());
  writer->close();
  delete writer;
  delete nSup;
}

void testSmilesWriterClose() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles_molwriter.csv";

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ", "");
  writer->setProps(propNames);

  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("Column_2", pval);
    mol->setProp(common_properties::_Name, "bogus");
    props.push_back(pval);
    writer->write(*mol);
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  writer->flush();
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname, " ", 0, -1);
  int i = 0;
  mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    delete mol;
    TEST_ASSERT(mname != "bogus");
    TEST_ASSERT(pval == props[i]);
    i++;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  TEST_ASSERT(writer->numMols() == nSup->length());
  writer->close();
  delete nSup;
  delete writer;
}

void testSDWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);

  std::string ofile =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/outNCI_few.sdf_molwriter.sdf";
  auto *writer = new SDWriter(ofile);

  STR_VECT names;

  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    names.push_back(mname);

    writer->write(*mol);
    delete mol;
  }
  writer->flush();
  CHECK_INVARIANT(writer->numMols() == 16, "");

  // make sure we can close() the writer and delete it:
  writer->close();
  delete writer;

  // now read in the file we just finished writing
  SDMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == names[i], "");

    delete mol;
    i++;
  }

  // now read in a file with aromatic information on the bonds
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
  SDMolSupplier nreader(infile);
  i = 0;
  while (!nreader.atEnd()) {
    ROMol *mol = nreader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == names[i], "");
    i++;

    delete mol;
  }
}

void testTDTWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);

  std::string ofile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few_molwriter.tdt";
  auto *writer = new TDTWriter(ofile);

  STR_VECT names;
  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp("CAS_RN", mname);
    names.push_back(mname);

    writer->write(*mol);
    delete mol;
  }
  TEST_ASSERT(writer->numMols() == 16);
  writer->close();
  delete writer;

  // now read in the file we just finished writing
  TDTMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    if (mol) {
      std::string mname;
      mol->getProp("CAS_RN", mname);
      CHECK_INVARIANT(mname == names[i], "");
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i == 16);
}

void testSmilesWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles_molwriter.csv";
  auto *oStream = new std::ofstream(oname.c_str());

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oStream, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    names.push_back(mname);
    props.push_back(pval);
    writer->write(*mol);
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  writer->close();
  delete writer;
  delete nSup;
  delete oStream;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname);
  int i = 0;
  mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    CHECK_INVARIANT(mname == names[i], "");
    CHECK_INVARIANT(pval == props[i], "");
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  delete nSup;
}

void testSDWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  {
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);

    std::string ofile =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/outNCI_few_molwriter.sdf";
    auto *oStream = new std::ofstream(ofile.c_str());

    auto *writer = new SDWriter(oStream);

    STR_VECT names;

    while (!sdsup.atEnd()) {
      ROMol *mol = sdsup.next();
      std::string mname;
      mol->getProp(common_properties::_Name, mname);
      names.push_back(mname);

      writer->write(*mol);
      delete mol;
    }
    writer->flush();
    CHECK_INVARIANT(writer->numMols() == 16, "");
    writer->close();
    delete writer;
    delete oStream;

    // now read in the file we just finished writing
    SDMolSupplier reader(ofile);
    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      std::string mname;
      mol->getProp(common_properties::_Name, mname);
      CHECK_INVARIANT(mname == names[i], "");

      delete mol;
      i++;
    }
  }

  {
    // now read in a file with aromatic information on the bonds
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
    SDMolSupplier nreader(infile);
    unsigned int i = 0;
    while (!nreader.atEnd()) {
      ROMol *mol = nreader.next();
      TEST_ASSERT(mol);
      ++i;
      delete mol;
    }
    TEST_ASSERT(i == 16);
  }
}

void testTDTWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);

  std::string ofile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few_molwriter.tdt";
  auto *oStream = new std::ofstream(ofile.c_str());
  auto *writer = new TDTWriter(oStream);

  STR_VECT names;

  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp("CAS_RN", mname);
    names.push_back(mname);

    writer->write(*mol);
    delete mol;
  }
  writer->flush();
  TEST_ASSERT(writer->numMols() == 16);
  writer->close();
  delete writer;
  delete oStream;

  // now read in the file we just finished writing
  TDTMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    if (mol) {
      std::string mname;
      mol->getProp("CAS_RN", mname);
      CHECK_INVARIANT(mname == names[i], "");
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i == 16);
}

void testSDMemoryCorruption() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Data/NCI/first_200.props.sdf";
  SDMolSupplier sdsup(fname, true);
  std::string ofile = rdbase +
                      "/Code/GraphMol/FileParsers/test_data/"
                      "outNCI_first_200.props_molwriter.sdf";
  std::ostream *os = new std::ofstream(ofile.c_str());
  // std::ostream *os=new std::stringstream();
  auto *writer = new SDWriter(os, false);

  STR_VECT names;
#if 1
  ROMol *m1 = sdsup.next();
  MolOps::sanitizeMol(*(RWMol *)m1);
  delete m1;
#else
  ROMol *m1 = SmilesToMol("C1CC1");
  TEST_ASSERT(m1);
#endif
  sdsup.reset();
  int nDone = 0;
  while (!sdsup.atEnd()) {
    // std::cerr<<nDone<<std::endl;
    ROMol *mol = sdsup.next();
    // std::cerr<<"m:"<<mol<<std::endl;
    TEST_ASSERT(mol);
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    names.push_back(mname);
    // std::cerr<<"  w"<<std::endl;
    writer->write(*mol);

    // std::cerr<<"  ok"<<std::endl;

    delete mol;
    nDone++;
  }
  CHECK_INVARIANT(nDone == 200, "");
  writer->flush();
  CHECK_INVARIANT(writer->numMols() == 200, "");
  writer->close();

  delete writer;
  delete os;
#if 1
  // now read in the file we just finished writing
  SDMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == names[i], "");

    delete mol;
    i++;
  }
#endif
}

void testIssue3525000() {
  {
    std::string rdbase = getenv("RDBASE");
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525000.sdf";
    RWMol *mol = MolFileToMol(fname);
    TEST_ASSERT(mol);

    std::string cip;
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(8)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(8)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(9)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(9)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(10)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(10)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(14)->hasProp(common_properties::_CIPCode));
    // FIX: Marvin disagrees about this one:
    mol->getAtomWithIdx(14)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    std::string mb = MolToMolBlock(*mol);
    delete mol;
    mol = MolBlockToMol(mb);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(8)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(8)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(9)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(9)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(10)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(10)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(14)->hasProp(common_properties::_CIPCode));
    // FIX: Marvin disagrees about this one:
    mol->getAtomWithIdx(14)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(15)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete mol;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string fname =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525000b.sdf";
    RWMol *mol = MolFileToMol(fname);
    TEST_ASSERT(mol);
    MolOps::assignChiralTypesFrom3D(*mol);
    MolOps::assignStereochemistry(*mol);
    std::string cip;
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    std::string mb = MolToMolBlock(*mol);
    delete mol;
    mol = MolBlockToMol(mb);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete mol;
  }
}

void testIssue265() {
  {
    ROMol *m1 = SmilesToMol("C1ON1");
    TEST_ASSERT(m1);
    auto *conf = new Conformer(m1->getNumAtoms());
    RDGeom::Point3D p1(0, 0, 0);
    RDGeom::Point3D p2(1, 0, 0);
    RDGeom::Point3D p3(0, 1, 0);
    conf->setAtomPos(0, p1);
    conf->setAtomPos(1, p2);
    conf->setAtomPos(2, p3);
    m1->addConformer(conf);

    std::stringstream sstream;
    TDTWriter writer(&sstream);
    writer.write(*m1);
    std::string otext = sstream.str();
    TEST_ASSERT(otext == "$SMI<C1NO1>\n3D<0,0,0,0,1,0,1,0,0;>\n");
    writer.close();
    delete m1;
  }
}

void testMolFileChiralFlag() {
  {
    ROMol *m1 = SmilesToMol("C[C@H](Cl)F");
    TEST_ASSERT(m1);

    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1->hasProp(common_properties::_MolFileChiralFlag));
    TEST_ASSERT(m1->getProp<int>(common_properties::_MolFileChiralFlag) == 0);
    delete m1;
  }
  {
    ROMol *m1 = SmilesToMol("C[C@H](Cl)F");
    TEST_ASSERT(m1);
    m1->setProp(common_properties::_MolFileChiralFlag,
                static_cast<unsigned int>(1));
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1->hasProp(common_properties::_MolFileChiralFlag));
    delete m1;
  }
}

void testMolFileTotalValence() {
  BOOST_LOG(rdInfoLog) << "testing handling of mol file valence flags"
                       << std::endl;

  {
    RWMol *m1 = SmilesToMol("[Na]");
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 0);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m1;
  }
  {
    RWMol *m1 = SmilesToMol("[CH]");
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

    delete m1;
  }
  {
    RWMol *m1 = SmilesToMol("[CH2]");
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 2);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 2);
    delete m1;
  }
  {
    RWMol *m1 = SmilesToMol("[CH3]");
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileWithRxn() {
  BOOST_LOG(rdInfoLog) << "testing handling of mol files with reactions"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "rxn1.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 16);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->getProp<int>(common_properties::molRxnRole) == 1);
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnComponent) == 1);

    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(
        m->getAtomWithIdx(17)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnComponent) == 3);

    std::string mb = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 16);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->getProp<int>(common_properties::molRxnRole) == 1);
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnComponent) == 1);

    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(
        m->getAtomWithIdx(17)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnComponent) == 3);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSDWriterOptions() {
  BOOST_LOG(rdInfoLog) << "testing SDWriter options" << std::endl;
  {
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("V2000") != std::string::npos);
    TEST_ASSERT(txt.find("V3000") == std::string::npos);
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.setForceV3000(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("V2000") == std::string::npos);
    TEST_ASSERT(txt.find("V3000") != std::string::npos);
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.write(*mol);
    writer.setForceV3000(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("V2000") != std::string::npos);
    TEST_ASSERT(txt.find("V3000") != std::string::npos);
    TEST_ASSERT(txt.find("V2000") < txt.find("V3000"));
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.setForceV3000(true);
    writer.write(*mol);
    writer.setForceV3000(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("V2000") != std::string::npos);
    TEST_ASSERT(txt.find("V3000") != std::string::npos);
    TEST_ASSERT(txt.find("V2000") > txt.find("V3000"));
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("  1  2  2") != std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4") == std::string::npos);
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.setKekulize(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("  1  2  2") == std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4") != std::string::npos);
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.write(*mol);
    writer.setKekulize(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("  1  2  2") != std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4") != std::string::npos);
    TEST_ASSERT(txt.find("  1  2  2") < txt.find("  1  2  4"));
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss, false);
    RWMol *mol = SmilesToMol("c1ccccc1");
    writer.setKekulize(false);
    writer.write(*mol);
    writer.setKekulize(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");
    writer.close();

    std::string txt = ss.str();
    TEST_ASSERT(txt.find("  1  2  2") != std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4") != std::string::npos);
    TEST_ASSERT(txt.find("  1  2  2") > txt.find("  1  2  4"));
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testZBO() {
  BOOST_LOG(rdInfoLog) << "testing handling of ZBO specs" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName;
    fName = rdbase + "FeCO5.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType() == Bond::ZERO);

    std::string mb = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType() == Bond::ZERO);

    delete m;
  }

  {
    std::string fName;
    fName = rdbase + "H3BNH3.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getTotalNumHs() == 3);

    std::string mb = MolToMolBlock(*m);
    delete m;
    // std::cerr<<"MOLBLOCK:\n"<<mb<<"------\n";
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getTotalNumHs() == 3);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testV3000WriterDetails() {
  BOOST_LOG(rdInfoLog) << "testing details of v3000 writing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "chebi_57262.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 22);
    TEST_ASSERT(m->getNumBonds() == 21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope() == 2);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    TEST_ASSERT(mb.find("MASS=1") != std::string::npos);
    TEST_ASSERT(mb.find("MASS=2") != std::string::npos);

    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 22);
    TEST_ASSERT(m->getNumBonds() == 21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope() == 2);

    // repeat that one more time to make sure we're really solid:
    mb = MolToMolBlock(*m, true, -1, true, true);
    TEST_ASSERT(mb.find("MASS=1") != std::string::npos);
    TEST_ASSERT(mb.find("MASS=2") != std::string::npos);

    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 22);
    TEST_ASSERT(m->getNumBonds() == 21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope() == 2);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testV3000DoublePrecision() {
  BOOST_LOG(rdInfoLog) << "testing V3000 outputs coordinates at maximum robust double precision"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "precision.v3k.mol";
    RWMol *mol = MolFileToMol(fName);
    TEST_ASSERT(mol);
    size_t numAtoms = mol->getNumAtoms();
    TEST_ASSERT(numAtoms == 7);
    MolWriterParams params{true, true, true, 15};
    std::string molBlock = MolToMolBlock(*mol, params, -1);
    RWMol *readMol = MolBlockToMol(molBlock);
    TEST_ASSERT(numAtoms == readMol->getNumAtoms());
    const Conformer &conformer = mol->getConformer();
    const Conformer &readConformer = readMol->getConformer();
    for (size_t i = 0; i < numAtoms; i++) {
      std::cout << std::setprecision(15) << conformer.getAtomPos(i).x << ' ' << readConformer.getAtomPos(i).x << std::setprecision(6) << std::endl;
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).x - readConformer.getAtomPos(i).x) < 1e-15);
      std::cout << std::setprecision(15) << conformer.getAtomPos(i).y << ' ' << readConformer.getAtomPos(i).y << std::setprecision(6) << std::endl;
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).y - readConformer.getAtomPos(i).y) < 1e-15);
      std::cout << std::setprecision(15) << conformer.getAtomPos(i).z << ' ' << readConformer.getAtomPos(i).z << std::setprecision(6) << std::endl;
      TEST_ASSERT(std::abs(conformer.getAtomPos(i).z - readConformer.getAtomPos(i).z) < 1e-15);
    }
  }
}

void testGithub187() {
  BOOST_LOG(rdInfoLog) << "testing github issue 187: A not written to mol block"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "github187.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getNumBonds() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" A   0") != std::string::npos);

    // try the v3000 version:
    mb = MolToMolBlock(*m, true, -1, true, true);
    TEST_ASSERT(mb.find("V30 1 A 0") != std::string::npos);

    delete m;
  }

  {
    std::string fName = rdbase + "github187.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getNumBonds() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" A   0") != std::string::npos);

    // try the v3000 version:
    mb = MolToMolBlock(*m, true, -1, true, true);
    TEST_ASSERT(mb.find("V30 1 A 0") != std::string::npos);

    delete m;
  }

  {
    std::string fName = rdbase + "github187.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getNumBonds() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" Q   0") != std::string::npos);

    // try the v3000 version:
    mb = MolToMolBlock(*m, true, -1, true, true);
    TEST_ASSERT(mb.find("V30 1 \"NOT [C,H]\" 0") == std::string::npos);
    TEST_ASSERT(mb.find("1 Q 0") != std::string::npos);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub186() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 186: chiral S not written to ctab" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "github186.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getAtomWithIdx(6)->getChiralTag() != Atom::CHI_UNSPECIFIED &&
                m->getAtomWithIdx(6)->getChiralTag() != Atom::CHI_OTHER);

    std::string mb = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getAtomWithIdx(6)->getChiralTag() != Atom::CHI_UNSPECIFIED &&
                m->getAtomWithIdx(6)->getChiralTag() != Atom::CHI_OTHER);

    delete m;
  }
}

void testGithub189() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 189: Problems round-tripping Al2Cl6 via CTAB"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "github189.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getNumBonds() == 8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence() == 4);

    std::string mb = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getNumBonds() == 8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence() == 4);

    // try v3k
    mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getNumBonds() == 8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence() == 4);

    delete m;
  }
}

void testGithub266() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 266: Bond query information written to CTAB"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "bond-query.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrDoubleBond");

    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrDoubleBond");

    // try v3k
    mb = MolToMolBlock(*m, true, -1, true, true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrDoubleBond");

    delete m;
    delete m2;
  }

  {
    std::string fName = rdbase + "bond-query2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrAromaticBond");

    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrAromaticBond");

    // try v3k
    mb = MolToMolBlock(*m, true, -1, true, true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrAromaticBond");

    delete m;
    delete m2;
  }

  {
    std::string fName = rdbase + "bond-query3.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription() ==
                "DoubleOrAromaticBond");

    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "DoubleOrAromaticBond");

    // try v3k
    mb = MolToMolBlock(*m, true, -1, true, true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "DoubleOrAromaticBond");

    delete m;
    delete m2;
  }

  {
    ROMol *m = SmartsToMol("C-CN");
    TEST_ASSERT(m);
    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(!m2->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(!m2->getAtomWithIdx(2)->hasQuery());
    TEST_ASSERT(m2->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(m2->getAtomWithIdx(1)->getAtomicNum() == 6);
    TEST_ASSERT(m2->getAtomWithIdx(2)->getAtomicNum() == 7);
    TEST_ASSERT(!m2->getBondWithIdx(0)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "SingleOrAromaticBond");
    delete m;
    delete m2;
  }
}

void testGithub268() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 268: Bond topology information written to CTAB"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "bond-query4.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription() ==
                "BondAnd");
    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "BondAnd");

    // try v3k
    mb = MolToMolBlock(*m, true, -1, true, true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds() == 4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription() ==
                "BondAnd");

    delete m;
    delete m2;
  }
}

void testGithub357() {
  BOOST_LOG(rdInfoLog) << "testing github issue 357: Hydrogens in mol blocks "
                          "have a valence value set"
                       << std::endl;
  {
    ROMol *m1 = SmilesToMol("O");
    TEST_ASSERT(m1);
    ROMol *m2 = MolOps::addHs(*m1);
    TEST_ASSERT(m2);
    delete m1;
    std::string mb = MolToMolBlock(*m2);
    TEST_ASSERT(
        mb.find("    0.0000    0.0000    0.0000 H   0  0  0  0  0  1") ==
        std::string::npos);
    delete m2;
  }
}

void testNeedsUpdatePropertyCacheSDWriter() {
  BOOST_LOG(rdInfoLog)
      << "testing test needsUpdatePropertyCache functionality in SDwriter"
      << std::endl;
  {
    ROMol *m1 = SmilesToMol("c1ccccc1[NH]C(=O)", 0, false);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->needsUpdatePropertyCache() == true);
    std::string mb = MolToMolBlock(*m1);
    delete m1;
    ROMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2);
    delete m2;
  }
}

void testGithub488() {
  BOOST_LOG(rdInfoLog) << "testing github issue 488: SmilesWriter not creating "
                          "automatic name values for molecules read from CTABs"
                       << std::endl;
  {
    ROMol *m1 = SmilesToMol("O");
    TEST_ASSERT(m1);
    m1->setProp("_Name", "");
    std::stringstream ss;
    SmilesWriter w(&ss);
    w.write(*m1);
    m1->setProp("_Name", "foo");
    w.write(*m1);
    m1->clearProp("_Name");
    w.write(*m1);
    m1->setProp("_Name", " ");
    w.write(*m1);
    w.close();
    std::string txt = ss.str();
    TEST_ASSERT(txt.find("O 0") != std::string::npos);
    TEST_ASSERT(txt.find("O foo") != std::string::npos);
    TEST_ASSERT(txt.find("O 2") != std::string::npos);
    TEST_ASSERT(txt.find("O  \n") != std::string::npos);
    delete m1;
  }
}

void testGithub611() {
  BOOST_LOG(rdInfoLog) << "testing github issue 611: If wedged bonds are "
                          "already present, write them to mol blocks"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "Github611.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find("3  5  1  1") != std::string::npos);

    m->getBondWithIdx(2)->setBondDir(Bond::BEGINWEDGE);
    mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find("3  5  1  1") == std::string::npos);
    TEST_ASSERT(mb.find("3  4  1  1") != std::string::npos);
    delete m;
  }
}

void testGetSDText() {
  BOOST_LOG(rdInfoLog) << "testing SDWriter::getText()" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fname = rdbase + "NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);

    while (!sdsup.atEnd()) {
      ROMol *mol = sdsup.next();
      TEST_ASSERT(mol);
      std::string sdf = SDWriter::getText(*mol);
      SDMolSupplier tsupp;
      tsupp.setData(sdf);
      ROMol *mol2 = tsupp[0];
      TEST_ASSERT(mol2);
      std::string csmi1 = MolToSmiles(*mol, true);
      std::string csmi2 = MolToSmiles(*mol2, true);
      TEST_ASSERT(csmi1 == csmi2);
      STR_VECT pns = mol->getPropList(false, false);
      for (const auto &pn : pns) {
        TEST_ASSERT(mol2->hasProp(pn));
        TEST_ASSERT(mol->getProp<std::string>(pn) ==
                    mol2->getProp<std::string>(pn));
      }
      delete mol;
      delete mol2;
    }
  }
}

void testMolFileWriterDativeBonds() {
  BOOST_LOG(rdInfoLog) << "testing molfile writer dative bond support"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "dative_bonds_two.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondWithIdx(9)->getBondType() == Bond::DATIVE);

    std::string mb = MolToMolBlock(*m);
    // Bonds #9 and #10 (index 8 and 9) will have a bond type of 9
    // in the molfile.  Bond # --+ +-- Bond type.
    //                           | |
    TEST_ASSERT(mb.find("M  V30 9 9 5 11") != std::string::npos);
    TEST_ASSERT(mb.find("M  V30 10 9 10 11") != std::string::npos);

    // Roundtrip - can we read produced mol block above ?
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondWithIdx(9)->getBondType() == Bond::DATIVE);
    delete m;
    delete m2;
  }

  // Small molecules without dative bonds are output in V2000 format.
  {
    RWMol *m = SmilesToMol("CCC(=O)O[Cu]");
    TEST_ASSERT(m);
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find("0999 V2000") != std::string::npos);
    TEST_ASSERT(mb.find("0999 V3000") == std::string::npos);
    delete m;
  }
  // ... but molecules with dative bonds will always be
  // output in V3000 format.
  {
    RWMol *m = SmilesToMol("CCC(=O)O->[Cu]");
    TEST_ASSERT(m);
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find("0999 V2000") == std::string::npos);
    TEST_ASSERT(mb.find("0999 V3000") != std::string::npos);
    delete m;
  }
}

void testRGPMolFileWriterV2KV3K() {
  BOOST_LOG(rdInfoLog)
      << "testing that R groups do not lead either to M  ISO or MASS records"
      << std::endl;
  {
    auto m = R"CTAB(

     RDKit          2D
  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    2.4750    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    2.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    2.4750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
M  RGP  1   1   2
M  END)CTAB"_ctab;
    auto mbV2K = MolToMolBlock(*m);
    TEST_ASSERT(mbV2K.find("M  ISO") == std::string::npos);
    TEST_ASSERT(mbV2K.find("M  RGP") != std::string::npos);
    auto mbV3K = MolToV3KMolBlock(*m);
    TEST_ASSERT(mbV3K.find("MASS") == std::string::npos);
    TEST_ASSERT(mbV3K.find("RGROUPS") != std::string::npos);
  }
  {
    auto m = R"CTAB(
  MJ201100                      

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.5623    1.6625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2767    1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5623    2.4875    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  1  3  1  0  0  0  0
M  ISO  1   3   3
M  END)CTAB"_ctab;
    auto mbV2K = MolToMolBlock(*m);
    TEST_ASSERT(mbV2K.find("M  ISO") != std::string::npos);
    auto mbV3K = MolToV3KMolBlock(*m);
    TEST_ASSERT(mbV3K.find("MASS") != std::string::npos);
  }
}

int main() {
  RDLog::InitLogs();

#if 1
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriter()\n";
  testSmilesWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriter2()\n";
  testSmilesWriter2();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterNoNames()\n";
  testSmilesWriterNoNames();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterClose()\n";
  testSmilesWriterClose();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDWriter()\n";
  testSDWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testTDTWriter()\n";
  testTDTWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterStrm()\n";
  testSmilesWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDWriterStrm()\n";
  testSDWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testTDTWriterStrm()\n";
  testTDTWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDMemoryCorruption()\n";
  testSDMemoryCorruption();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testIssue265()\n";
  testIssue265();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testMolFileChiralFlag()\n";
  testMolFileChiralFlag();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testMolFileTotalValence();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testMolFileWithRxn();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testZBO();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";
  testSDWriterOptions();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testV3000WriterDetails();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testV3000DoublePrecision()\n";
  testV3000DoublePrecision();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub187();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub186();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub189();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub266();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub268();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub357();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testNeedsUpdatePropertyCacheSDWriter();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

#endif

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testIssue3525000()\n";
  testIssue3525000();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub488();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGithub611();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testGetSDText();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testMolFileWriterDativeBonds();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n";
  testRGPMolFileWriterV2KV3K();
  BOOST_LOG(rdInfoLog) << "-----------------------------------------\n\n";
}
