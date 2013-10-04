//  $Id$
// 
//   Copyright (C) 2002-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname = rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";
  
  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("_Name", mname);
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
  writer->flush();
  delete writer;
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname);
  int i = 0;
  mol = nSup->next();
  while (mol){
    std::string mname, pval;
    mol->getProp("_Name", mname);
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
}

void testSmilesWriter2() {
  {
    std::stringstream ss;
    SmilesWriter *writer = new SmilesWriter(&ss, " ","Name",false);
    RWMol *mol;

    mol = SmilesToMol("c1ccccc1");
    //MolOps::Kekulize(*mol);
    writer->write(*mol);
    delete mol;

    mol = SmilesToMol("F[C@H](Cl)Br");
    writer->write(*mol);
    delete mol;
    writer->flush();
    TEST_ASSERT(ss.str()=="c1ccccc1 0\nFC(Cl)Br 1\n");
    delete writer;
  }
  {
    std::stringstream ss;
    SmilesWriter *writer = new SmilesWriter(&ss, " ","Name",false,false,true);
    RWMol *mol;

    mol = SmilesToMol("c1ccccc1");
    MolOps::Kekulize(*mol);
    writer->write(*mol);
    delete mol;

    mol = SmilesToMol("F[C@H](Cl)Br");
    writer->write(*mol);
    delete mol;
    writer->flush();
    TEST_ASSERT(ss.str()=="C1=CC=CC=C1 0\nF[C@H](Cl)Br 1\n");
    delete writer;
  }
}

void testSmilesWriterNoNames() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname = rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";
  
  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname," ","");
  writer->setProps(propNames);

  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("Column_2", pval);
    mol->setProp("_Name","bogus");
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
  delete writer;
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname,",",0,-1);
  int i = 0;
  mol = nSup->next();
  while (mol){
    std::string mname, pval;
    mol->getProp("_Name", mname);
    TEST_ASSERT(mname!="bogus");
    mol->getProp("Column_2", pval);
    TEST_ASSERT(pval == props[i]);
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
}

void testSmilesWriterClose() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname = rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";
  
  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname," ","");
  writer->setProps(propNames);

  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("Column_2", pval);
    mol->setProp("_Name","bogus");
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
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname,",",0,-1);
  int i = 0;
  mol = nSup->next();
  while (mol){
    std::string mname, pval;
    mol->getProp("_Name", mname);
    TEST_ASSERT(mname!="bogus");
    mol->getProp("Column_2", pval);
    TEST_ASSERT(pval == props[i]);
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
}

void testSDWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);
  
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few.sdf";
  SDWriter *writer = new SDWriter(ofile);

  STR_VECT names;
 
  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp("_Name", mname);
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
    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == names[i], "");
    
    delete mol;
    i++;
  }

  // now read in a file with aromatic information on the bonds
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
  SDMolSupplier nreader(infile);
  i = 0;
  while (!nreader.atEnd()) {
    ROMol *mol = nreader.next();
    std::string mname;
    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == names[i], "");
    i++;
    
    delete mol;
  }

}

void testTDTWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);
  
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few.tdt";
  TDTWriter *writer = new TDTWriter(ofile);

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
  delete writer;

  // now read in the file we just finished writing
  TDTMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    if(mol){
      std::string mname;
      mol->getProp("CAS_RN", mname);
      CHECK_INVARIANT(mname == names[i], "");
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i==16);
}

void testSmilesWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname = rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";
  std::ofstream *oStream=new std::ofstream(oname.c_str());
  
  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oStream, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp("_Name", mname);
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
  writer->flush();
  delete writer;
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname);
  int i = 0;
  mol = nSup->next();
  while (mol){
    std::string mname, pval;
    mol->getProp("_Name", mname);
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

}


void testSDWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  {
    std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);
  
    std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few.sdf";
    std::ofstream *oStream=new std::ofstream(ofile.c_str());

    SDWriter *writer = new SDWriter(oStream);

    STR_VECT names;
 
    while (!sdsup.atEnd()) {
      ROMol *mol = sdsup.next();
      std::string mname;
      mol->getProp("_Name", mname);
      names.push_back(mname);

      writer->write(*mol);
      delete mol;
    }
    writer->flush();
    CHECK_INVARIANT(writer->numMols() == 16, "");
    delete writer;

    // now read in the file we just finished writing
    SDMolSupplier reader(ofile);
    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      std::string mname;
      mol->getProp("_Name", mname);
      CHECK_INVARIANT(mname == names[i], "");
    
      delete mol;
      i++;
    }
  }

  {
    // now read in a file with aromatic information on the bonds
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
    SDMolSupplier nreader(infile);
    unsigned int i = 0;
    while (!nreader.atEnd()) {
      ROMol *mol = nreader.next();
      TEST_ASSERT(mol);
      ++i;
      delete mol;
    }
    TEST_ASSERT(i==16);
  }
}

void testTDTWriterStrm() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);
  
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_few.tdt";
  std::ofstream *oStream=new std::ofstream(ofile.c_str());
  TDTWriter *writer = new TDTWriter(oStream);

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
  delete writer;

  // now read in the file we just finished writing
  TDTMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    if(mol){
      std::string mname;
      mol->getProp("CAS_RN", mname);
      CHECK_INVARIANT(mname == names[i], "");
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i==16);
}

void testSDMemoryCorruption() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Data/NCI/first_200.props.sdf";
  SDMolSupplier sdsup(fname,true);
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_first_200.props.sdf";
  std::ostream *os=new std::ofstream(ofile.c_str());
  //std::ostream *os=new std::stringstream();
  SDWriter *writer = new SDWriter(os,false);

  STR_VECT names;
#if 1
  ROMol *m1=sdsup.next();
  MolOps::sanitizeMol(*(RWMol *)m1);
#else
  ROMol *m1=SmilesToMol("C1CC1");
  TEST_ASSERT(m1);
#endif
  sdsup.reset();
  int nDone=0;
  while (!sdsup.atEnd()) {
    //std::cerr<<nDone<<std::endl;
    ROMol *mol = sdsup.next();
    //std::cerr<<"m:"<<mol<<std::endl;
    TEST_ASSERT(mol);
    std::string mname;
    mol->getProp("_Name", mname);
    names.push_back(mname);
    //std::cerr<<"  w"<<std::endl;
    writer->write(*mol);

    //std::cerr<<"  ok"<<std::endl;
    
    delete mol;
    nDone++;
  }
  CHECK_INVARIANT(nDone == 200, "");
  writer->flush();
  CHECK_INVARIANT(writer->numMols() == 200, "");

  delete writer;
#if 1
  // now read in the file we just finished writing
  SDMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == names[i], "");

    
    delete mol;
    i++;
  }
#endif
}


void testIssue3525000() {
  {
    std::string rdbase = getenv("RDBASE");
    std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525000.sdf";
    RWMol *mol = MolFileToMol(fname);
    TEST_ASSERT(mol);
    std::string cip;
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(6)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(8)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(8)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(9)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(9)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(10)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(10)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(14)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(14)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(15)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");

    std::string mb=MolToMolBlock(*mol);
    delete mol;
    mol = MolBlockToMol(mb);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(6)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(6)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(8)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(8)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(9)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(9)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(10)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(10)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(14)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(14)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(15)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(15)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
  }  
  {
    std::string rdbase = getenv("RDBASE");
    std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525000b.sdf";
    RWMol *mol = MolFileToMol(fname);
    TEST_ASSERT(mol);
    MolOps::assignChiralTypesFrom3D(*mol);
    MolOps::assignStereochemistry(*mol);
    std::string cip;
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");

    std::string mb=MolToMolBlock(*mol);
    delete mol;
    mol = MolBlockToMol(mb);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(0)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(1)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(2)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(3)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="R");
    TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp("_CIPCode"));
    mol->getAtomWithIdx(4)->getProp("_CIPCode",cip);
    TEST_ASSERT(cip=="S");
  }  
}

void testIssue265() {

  {
    ROMol *m1=SmilesToMol("C1ON1");
    TEST_ASSERT(m1);
    Conformer *conf=new Conformer(m1->getNumAtoms());
    RDGeom::Point3D p1(0,0,0);
    RDGeom::Point3D p2(1,0,0);
    RDGeom::Point3D p3(0,1,0);
    conf->setAtomPos(0,p1);
    conf->setAtomPos(1,p2);
    conf->setAtomPos(2,p3);
    m1->addConformer(conf);

    std::stringstream sstream;
    TDTWriter writer(&sstream);
    writer.write(*m1);
    writer.flush();
    std::string otext=sstream.str();
    TEST_ASSERT(otext=="$SMI<C1NO1>\n3D<0,0,0,0,1,0,1,0,0;>\n");
  }
}

void testMolFileChiralFlag() {

  {
    ROMol *m1=SmilesToMol("C[C@H](Cl)F");
    TEST_ASSERT(m1);

    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(!m1->hasProp("_MolFileChiralFlag"));
  }
  {
    ROMol *m1=SmilesToMol("C[C@H](Cl)F");
    TEST_ASSERT(m1);
    m1->setProp("_MolFileChiralFlag",static_cast<unsigned int>(1));
    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1->hasProp("_MolFileChiralFlag"));
  }
}

void testMolFileTotalValence(){
  BOOST_LOG(rdInfoLog) << "testing handling of mol file valence flags" << std::endl;

  {
    RWMol *m1=SmilesToMol("[Na]");
    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms()==1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs()==0);    
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m1;
  }
  {
    RWMol *m1=SmilesToMol("[CH]");
    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms()==1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs()==1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons()==1);

    delete m1;
  }
  {
    RWMol *m1=SmilesToMol("[CH2]");
    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms()==1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs()==2);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons()==2);
    delete m1;
  }
  {
    RWMol *m1=SmilesToMol("[CH3]");
    std::string mb=MolToMolBlock(*m1);
    delete m1;
    m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms()==1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs()==3);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons()==1);
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


int main() {
  RDLog::InitLogs();
#if 1
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriter()\n";
  testSmilesWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriter2()\n";
  testSmilesWriter2();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterNoNames()\n";
  testSmilesWriterNoNames();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterClose()\n";
  testSmilesWriterNoNames();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDWriter()\n";
  testSDWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testTDTWriter()\n";
  testTDTWriter();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSmilesWriterStrm()\n";
  testSmilesWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDWriterStrm()\n";
  testSDWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testTDTWriterStrm()\n";
  testTDTWriterStrm();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testSDMemoryCorruption()\n";
  testSDMemoryCorruption();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testIssue3525000()\n";
  testIssue3525000();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testIssue265()\n";
  testIssue265();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "Running testMolFileChiralFlag()\n";
  testMolFileChiralFlag();
  BOOST_LOG(rdInfoLog) << "Finished\n";
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
#endif

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testMolFileTotalValence();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
}
