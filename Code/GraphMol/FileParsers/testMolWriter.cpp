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

void testMolFileWithRxn(){
  BOOST_LOG(rdInfoLog) << "testing handling of mol files with reactions" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase+"rxn1.mol";
    ROMol *m=MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==18);
    TEST_ASSERT(m->getNumBonds()==16);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("molRxnRole"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>("molRxnRole")==1);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("molRxnComponent"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>("molRxnComponent")==1);

    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp("molRxnRole"));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>("molRxnRole")==2);
    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp("molRxnComponent"));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>("molRxnComponent")==3);

    std::string mb=MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==18);
    TEST_ASSERT(m->getNumBonds()==16);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("molRxnRole"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>("molRxnRole")==1);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("molRxnComponent"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>("molRxnComponent")==1);

    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp("molRxnRole"));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>("molRxnRole")==2);
    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp("molRxnComponent"));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>("molRxnComponent")==3);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSDWriterOptions() {
  BOOST_LOG(rdInfoLog) << "testing SDWriter options" << std::endl;
  {
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("V2000")!=std::string::npos);
    TEST_ASSERT(txt.find("V3000")==std::string::npos);
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.setForceV3000(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("V2000")==std::string::npos);
    TEST_ASSERT(txt.find("V3000")!=std::string::npos);
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.write(*mol);
    writer.setForceV3000(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("V2000")!=std::string::npos);
    TEST_ASSERT(txt.find("V3000")!=std::string::npos);
    TEST_ASSERT(txt.find("V2000")<txt.find("V3000"));
  }
  {
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.setForceV3000(true);
    writer.write(*mol);
    writer.setForceV3000(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("V2000")!=std::string::npos);
    TEST_ASSERT(txt.find("V3000")!=std::string::npos);
    TEST_ASSERT(txt.find("V2000")>txt.find("V3000"));
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("  1  2  2")!=std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4")==std::string::npos);
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.setKekulize(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 1, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("  1  2  2")==std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4")!=std::string::npos);
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.write(*mol);
    writer.setKekulize(false);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("  1  2  2")!=std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4")!=std::string::npos);
    TEST_ASSERT(txt.find("  1  2  2")<txt.find("  1  2  4"));
  }
  {
    // kekulization
    std::stringstream ss;
    SDWriter writer(&ss,false);
    RWMol *mol=SmilesToMol("c1ccccc1");
    writer.setKekulize(false);
    writer.write(*mol);
    writer.setKekulize(true);
    writer.write(*mol);
    delete mol;
    writer.flush();
    CHECK_INVARIANT(writer.numMols() == 2, "");

    std::string txt=ss.str();
    TEST_ASSERT(txt.find("  1  2  2")!=std::string::npos);
    TEST_ASSERT(txt.find("  1  2  4")!=std::string::npos);
    TEST_ASSERT(txt.find("  1  2  2")>txt.find("  1  2  4"));
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testZBO(){
  BOOST_LOG(rdInfoLog) << "testing handling of ZBO specs" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName;
    fName = rdbase+"FeCO5.mol";
    ROMol *m=MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);
    TEST_ASSERT(m->getNumBonds()==10);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType()==Bond::ZERO);

    std::string mb=MolToMolBlock(*m);
    delete m;
    std::cerr<<"MOLBLOCK:\n"<<mb<<"------\n";
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);
    TEST_ASSERT(m->getNumBonds()==10);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType()==Bond::ZERO);

    delete m;
  }

  {
    std::string fName;
    fName = rdbase+"H3BNH3.mol";
    ROMol *m=MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==2);
    TEST_ASSERT(m->getNumBonds()==1);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs()==3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs()==3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getTotalNumHs()==3);

    std::string mb=MolToMolBlock(*m);
    delete m;
    //std::cerr<<"MOLBLOCK:\n"<<mb<<"------\n";
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==2);
    TEST_ASSERT(m->getNumBonds()==1);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType()==Bond::ZERO);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs()==3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs()==3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs()==3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getTotalNumHs()==3);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testV3000WriterDetails(){
  BOOST_LOG(rdInfoLog) << "testing details of v3000 writing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "chebi_57262.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==22);
    TEST_ASSERT(m->getNumBonds()==21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope()==1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope()==2);

    std::string mb=MolToMolBlock(*m,true,-1,true,true);
    TEST_ASSERT(mb.find("MASS=1")!=std::string::npos);
    TEST_ASSERT(mb.find("MASS=2")!=std::string::npos);

    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==22);
    TEST_ASSERT(m->getNumBonds()==21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope()==1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope()==2);

    // repeat that one more time to make sure we're really solid:
    mb=MolToMolBlock(*m,true,-1,true,true);
    TEST_ASSERT(mb.find("MASS=1")!=std::string::npos);
    TEST_ASSERT(mb.find("MASS=2")!=std::string::npos);

    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==22);
    TEST_ASSERT(m->getNumBonds()==21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope()==1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum()==0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope()==2);
    
    delete m;
  }
  
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub187(){
  BOOST_LOG(rdInfoLog) << "testing github issue 187: A not written to mol block" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  
  {
    std::string fName = rdbase + "github187.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==1);
    TEST_ASSERT(m->getNumBonds()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb=MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" A   0")!=std::string::npos);

    // try the v3000 version:
    mb=MolToMolBlock(*m,true,-1,true,true);
    TEST_ASSERT(mb.find("V30 1 A 0")!=std::string::npos);
    
    delete m;
  }
  
  {
    std::string fName = rdbase + "github187.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==1);
    TEST_ASSERT(m->getNumBonds()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb=MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" A   0")!=std::string::npos);

    // try the v3000 version:
    mb=MolToMolBlock(*m,true,-1,true,true);
    TEST_ASSERT(mb.find("V30 1 A 0")!=std::string::npos);
    
    delete m;
  }

  {
    std::string fName = rdbase + "github187.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==1);
    TEST_ASSERT(m->getNumBonds()==0);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    std::string mb=MolToMolBlock(*m);
    TEST_ASSERT(mb.find(" Q   0")!=std::string::npos);

    // try the v3000 version:
    mb=MolToMolBlock(*m,true,-1,true,true);
    TEST_ASSERT(mb.find("V30 1 \"NOT [C,H]\" 0")!=std::string::npos);
    
    delete m;
  }
  
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub186(){
  BOOST_LOG(rdInfoLog) << "testing github issue 186: chiral S not written to ctab" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  
  {
    std::string fName = rdbase + "github186.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);
    TEST_ASSERT(m->getNumBonds()==10);
    TEST_ASSERT(m->getAtomWithIdx(6)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                m->getAtomWithIdx(6)->getChiralTag()!=Atom::CHI_OTHER
                );

    std::string mb=MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==11);
    TEST_ASSERT(m->getNumBonds()==10);
    TEST_ASSERT(m->getAtomWithIdx(6)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                m->getAtomWithIdx(6)->getChiralTag()!=Atom::CHI_OTHER
                );

    delete m;
  }
  
}

void testGithub189(){
  BOOST_LOG(rdInfoLog) << "testing github issue 189: Problems round-tripping Al2Cl6 via CTAB" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  
  {
    std::string fName = rdbase + "github189.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    TEST_ASSERT(m->getNumBonds()==8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence()==4);


    std::string mb=MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    TEST_ASSERT(m->getNumBonds()==8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence()==4);

    // try v3k
    mb=MolToMolBlock(*m,true,-1,true,true);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==8);
    TEST_ASSERT(m->getNumBonds()==8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(2)->getTotalValence()==4);

    delete m;
  }
  
}

void testGithub266(){
  BOOST_LOG(rdInfoLog) << "testing github issue 266: Bond query information written to CTAB" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  
  {
    std::string fName = rdbase + "bond-query.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds()==4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    std::string mb=MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    // try v3k
    mb=MolToMolBlock(*m,true,-1,true,true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");

    delete m;
  }
  
  {
    std::string fName = rdbase + "bond-query2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds()==4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    std::string mb=MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    // try v3k
    mb=MolToMolBlock(*m,true,-1,true,true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");

    delete m;
  }
  
  {
    std::string fName = rdbase + "bond-query3.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds()==4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    std::string mb=MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");
    
    // try v3k
    mb=MolToMolBlock(*m,true,-1,true,true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");

    delete m;
  }

  {
    ROMol *m=SmartsToMol("C-CN");
    TEST_ASSERT(m);
    std::string mb=MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(!m2->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(!m2->getAtomWithIdx(2)->hasQuery());
    TEST_ASSERT(m2->getAtomWithIdx(0)->getAtomicNum()==6);
    TEST_ASSERT(m2->getAtomWithIdx(1)->getAtomicNum()==6);
    TEST_ASSERT(m2->getAtomWithIdx(2)->getAtomicNum()==7);
    TEST_ASSERT(!m2->getBondWithIdx(0)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondOr");

  }

  
}


void testGithub268(){
  BOOST_LOG(rdInfoLog) << "testing github issue 268: Bond topology information written to CTAB" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  
  {
    std::string fName = rdbase + "bond-query4.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds()==4);
    TEST_ASSERT(m->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getBondWithIdx(1)->getQuery()->getDescription()=="BondAnd");
    
    std::string mb=MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondAnd");
    
    // try v3k
    mb=MolToMolBlock(*m,true,-1,true,true);
    delete m2;
    m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2->getNumBonds()==4);
    TEST_ASSERT(m2->getBondWithIdx(1)->hasQuery());
    TEST_ASSERT(m2->getBondWithIdx(1)->getQuery()->getDescription()=="BondAnd");

    delete m;
  }
}  

void testGithub357(){
  BOOST_LOG(rdInfoLog) << "testing github issue 357: Hydrogens in mol blocks have a valence value set" << std::endl;
  {
    ROMol *m1=SmilesToMol("O");
    TEST_ASSERT(m1);
    ROMol *m2=MolOps::addHs(*m1);
    TEST_ASSERT(m2);
    delete m1;
    std::string mb=MolToMolBlock(*m2);
    TEST_ASSERT(mb.find("    0.0000    0.0000    0.0000 H   0  0  0  0  0  1")==std::string::npos);
  }
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
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testMolFileWithRxn();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testZBO();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  testSDWriterOptions();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testV3000WriterDetails();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub187();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub186();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub189();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub266();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub268();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";

  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n";
  testGithub357();
  BOOST_LOG(rdInfoLog) <<  "-----------------------------------------\n\n";
  
}
