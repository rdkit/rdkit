//  $Id$
// 
//   Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <fstream>

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


int main() {
  RDLog::InitLogs();

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testSmilesWriter()\n";
  testSmilesWriter();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testSDWriter()\n";
  testSDWriter();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testTDTWriter()\n";
  testTDTWriter();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testSmilesWriterStrm()\n";
  testSmilesWriterStrm();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testSDWriterStrm()\n";
  testSDWriterStrm();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testTDTWriterStrm()\n";
  testTDTWriterStrm();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";

  std::cout <<  "-----------------------------------------\n";
  std::cout << "Running testSDMemoryCorruption()\n";
  testSDMemoryCorruption();
  std::cout << "Finished\n";
  std::cout <<  "-----------------------------------------\n\n";
}
