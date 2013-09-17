//  $Id$
// 
//   Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
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
#include <map>

#include "MolSupplier.h"
#include "MolWriters.h"
#include "FileParsers.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>

#ifdef TEST_GZIP_SD
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
namespace io=boost::iostreams;
#endif

using namespace RDKit;

int testMolSup() {
  
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  
  {
    SDMolSupplier sdsup(fname);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i==16);
  }
  {
    SDMolSupplier sdsup(fname);
    for(unsigned int i=0;i<16;++i){
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
    }
    // test issue 3524949:
    TEST_ASSERT(sdsup.atEnd());
    bool ok=false; 
    try{
      ROMol *mol = sdsup.next();
    } catch (FileParseException &) {
      ok=true;
    }
    TEST_ASSERT(ok);
  }
  {
    std::ifstream strm(fname.c_str());
    SDMolSupplier sdsup(&strm,false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i==16);
  }
  {
    std::ifstream *strm=new std::ifstream(fname.c_str());
    SDMolSupplier sdsup(strm,true);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i==16);
  }
  return 1;
}

void testRandMolSup() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  //std::string fname("../test_data/NCI_aids_few.sdf");
  SDMolSupplier sdsup(fname);
  
  ROMol *tmol = sdsup[7];
  delete tmol;

  CHECK_INVARIANT(sdsup.length() == 16, "");

  STR_VECT names;
  names.push_back(std::string("48"));
  names.push_back(std::string("128"));
  names.push_back(std::string("164"));
  names.push_back(std::string("180"));
  names.push_back(std::string("192"));
  names.push_back(std::string("210"));
  names.push_back(std::string("213"));
  names.push_back(std::string("229"));

  int i;
  for (i = 0; i < 8; i++) {
    ROMol *mol = sdsup[2*i];
    std::string mname;
    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == names[i], "");
    delete mol;
  }
  
  // get a random molecule
  ROMol *mol = sdsup[5];
  TEST_ASSERT(mol);
  std::string mname;
  mol->getProp("_Name", mname);
  delete mol;
  CHECK_INVARIANT(mname == "170", "");

  // get the last molecule:
  mol = sdsup[15];
  TEST_ASSERT(mol);
  delete mol;

  // and make sure we're at the end:
  TEST_ASSERT(sdsup.atEnd());
  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  mol = sdsup[0];
  TEST_ASSERT(mol);
  delete mol;


  // Issue 113: calling length before grabbing a molecule results in crashes:
  SDMolSupplier sdsup2(fname);
  CHECK_INVARIANT(sdsup2.length() == 16, "");
  
}

void testSmilesSup() {
  std::string mname;
  std::string fname;
  ROMol *mol;

  std::string rdbase = getenv("RDBASE");
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);
    TEST_ASSERT(nSup2.length() == 10);
  }
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);

    mol = nSup2[3];
    TEST_ASSERT(!nSup2.atEnd())
    TEST_ASSERT(nSup2.length() == 10);

    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == "4", "");
    mol->getProp("TPSA", mname);
    CHECK_INVARIANT(mname == "82.78", "");
    delete mol;

    mol = nSup2[9];
    TEST_ASSERT(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    TEST_ASSERT(mol);
    delete mol;
  }
  {
    std::ifstream strm(fname.c_str(),std::ios_base::binary);
    SmilesMolSupplier nSup2(&strm,false,",", 1, 0, true);

    mol = nSup2[3];
    CHECK_INVARIANT(nSup2.length() == 10, "");

    mol->getProp("_Name", mname);
    CHECK_INVARIANT(mname == "4", "");
    mol->getProp("TPSA", mname);
    CHECK_INVARIANT(mname == "82.78", "");
    delete mol;

    mol = nSup2[9];
    TEST_ASSERT(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    TEST_ASSERT(mol);
    delete mol;
  }


  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  SmilesMolSupplier smiSup(fname, ",", 0, -1);
  
  mol = smiSup[16];
  
  mol->getProp("TPSA", mname);
  CHECK_INVARIANT(mname == "46.25", "");
  delete mol;

  mol = smiSup[8];
  mol->getProp("TPSA", mname);
  CHECK_INVARIANT(mname == "65.18", "");
  delete mol;

  int len = smiSup.length();
  CHECK_INVARIANT(len == 200, "");

  smiSup.reset();
  int i = 0;
  mol = smiSup.next();
  while (1) {
    std::string mname;
    mol->getProp("_Name", mname);
    i++;
    delete mol;
    try {
      mol = smiSup.next();
    } catch (FileParseException &) {
      break;
    }
  }
  
  CHECK_INVARIANT(i == 200, "");

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);

  // check the length before we read anything out...
  //  this was a problem at one point (Issue 113)
  CHECK_INVARIANT(nSup->length() == 10, "");
  mol = (*nSup)[3];
  

  mol->getProp("_Name", mname);
  CHECK_INVARIANT(mname == "4", "");
  mol->getProp("Column_2", mname);
  CHECK_INVARIANT(mname == "82.78", "");

  delete nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  unsigned int nRead=0;
  while(!nSup->atEnd()){
    mol=nSup->next();
    TEST_ASSERT(mol);
    delete mol;
    nRead++;
  }
  TEST_ASSERT(nSup->length() == 10);
  TEST_ASSERT(nRead == 10);

  delete nSup;

}

void testSmilesSupFromText() {
  std::string mname;
  std::string fname;
  ROMol *mol;

  SmilesMolSupplier nSup2;
  std::string text;
  bool failed;
  int nAts;

  // this was a delightful boundary condition:
  BOOST_LOG(rdErrorLog) << "------------------------------------------------------" << std::endl;
  text="CC\n"
    "CCC\n"
    "CCOC\n"
    "CCCCOC"
    ;
  {
    nSup2.setData(text," ",0,-1,false,true);
    //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
    mol = nSup2.next();
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==2);
    mol = nSup2[3];
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==6);
    TEST_ASSERT(nSup2.length()==4);
    failed=false;
    try {
      mol = nSup2[4];
    } catch (FileParseException &) {
      failed=true;
    }
    TEST_ASSERT(failed);
    mol = nSup2[2];
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==4);
    TEST_ASSERT(mol->hasProp("_Name"));
    mol->getProp("_Name",mname);
    TEST_ASSERT(mname=="2");
  }
  {
    nSup2.setData(text," ",0,-1,false,true);
    mol = nSup2[2];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==4);
    TEST_ASSERT(mol->hasProp("_Name"));
    mol->getProp("_Name",mname);
    TEST_ASSERT(mname=="2");
    mol = nSup2[3];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==6);
    TEST_ASSERT(mol->hasProp("_Name"));
    mol->getProp("_Name",mname);
    TEST_ASSERT(mname=="3");
  }
  {
    nSup2.setData(text," ",0,-1,false,true);
    mol = nSup2[3];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==6);
    TEST_ASSERT(mol->hasProp("_Name"));
    mol->getProp("_Name",mname);
    TEST_ASSERT(mname=="3");
    mol = nSup2[2];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts==4);
    TEST_ASSERT(mol->hasProp("_Name"));
    mol->getProp("_Name",mname);
    TEST_ASSERT(mname=="2");
  }
  // --------------
  // basics:
  text="Id SMILES Column_2\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-3 CCC 9.0\n"
    "mol-4 CCCC 16.0\n"
    ;
#if 1
  nSup2.setData(text," ",1,0,true,true);
  mol = nSup2[3];
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  CHECK_INVARIANT(nSup2.length() == 4, "");
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-4");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="16.0");

  // ensure that we can call setData a second time:
  text="Id SMILES Column_2\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-3 CCC 9.0\n"
    ;
  nSup2.setData(text," ",1,0,true,true);
  CHECK_INVARIANT(nSup2.length() == 3, "");
  mol = nSup2[2];
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-3");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="9.0");


  // now test for failure handling:
  text="Id SMILES Column_2\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-3 fail 9.0\n"
    "mol-4 CCCC 16.0\n"
    ;
  nSup2.setData(text," ",1,0,true,true);
  mol = nSup2[3];
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 4);
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-4");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="16.0");
  // failures should give null molecules:
  mol = nSup2[2];
  TEST_ASSERT(!mol);
#endif
  // issue 114, no \n at EOF:
  text="Id SMILES Column_2\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-4 CCCC 16.0\n"
    ;
  nSup2.setData(text," ",1,0,true,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-4");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="16.0");
  TEST_ASSERT(nSup2.atEnd());

  text="Id SMILES Column_2\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-4 CCCC 16.0"
    ;
  nSup2.setData(text," ",1,0,true,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-4");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="16.0");
  TEST_ASSERT(nSup2.atEnd());

  try {
    mol = nSup2[3];
  } catch (FileParseException &) {
    failed=true;
  }
  TEST_ASSERT(failed);

  text="mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-4 CCCC 16.0"
    ;
  nSup2.setData(text," ",1,0,false,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-4");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="16.0");

  text="C\n"
    "CC\n"
    "CCCC"
    ;
  nSup2.setData(text," ",0,-1,false,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==4);


  // this was a delightful boundary condition:
  BOOST_LOG(rdErrorLog) << "------------------------------------------------------" << std::endl;
  text="CC\n"
    "CCC\n"
    "CCOC\n"
    "CCCCOC"
    ;
  nSup2.setData(text," ",0,-1,false,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  nSup2.next();
  mol = nSup2[3];
  TEST_ASSERT(nSup2.length()==4);
  failed=false;
  try {
    mol = nSup2[4];
  } catch (FileParseException &) {
    failed=true;
  }
  TEST_ASSERT(failed);

  BOOST_LOG(rdErrorLog) << "------------------------------------------------------" << std::endl;
  // this was a delightful boundary condition:
  text="CC\n"
    "CCC\n"
    "CCOC\n"
    "CCCCOC"
    ;
  nSup2.setData(text," ",0,-1,false,true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  failed=false;
  try {
    mol = nSup2[4];
  } catch (FileParseException &) {
    failed=true;
  }
  TEST_ASSERT(failed);
  BOOST_LOG(rdErrorLog) << ">>> This may result in an infinite loop.  It should finish almost immediately:" << std::endl;
  TEST_ASSERT(nSup2.length()==4);
  BOOST_LOG(rdErrorLog) << "<<< done." << std::endl;

  nSup2.reset();
  unsigned int nDone=0;
  while(!nSup2.atEnd()){
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  TEST_ASSERT(nDone==nSup2.length());

  
  // ensure that we can call setData a second time:
  text="Id SMILES Column_2\n"
    "# comment, ignore\n"
    "mol-1 C 1.0\n"
    "mol-2 CC 4.0\n"
    "mol-3 CCC 9.0\n"
    "mol-4 CCCC 16.0\n"
    ;
  nSup2.setData(text," ",1,0,true,true);
  mol = nSup2[2];
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-3");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="9.0");
  mol = nSup2[1];
  mol->getProp("_Name",mname);
  TEST_ASSERT(mname=="mol-2");
  mol->getProp("Column_2",mname);
  TEST_ASSERT(mname=="4.0");

  // this was a delightful boundary condition:
  text="CC\n"
    "CCC\n"
    "CCOC\n"
    "CCCCOC\n"
    "\n"
    "\n"
    ;
  nSup2.setData(text," ",0,-1,false,true);
  TEST_ASSERT(nSup2.length()==4);
  nSup2.reset();
  nDone=0;
  while(!nSup2.atEnd()){
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  TEST_ASSERT(nDone==nSup2.length());

};

void testSmilesWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  //std::string fname = "../test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname = rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles.csv";
  //std::string oname = "../test_data/outSmiles.csv";
  
  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  //BOOST_LOG(rdErrorLog) << "WRITING" << std::endl;
  while (mol) {
    //BOOST_LOG(rdErrorLog) << "MOL: " << MolToSmiles(*mol) << std::endl;
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
    BOOST_LOG(rdInfoLog) << mname << "\n";
    //CHECK_INVARIANT(mname == names[i], "");
    
    delete mol;
    i++;
  }
  BOOST_LOG(rdInfoLog) << i << "\n";
  /*
  // now read in a file with aromatic information on the bonds
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
  SDMolSupplier nreader(infile);
  i = 0;
  while (!nreader.atEnd()) {
    ROMol *mol = nreader.next();
    std::string mname;
    mol->getProp("_Name", mname);
    BOOST_LOG(rdInfoLog) << mname << "\n";
    //CHECK_INVARIANT(mname == names[i], "");
    i++;
    
    delete mol;
    }*/

}

void testSDSupplierEnding() {
  std::string rdbase = getenv("RDBASE");
  // test the SD supplier to check if it properly handle the end of sd file 
  // conditions
  // should work fine if the sd file end with  a $$$$ follwed by blank line or no 
  // no blank lines
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/esters_end.sdf";
  int i = 0;
  SDMolSupplier reader(infile);
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp("_Name", mname);
    i++;
    delete mol;
  }
  CHECK_INVARIANT(i == 6, "");
}

void testSuppliersEmptyFile() {
  std::string rdbase = getenv("RDBASE");
  {
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/empty.sdf";
    SDMolSupplier reader(infile);
    TEST_ASSERT(reader.atEnd());
  }
  {
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/empty.smi";
    SmilesMolSupplier smiSup(infile, ",", 0, -1);
    TEST_ASSERT(smiSup.atEnd());
  }
  // tests for GitHub issue 19:
  {
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/empty2.sdf";
    SDMolSupplier reader(infile);
    TEST_ASSERT(reader.length()==0);
  }
  {
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/empty2.sdf";
    SDMolSupplier reader(infile);
    TEST_ASSERT(reader.atEnd());
    bool failed=false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed=true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length()==0);
  }
  {
    std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/empty2.sdf";
    SDMolSupplier reader(infile);
    bool failed=false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed=true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length()==0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    TEST_ASSERT(reader.atEnd());
    bool failed=false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed=true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length()==0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    bool failed=false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed=true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length()==0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    TEST_ASSERT(reader.length()==0);
  }
  
  
}

void testCisTrans() {
  std::string text;
  text = "mol-1 ClC(C)=C(Br)C\n"
    "mol-2 C1=COC=CC1C(Cl)=C(Br)C\n"
    "mol-3 C1=COC=CC1\\C(Cl)=C(Br)\\C";
  SmilesMolSupplier smiSup;
  smiSup.setData(text," ",1,0,false,true);
  
  std::string ofile = "cisTrans.sdf";
  SDWriter writer(ofile);
  while(!smiSup.atEnd()){
    ROMol *mol = smiSup.next();
    TEST_ASSERT(mol);
    RDDepict::compute2DCoords(*mol);
    writer.write(*mol);
    delete mol;
  }
  // do the round test
  // parse the sd file and write it out to smiles
  
  SDMolSupplier *reader;
  try{
    reader = new SDMolSupplier("cisTrans.sdf");
  } catch (FileParseException &) {
    reader=0;
  }
  TEST_ASSERT(reader);
  while (!reader->atEnd()) {
    ROMol *mol = reader->next();
    std::string mname;
    mol->getProp("_Name", mname);
    BOOST_LOG(rdInfoLog) << mname << " ";
    BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, 1) << "\n";
    delete mol;
  }
  delete reader;
}

void testStereoRound() {
  // - we will read ina bunch of cdk2 smiles with stereo on them
  // - generate the canonical smiles for each one
  // - generate 2D coordinates, write to an sdf file
  // - read the sdf file back in and compare the canonical smiles
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/cdk2_stereo.csv";
  SmilesMolSupplier *smiSup;
  try{
    smiSup = new SmilesMolSupplier(infile, ",", 0, 1, false, true);
  } catch (FileParseException &){
    smiSup=0;
  }
  TEST_ASSERT(smiSup)
  std::map<std::string, std::string> nameSmi;
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/cdk2_stereo.sdf";
  SDWriter* writer = new SDWriter(ofile);
  int count = 0;
  
  while(!smiSup->atEnd()) {
    ROMol *mol = smiSup->next();
    //mol->debugMol(std::cout);
    std::string mname;
    mol->getProp("_Name", mname);
    nameSmi[mname] = MolToSmiles(*mol, 1);
      
    ROMol *nmol = SmilesToMol(nameSmi[mname]);
    //nmol->debugMol(std::cout);
      
    std::string nsmi = MolToSmiles(*nmol, 1);
    //BOOST_LOG(rdErrorLog) << mname << "\n";
    if (nameSmi[mname] != nsmi) {
      BOOST_LOG(rdInfoLog) << mname << " " << nameSmi[mname] << " " << nsmi << "\n";
    }
    RDDepict::compute2DCoords(*mol);
    writer->write(*mol);
    count++;
    delete mol;
    delete nmol;
      
    if (count%50 == 0) {
      BOOST_LOG(rdInfoLog) << count << " " << mname << "\n";
    }
   }
  delete smiSup;
  delete writer;
  
  // now read the SD file back in check if the canonical smiles are the same
  SDMolSupplier *reader;
  try{
    reader = new SDMolSupplier(ofile);
  } catch (FileParseException &) {
    reader=0;
  }
  TEST_ASSERT(reader);
  count = 0;
  
  while (!reader->atEnd()) {
    ROMol *mol = reader->next();
    //mol->debugMol(std::cout);
    std::string smiles = MolToSmiles(*mol, 1);
    std::string mname;
    mol->getProp("_Name", mname);
    if (nameSmi[mname] != smiles) {
      BOOST_LOG(rdInfoLog) << mname << " " << nameSmi[mname] << " " << smiles << "\n";
    }
    count++;
  }
  delete reader;

}

void testIssue226() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue226.sdf";
  SDMolSupplier sdsup(fname);
  
  ROMol *mol;

  mol = sdsup.next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("E1"));
  TEST_ASSERT(mol->hasProp("E2"));
  delete mol;
  mol = sdsup.next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("E1"));
  TEST_ASSERT(mol->hasProp("E2"));
}


int testTDTSupplier1() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  {
    TDTMolSupplier suppl(fname,"PN");
    unsigned int i=0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1,prop2;
        TEST_ASSERT(nmol->getNumAtoms()>0);
        TEST_ASSERT(nmol->hasProp("PN"));
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("MFCD"));

        nmol->getProp("PN",prop1);
        nmol->getProp("_Name",prop2);
        TEST_ASSERT(prop1==prop2);
      
        // we didn't ask for 2D conformers, so there should be a property 2D:
        TEST_ASSERT(nmol->hasProp("2D"));
        // and no conformer:
        TEST_ASSERT(!nmol->getNumConformers());
      
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i==10);
  }
  {
    std::ifstream strm(fname.c_str());
    TDTMolSupplier suppl(&strm,false,"PN");
    unsigned int i=0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1,prop2;
        TEST_ASSERT(nmol->getNumAtoms()>0);
        TEST_ASSERT(nmol->hasProp("PN"));
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("MFCD"));

        nmol->getProp("PN",prop1);
        nmol->getProp("_Name",prop2);
        TEST_ASSERT(prop1==prop2);
      
        // we didn't ask for 2D conformers, so there should be a property 2D:
        TEST_ASSERT(nmol->hasProp("2D"));
        // and no conformer:
        TEST_ASSERT(!nmol->getNumConformers());
      
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i==10);
  }
  return 1;
}
int testTDTSupplier2() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  int i;
  std::string prop1,prop2;
  
  TDTMolSupplier suppl(fname,"PN",2);
  i=0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms()>0);
      TEST_ASSERT(nmol->hasProp("PN"));
      TEST_ASSERT(nmol->hasProp("_Name"));
      TEST_ASSERT(nmol->hasProp("MFCD"));

      nmol->getProp("PN",prop1);
      nmol->getProp("_Name",prop2);
      TEST_ASSERT(prop1==prop2);
      
      // we asked for 2D conformers, so there should be no property 2D:
      TEST_ASSERT(!nmol->hasProp("2D"));
      // and a conformer:
      TEST_ASSERT(nmol->getNumConformers()==1);
      // with id "2":
      TEST_ASSERT(nmol->beginConformers()->get()->getId()==2);
      
      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i==10);
  return 1;
}
int testTDTSupplier3() {
  std::string data;
  int i;
  std::string prop1,prop2;
  
  TDTMolSupplier suppl;

  data="$SMI<Cc1nnc(N)nc1C>\n"
    "CAS<17584-12-2>\n"
    "|\n"
    "$SMI<Cc1n[nH]c(=O)nc1N>\n"
    "CAS<~>\n"
    "|\n"
    "$SMI<Cc1n[nH]c(=O)[nH]c1=O>\n"
    "CAS<932-53-6>\n"
    "|\n"
    "$SMI<Cc1nnc(NN)nc1O>\n"
    "CAS<~>\n"
    "|\n";
  suppl.setData(data,"CAS");

  i=0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms()>0);
      TEST_ASSERT(nmol->hasProp("CAS"));
      TEST_ASSERT(nmol->hasProp("_Name"));

      nmol->getProp("CAS",prop1);
      nmol->getProp("_Name",prop2);
      TEST_ASSERT(prop1==prop2);
      
      // no conformers should have been read:
      TEST_ASSERT(nmol->getNumConformers()==0);

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i==4);
  TEST_ASSERT(suppl.length()==4);

  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  ROMol *mol = suppl[0];
  TEST_ASSERT(mol);
  delete mol;

  // make sure we can reset the supplier and still process it properly;
  suppl.setData(data,"CAS");

  i=0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms()>0);
      TEST_ASSERT(nmol->hasProp("CAS"));
      TEST_ASSERT(nmol->hasProp("_Name"));

      nmol->getProp("CAS",prop1);
      nmol->getProp("_Name",prop2);
      TEST_ASSERT(prop1==prop2);
      
      // no conformers should have been read:
      TEST_ASSERT(nmol->getNumConformers()==0);

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i==4);

  return 1;
}


void testSDSupplierFromText() {
  std::string text;
  int i = 0;
  SDMolSupplier reader;

  text="Structure1\n"
    "csChFnd70/05230312262D\n"
    "\n"
    "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
    "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  4  2  0  0  0  0\n"
    "  1  5  1  0  0  0  0\n"
    "M  END\n"
    ">  <ID> (3)\n"
    "Lig1\n"
    "\n"
    "$$$$\n"
    "Structure1\n"
    "csChFnd70/05230312262D\n"
    "\n"
    "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
    "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.8477    0.6988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  4  2  0  0  0  0\n"
    "  1  5  1  0  0  0  0\n"
    "  3  6  1  0  0  0  0\n"
    "M  END\n"
    ">  <ID> (4)\n"
    "Lig2\n"
    "\n"
    "$$$$\n";
  reader.setData(text);
  
  i=0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    TEST_ASSERT(mol->hasProp("_Name"));
    TEST_ASSERT(mol->hasProp("ID"));
    i++;
    delete mol;
  }
  TEST_ASSERT(i==2);
}

void testIssue265() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NotThere.sdf";
  bool ok=false;
  try {
    SDMolSupplier reader(fname);
    ok=false;
  } catch (BadFileException &) {
    ok=true;
  }
  TEST_ASSERT(ok);

  try {
    SmilesMolSupplier reader(fname);
    ok = false;
  } catch (BadFileException &) {
    ok=true;
  }
  TEST_ASSERT(ok);

  try {
    TDTMolSupplier reader(fname);
    ok = false;
  } catch (BadFileException &) {
    ok=true;
  }
  TEST_ASSERT(ok);
}

void testSDErrorHandling() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors1.sdf";
  SDMolSupplier *sdsup;
  ROMol *nmol=0;  

  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(!nmol->hasProp("ID"));
  delete sdsup;
   
  // case 2: can't be sanitized
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors2.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;

  // entry 3: bad number of atoms
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors3.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;
  
  // entry 4: bad number of bonds
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors4.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;
}

void testIssue381() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue381.sdf";
  SDMolSupplier *sdsup;

  ROMol *nmol=0;  
  int count;
  
  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  count = 0;
  while(!sdsup->atEnd()){
    nmol = sdsup->next();
    if(nmol) delete nmol;
    count++;
  }
  TEST_ASSERT(sdsup->atEnd());
  TEST_ASSERT(count==9);
  
  TEST_ASSERT(sdsup->length()==9);

  delete sdsup;
}

void testSetStreamIndices() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  int tmpIndices[16]={0, 2136, 6198, 8520, 11070, 12292, 14025, 15313, 17313, 20125, 22231,
				     23729, 26147, 28331, 32541, 33991};
  std::vector<std::streampos> indices(tmpIndices,tmpIndices+16);
  SDMolSupplier *sdsup;

  ROMol *nmol=0;  
  int count;
  
  sdsup = new SDMolSupplier(fname);
  sdsup->setStreamIndices(indices);
  TEST_ASSERT(!sdsup->atEnd());
  TEST_ASSERT(sdsup->length()==16);
  
  count = 0;
  while(!sdsup->atEnd()){
    nmol = sdsup->next();
    if(nmol) delete nmol;
    count++;
  }
  TEST_ASSERT(sdsup->atEnd());
  TEST_ASSERT(count==16);
  
  TEST_ASSERT(sdsup->length()==16);

  delete sdsup;
}

int testMixIterAndRandom() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/esters.sdf";
  bool ok;

  SDMolSupplier *sdsup;
  ROMol *mol;
  std::string name;
  
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(sdsup);
  unsigned int i = 0;
  while (!sdsup->atEnd()) {
    mol = sdsup->next();
    if (mol) {
      TEST_ASSERT(mol->hasProp("ID"));
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i==6);
  TEST_ASSERT(sdsup->length()==6);

  delete sdsup;
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(sdsup);
  TEST_ASSERT(sdsup->length()==6);

  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID",name);
  TEST_ASSERT(name=="Lig1");    
  delete mol;
  
  mol = (*sdsup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID",name);
  TEST_ASSERT(name=="Lig1");    
  delete mol;
  
  sdsup->reset();
  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID",name);
  TEST_ASSERT(name=="Lig1");    
  delete mol;
  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID",name);
  TEST_ASSERT(name=="Lig2");    
  delete mol;
  delete sdsup;
  
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  TEST_ASSERT(nSup->length()==10);
  mol=(*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol=(*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(nSup->length()==10);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol=nSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(nSup->length()==10);
  delete mol;
  mol=(*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(nSup->length()==10);
  delete mol;
  mol=nSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==20);
  delete nSup;
  delete mol;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol=0;
  try {
    mol=(*nSup)[20];
    ok = false;
  } catch (FileParseException &) {
    ok=true;
  }
  TEST_ASSERT(ok);
  delete nSup;


  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  TDTMolSupplier *tSup;
  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  TEST_ASSERT(tSup->length()==10);
  mol=(*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol=(*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(tSup->length()==10);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol=tSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(tSup->length()==10);
  delete mol;
  mol=(*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==9);
  TEST_ASSERT(tSup->length()==10);
  delete mol;
  mol=tSup->next();
  TEST_ASSERT(mol);
  mol=tSup->next();
  TEST_ASSERT(mol);
  mol=tSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms()==10);
  delete tSup;
  delete mol;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol=0;
  try {
    mol=(*tSup)[20];
    ok = false;
  } catch (FileParseException &) {
    ok=true;
  }
  TEST_ASSERT(ok);
  delete tSup;


    
  return 1;
}


int testRemoveHs() {
  
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/withHs.sdf";
  
  SDMolSupplier sdsup(fname);
  ROMol *nmol;

  nmol = sdsup.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==23);
  delete nmol;
  nmol = sdsup.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==28);
  delete nmol;

  std::cerr<<"build:"<<std::endl;
  SDMolSupplier sdsup2(fname,true,false);
  nmol = sdsup2.next();
  TEST_ASSERT(nmol);
  //std::cerr<<" count: "<<nmol->getNumAtoms()<<std::endl;
  TEST_ASSERT(nmol->getNumAtoms()==39);
  delete nmol;
  nmol = sdsup2.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==30);
  delete nmol;


  return 1;
}

void testGetItemText() {
  std::string rdbase = getenv("RDBASE");
  std::string fname;

  ROMol *mol1,*mol2;
  std::string molB,smiles;
  bool ok;

  {
    fname=rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);
    TEST_ASSERT(sdsup.length()==16);
  
    molB = sdsup.getItemText(0);
    mol1 = sdsup[0];
    TEST_ASSERT(mol1);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms()==mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    // make sure getItemText() doesn't screw up the current position:
    molB = sdsup.getItemText(10);
    mol1 = sdsup.next();
    molB = sdsup.getItemText(1);
    TEST_ASSERT(mol1);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms()==mol1->getNumAtoms());
    delete mol1;
    delete mol2;
  
    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = sdsup.getItemText(15);
    mol1 = sdsup[15];
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms()==mol1->getNumAtoms());
    delete mol1;
    delete mol2;
  
    try {
      molB=sdsup.getItemText(16);
      ok = false;
    } catch (FileParseException &) {
      ok=true;
    }
    TEST_ASSERT(ok);

    try {
      molB=sdsup.getItemText(20);
      ok = false;
    } catch (FileParseException &) {
      ok=true;
    }
    TEST_ASSERT(ok);
  }

  {
    fname=rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);

    // make sure getItemText() works if we haven't read at all from the supplier:
    // (this was sf.net issue 2632960)
    molB = sdsup.getItemText(0);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    mol1 = sdsup[0];
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol2->getNumAtoms()==mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    molB = sdsup.getItemText(5);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms()==16);
    mol1 = sdsup[5];
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol2->getNumAtoms()==mol1->getNumAtoms());
    delete mol1;
    delete mol2;

  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname,",",1,0,false);
    TEST_ASSERT(smisup.length()==10);
  
    molB = smisup.getItemText(0);
    TEST_ASSERT(molB=="1, CC1=CC(=O)C=CC1=O, 34.14");
    mol1 = smisup[0];
    TEST_ASSERT(mol1);
    delete mol1;

    molB = smisup.getItemText(5);
    TEST_ASSERT(molB=="6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");
    mol1 = smisup.next();
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol1->getNumAtoms()==20);
    delete mol1;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = smisup.getItemText(8);
    TEST_ASSERT(molB=="9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB=="10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");

    mol1 = smisup[0];
    TEST_ASSERT(mol1);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles=="CC1=CC(=O)C=CC1=O");
    TEST_ASSERT(mol1->getNumAtoms()==9);
    delete mol1;

  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname,",",1,0,false);

    // make sure getItemText() works if we haven't read at all from the supplier:
    // (this was sf.net issue 2632960)
    molB = smisup.getItemText(0);
    TEST_ASSERT(molB=="1, CC1=CC(=O)C=CC1=O, 34.14");

    molB = smisup.getItemText(5);
    TEST_ASSERT(molB=="6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");

    molB = smisup.getItemText(8);
    TEST_ASSERT(molB=="9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB=="10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname,",",1,0,false);

    // make sure getItemText() flags EOF
    // (this was sf.net issue 3299878)
    molB = smisup.getItemText(0);
    TEST_ASSERT(molB=="1, CC1=CC(=O)C=CC1=O, 34.14");

    ROMol *m=smisup[9];
    TEST_ASSERT(m);
    delete m;
    TEST_ASSERT(smisup.atEnd());
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB=="10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
    TEST_ASSERT(smisup.atEnd());
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    // make sure getItemText() works if we haven't read at all from the supplier:
    // (this was sf.net issue 2632960)
    molB = tdtsup.getItemText(0);
    TEST_ASSERT(molB!="");
  
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    TEST_ASSERT(tdtsup.length()==10);
  
    molB = tdtsup.getItemText(0);
    TEST_ASSERT(molB!="");

    mol1 = tdtsup[0];
    TEST_ASSERT(mol1);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles=="Cc1nnc(N)nc1C");
    TEST_ASSERT(mol1->getNumAtoms()==9);
    delete mol1;

    // make sure getItemText doesn't screw up next()
    molB = tdtsup.getItemText(5);
    mol1 = tdtsup.next();
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol1->getNumAtoms()==9);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles=="Cc1n[nH]c(=O)nc1N");
    delete mol1;


    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = tdtsup.getItemText(9);
    TEST_ASSERT(molB!="");
    TEST_ASSERT(molB.substr(0,12)=="$SMI<Cc1n[nH");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    TEST_ASSERT(tdtsup.length()==10);
  
    ROMol *mol=tdtsup[9];
    TEST_ASSERT(mol);
    delete mol;
    TEST_ASSERT(tdtsup.atEnd());

    // (this was sf.net issue 3299878
    molB = tdtsup.getItemText(9);
    TEST_ASSERT(molB!="");
    TEST_ASSERT(molB.substr(0,12)=="$SMI<Cc1n[nH");
    TEST_ASSERT(tdtsup.atEnd());
  }

}


int testForwardSDSupplier() {
  
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  std::string fname2 = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf.gz";
  
  {
    std::ifstream strm(fname.c_str());
    ForwardSDMolSupplier sdsup(&strm,false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      TEST_ASSERT(nmol || sdsup.atEnd());
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i==16);
  }
#ifdef TEST_GZIP_SD
  // make sure the boost::iostreams are working
  {
    io::filtering_istream strm;
    strm.push(io::file_source(fname));

    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm,line);
      ++i;
      if(i>1000) break;
    }
    TEST_ASSERT(i==998);
  }
  {
    io::filtering_istream strm;
    strm.push(io::gzip_decompressor());
    strm.push(io::file_source(fname2));

    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm,line);
      ++i;
      if(i>1000) break;
    }
    TEST_ASSERT(i==998);
  }
  // looks good, now do a supplier:
  {
    io::filtering_istream strm;
    strm.push(io::gzip_decompressor());
    strm.push(io::file_source(fname2));

    ForwardSDMolSupplier sdsup(&strm,false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp("_Name"));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i==16);
  }
#endif
  return 1;
}

void testMissingCRSDSupplier() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/missingCR.sdf";
  SDMolSupplier reader(infile);
  reader.next();
  TEST_ASSERT(reader.atEnd());
}


void testIssue3482695() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3482695.sdf";
  SDMolSupplier reader(infile);
  ROMol *nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==0);
  TEST_ASSERT(nmol->hasProp("test"));
  delete nmol;
}

void testIssue3525673() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525673.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins,false);
  ROMol *nmol;

  nmol = reader.next();
  TEST_ASSERT(!nmol);

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==37);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(!nmol);

  nmol = reader.next();
  TEST_ASSERT(!nmol);

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==58);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(!nmol);

  nmol = reader.next();
  TEST_ASSERT(!nmol);

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==58);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(!nmol);

}

void testBlankLinesInProps() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/BlankPropLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins,false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==19);
  TEST_ASSERT(nmol->hasProp("MultiLineProperty1"));
  nmol->getProp("MultiLineProperty1",pval);
  TEST_ASSERT(pval=="foo\nbar\n \nbaz");
  TEST_ASSERT(nmol->hasProp("MultiLineProperty2"));
  TEST_ASSERT(!(nmol->hasProp("fooprop")));
  nmol->getProp("MultiLineProperty2",pval);
  TEST_ASSERT(pval=="foo\n>  <fooprop>\nbaz\n ");
  delete nmol;

}

void testSkipLines() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/SkipLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins,false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==1);
  TEST_ASSERT(nmol->hasProp("prop1"));
  delete nmol;
}


void testGitHub23() {
  std::string rdbase = getenv("RDBASE");
  std::string ofile = rdbase + "/Code/GraphMol/FileParsers/test_data/blah.sdf";
  SDWriter *writer = new SDWriter(ofile);

  ROMol *mol = SmilesToMol("CCCC");
  INT_VECT iv;
  iv.push_back(1);
  iv.push_back(2);  
  mol->setProp("pval",iv);
  writer->write(*mol);
  delete mol;

  writer->flush();
  delete writer;
}

void testGitHub88() {
  std::string rdbase = getenv("RDBASE");
  std::string infile = rdbase + "/Code/GraphMol/FileParsers/test_data/github88.v3k.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins,false);
  ROMol *nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms()==8);
  TEST_ASSERT(nmol->hasProp("prop1"));
  std::string pval;
  nmol->getProp("prop1",pval);
  TEST_ASSERT(pval=="4");
  delete nmol;
}

int main() {
  RDLog::InitLogs();

#if 1
  BOOST_LOG(rdErrorLog) <<"\n-----------------------------------------\n";
  testMolSup();
  BOOST_LOG(rdErrorLog) <<"Finished: testMolSup()\n";
  BOOST_LOG(rdErrorLog) <<"-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testRandMolSup();
  BOOST_LOG(rdErrorLog) <<"Finished: testRandMolSup()\n";
  BOOST_LOG(rdErrorLog) <<"-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesSup();
  BOOST_LOG(rdErrorLog) <<"Finished: testSmilesSup()\n";
  BOOST_LOG(rdErrorLog) <<"-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesSupFromText();
  BOOST_LOG(rdErrorLog) <<"Finished: testSmilesSupFromText()\n";
  BOOST_LOG(rdErrorLog) <<"-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesWriter();
  BOOST_LOG(rdErrorLog) <<"Finished: testSmilesWriter()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDWriter();
  BOOST_LOG(rdErrorLog) <<"Finished: testSDWriter()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierEnding();
  BOOST_LOG(rdErrorLog) <<"Finished: testSDSupplierEnding()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSuppliersEmptyFile();
  BOOST_LOG(rdErrorLog) <<"Finished: testSuppliersEmptyFile()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
  
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testCisTrans();
  BOOST_LOG(rdErrorLog) <<"Finished: testCisTrans()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testStereoRound();
  BOOST_LOG(rdErrorLog) <<"Finished: testStereoRound()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue226();
  BOOST_LOG(rdErrorLog) <<"Finished: testIssue226()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier1();
  BOOST_LOG(rdErrorLog) <<"Finished: testTDTSupplier1()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier2();
  BOOST_LOG(rdErrorLog) <<"Finished: testTDTSupplier2()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier3();
  BOOST_LOG(rdErrorLog) <<"Finished: testTDTSupplier3()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierFromText();
  BOOST_LOG(rdErrorLog) <<"Finished: testSDSupplierFromText()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";


  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue265();
  BOOST_LOG(rdErrorLog) <<"Finished: testIssue265()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDErrorHandling();
  BOOST_LOG(rdErrorLog) <<"Finished: testSDErrorHandling()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue381();
  BOOST_LOG(rdErrorLog) <<"Finished: testIssue381()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSetStreamIndices();
  BOOST_LOG(rdErrorLog) <<"Finished: testSetStreamIndices()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testMixIterAndRandom();
  BOOST_LOG(rdErrorLog) <<"Finished: testMixIterAndRandom()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testRemoveHs();
  BOOST_LOG(rdErrorLog) <<"Finished: testRemoveHs()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGetItemText();
  BOOST_LOG(rdErrorLog) <<"Finished: testGetItemText()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testForwardSDSupplier();
  BOOST_LOG(rdErrorLog) <<"Finished: testForwardSDSupplier()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testMissingCRSDSupplier();
  BOOST_LOG(rdErrorLog) <<"Finished: testMissingCRSDSupplier()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue3482695();
  BOOST_LOG(rdErrorLog) <<"Finished: testIssue3482695()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
#endif
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue3525673();
  BOOST_LOG(rdErrorLog) <<"Finished: testIssue3525673()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testBlankLinesInProps();
  BOOST_LOG(rdErrorLog) <<"Finished: testBlankLinesInProps()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSkipLines();
  BOOST_LOG(rdErrorLog) <<"Finished: testSkipLines()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub23();
  BOOST_LOG(rdErrorLog) <<"Finished: testGitHub23()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub88();
  BOOST_LOG(rdErrorLog) <<"Finished: testGitHub88()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  return 0;
}
