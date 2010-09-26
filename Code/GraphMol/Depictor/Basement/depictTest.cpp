// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "Depictor.h"
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>


using namespace RDKit;

#ifdef WIN32_DLLBUILD
void test1(){
  BOOST_LOG(rdInfoLog)<< "------- Test1: round-trip via Mol file" << std::endl;
  std::string smi="CC(=O)O";
  std::string fName="out.mol";
  SmilesToMolFileDLL(smi,fName);
  ROMol *mol1 = (ROMol *)SmilesToMol(smi);
  ROMol *mol2 = (ROMol *)MolFileToMol(fName);

  std::string refSmi=MolToSmiles(*mol1);
  std::string testSmi=MolToSmiles(*mol2);
  TEST_ASSERT(refSmi==testSmi);
  delete mol1;
  delete mol2;

  BOOST_LOG(rdInfoLog)<< "------- Done" << std::endl;
}
#endif

void test2(){
  BOOST_LOG(rdInfoLog)<< "------- Test2: adding 2D coords" << std::endl;
  std::string smi="CC(=O)O";
  ROMol *mol = (ROMol *)SmilesToMol(smi);
  std::string refSmi=MolToSmiles(*mol);
  bool success;
  success=Add2DCoordsToMol(*mol);
  TEST_ASSERT(success);

  // test handling of chirality:
  delete mol;
  smi="F[C@H](Cl)Br";
  mol = (ROMol *)SmilesToMol(smi);
  success=Add2DCoordsToMol(*mol);
  TEST_ASSERT(success);
  delete mol;
  BOOST_LOG(rdInfoLog)<< "------- Done" << std::endl;
}


void main(){
  RDLog::InitLogs();
#ifdef WIN32_DLLBUILD
  test1();
#endif
  test2();
}
