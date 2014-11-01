// $Id$
//
//  Copyright (C) 2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>

#include <iostream>
using namespace std;
using namespace RDKit;


void testBasics()
{
  BOOST_LOG(rdInfoLog)  << "-----------------------\n Basic Allocations" << std::endl;
  Atom *a1 = new Atom(6);
  Bond *b1 = new Bond();
  ROMol *m1 = new ROMol();

  a1 = NULL; // intentional leak
  BOOST_LOG(rdInfoLog)  << "Finished" << std::endl;
}

void testSMILES()
{
  BOOST_LOG(rdInfoLog)  << "-----------------------\n SMILES Read" << std::endl;
  string smi="CCOC";
  ROMol *m = SmilesToMol(smi);
  smi="C1COC1";
  ROMol *m2 = SmilesToMol(smi);

  BOOST_LOG(rdInfoLog)  << "Finished" << std::endl;
}


// -------------------------------------------------------------------
int main()
{
  RDLog::InitLogs();
  //boost::logging::enable_logs("rdApp.info");
  // test1();  // <- this doesn't seem to actually do anything
#if 1
  testBasics();
  testSMILES();
#endif

  return 0;
}
