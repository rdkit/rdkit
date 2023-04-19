//
//  Copyright (C) 2023 Rocco Moretti and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <GraphMol/Fingerprints/NAMS.h>

#include <vector>
#include <string>

using namespace RDKit;
using namespace RDKit::NAMS;

std::vector< std::string > SMILES = {
  "C=CCC1CCCCC1",
  "C1CCCCC1",
  "CC=C",
  "OC=C",
  "CC",
  "CCCCC",
  "c1ccccc1",
  "Cc1ccccc1",
  "CCc1ccccc1",
  "Cc1ccccc1O",
  "CN(C)CCc1c[nH]c2ccccc12",
  "NCCc1c[nH]c2ccccc12"
};

std::vector< float > SSIM = {
  81.009,
  32.279,
  5.780,
  5.780,
  2.000,
  17.979,
  32.279,
  43.199,
  55.009,
  161.679,
  125.859
};

void testNAMSMolInfo() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test NAMS MolInfo loading"
                        << std::endl;

  std::vector< ROMOL_SPTR > mols;
  for ( std::string const & smi: SMILES ) {
    mols.push_back( ROMOL_SPTR(SmilesToMol(smi)) );
  }

  std::vector< NAMSMolInfo > molinfos;
  for ( ROMOL_SPTR const & mol: mols ) {
    molinfos.emplace_back( *mol ); // Rely on constructor that takes a ROMol object
  }

  // We're mainly testing that things don't crash/segfault,
  // but we can also check that we've got all the atoms.

  for ( unsigned int ii=0; ii< mols.size(); ++ii ) {
    TEST_ASSERT( molinfos[ii].natoms() == mols[ii]->getNumAtoms() );
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testNAMSSelfSimilarity() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test NAMS self similarity"
                        << std::endl;

  std::vector< ROMOL_SPTR > mols;
  for ( std::string const & smi: SMILES ) {
    mols.push_back( ROMOL_SPTR(SmilesToMol(smi)) );
  }

  std::vector< NAMSMolInfo > molinfos;
  for ( ROMOL_SPTR const & mol: mols ) {
    molinfos.emplace_back( *mol ); // Rely on constructor that takes a ROMol object
  }

  NAMSParameters params;

  for ( unsigned int ii=0; ii< mols.size(); ++ii ) {
    float ssim = calcSelfSimilarity(molinfos[ii], params)/10000.0f;
    TEST_ASSERT( ssim == SSIM[ii] );
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  testNAMSMolInfo();
  testNAMSSelfSimilarity();

  return 0;
}
