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
#include <chrono>

using namespace RDKit;
using namespace RDKit::NAMS;

std::vector< std::string > SMILES = {
  "OC=CCC1CCCCC1",
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

// Reference values computed by the official NAMS code (https://github.com/aofalcao/nams-docker/)

const std::vector< float > SSIM = {
  81.009,
  32.279,
  5.780,
  5.780,
  2.000,
  17.979,
  32.279,
  43.199,
  55.009,
  55.719,
  161.679,
  125.859
};

const std::vector< std::vector< float > > CROSSSIM = {
  { 32.039,  5.540,  5.540,  1.820, 14.426, 25.583, 35.569, 46.561, 42.365, 52.321, 52.495 },
  {  3.001,  1.498,  1.040, 10.068, 25.823, 25.763, 25.703, 25.703, 22.711, 22.831 },
  {  3.179,  1.980,  5.431,  2.539,  3.491,  4.791,  3.461,  4.581,  4.641 },
  {  1.780,  2.830,  1.325,  2.046,  2.652,  3.946,  2.442,  2.502 },
  {  1.940,  0.820,  1.480,  1.860,  1.460,  1.720,  1.760 },
  {  8.106,  9.650, 12.050,  9.600, 11.413, 11.513 },
  { 32.219, 32.159, 32.159, 28.471, 28.591 },
  { 43.129, 43.129, 38.113, 38.253 },
  { 48.153, 48.981, 49.141 },
  { 43.806, 43.966 },
  { 125.619 }
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
    //BOOST_LOG(rdErrorLog) << "For " << molinfos[ii].smiles << " found self similarity " << ssim << " expected " << SSIM[ii] << '\n';
    TEST_ASSERT ( ssim < SSIM[ii] + 0.001 && ssim > SSIM[ii] - 0.001 );
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testNAMSCrossSimiliarity() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test NAMS cross similarity"
                        << std::endl;

  std::vector< ROMOL_SPTR > mols;
  for ( std::string const & smi: SMILES ) {
    mols.push_back( ROMOL_SPTR(SmilesToMol(smi)) );
  }

  auto fullstart = std::chrono::high_resolution_clock::now();

  std::vector< NAMSMolInfo > molinfos;
  for ( ROMOL_SPTR const & mol: mols ) {
    molinfos.emplace_back( *mol ); // Rely on constructor that takes a ROMol object
  }

  NAMSParameters params;

  auto start = std::chrono::high_resolution_clock::now();

  int count = 0;
  for ( unsigned int ii=0; ii< mols.size(); ++ii ) {
    for ( unsigned int jj=ii+1; jj <mols.size(); ++jj ) {
      unsigned int offset( jj-ii-1 );

      float sim = nams_runner(molinfos[ii], molinfos[jj], params)/10000.0f;
      //BOOST_LOG(rdErrorLog) << "For " << molinfos[ii].smiles << " against " << molinfos[jj].smiles << " found similarity " << sim << " expected " << CROSSSIM[ii][offset] << '\n';
      TEST_ASSERT ( sim < CROSSSIM[ii][offset] + 0.001 && sim > CROSSSIM[ii][offset] - 0.001 );
      ++count;
    }
  }

  auto stop = std::chrono::high_resolution_clock::now();
  BOOST_LOG(rdErrorLog) << "Calculated similarity of " << count << " pairings in " << std::chrono::duration<double>(stop - start).count() << " seconds\n";
  BOOST_LOG(rdErrorLog) << "    Including generating " << molinfos.size() << " data objects it took " << std::chrono::duration<double>(stop - fullstart).count() << " seconds\n";

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  testNAMSMolInfo();
  testNAMSSelfSimilarity();
  testNAMSCrossSimiliarity();

  return 0;
}
