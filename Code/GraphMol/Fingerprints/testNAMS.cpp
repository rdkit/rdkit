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
  "CN(CCc1c[nH]c2c1cccc2)C", // The canonicalized order from NAMS
  "NCCc1c[nH]c2c1cccc2" // The canonicalized order from NAMS
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
    double ssim = calcSelfSimilarity(molinfos[ii], params);
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

      std::unique_ptr< NAMSResult > result( getNAMSResult(molinfos[ii], molinfos[jj], params) );
      //BOOST_LOG(rdErrorLog) << "For " << molinfos[ii].smiles << " against " << molinfos[jj].smiles << " found similarity " << result->similarity << " expected " << CROSSSIM[ii][offset] << '\n';
      TEST_ASSERT ( result->similarity < CROSSSIM[ii][offset] + 0.001 && result->similarity > CROSSSIM[ii][offset] - 0.001 );
      ++count;
    }
  }

  auto stop = std::chrono::high_resolution_clock::now();
  BOOST_LOG(rdErrorLog) << "Calculated similarity of " << count << " pairings in " << std::chrono::duration<double>(stop - start).count() << " seconds\n";
  BOOST_LOG(rdErrorLog) << "    Including generating " << molinfos.size() << " data objects it took " << std::chrono::duration<double>(stop - fullstart).count() << " seconds\n";

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

template< typename T >
std::string print_vector( const std::vector< T > & vec ) {
  std::stringstream ss;
  ss << "[ ";
  for ( auto const & entry: vec ) {
    ss << entry << ", ";
  }
  ss << ']';
  return ss.str();
}

void compare_vectors( const std::vector< int > & test, const std::vector< int > & ref, const std::string & tag ) {
  if ( test.size() != ref.size() ) {
    BOOST_LOG(rdErrorLog) << "For " << tag << " vector size of " << test.size() << " doesn't match reference size of " << ref.size() << "\n";
    BOOST_LOG(rdErrorLog) << '\t' << print_vector( test ) << '\n';
    BOOST_LOG(rdErrorLog) << '\t' << print_vector( ref ) << '\n';
  }
  TEST_ASSERT( test.size() == ref.size() );
  for ( unsigned int ii=0; ii < test.size(); ++ii ) {
    if ( test[ii] != ref[ii] ) {
      BOOST_LOG(rdErrorLog) << "For " << tag << " vector element " << ii << " does not match" << "\n";
      BOOST_LOG(rdErrorLog) << '\t' << print_vector( test ) << '\n';
      BOOST_LOG(rdErrorLog) << '\t' << print_vector( ref ) << '\n';
      TEST_ASSERT( test[ii] == ref[ii] );
    }
  }
}

void compare_vectors( const std::vector< double > & test, const std::vector< double > & ref, const std::string & tag, double delta = 0.01 ) {
  if ( test.size() != ref.size() ) {
    BOOST_LOG(rdErrorLog) << "For " << tag << " vector size of " << test.size() << " doesn't match reference size of " << ref.size() << "\n";
    BOOST_LOG(rdErrorLog) << '\t' << print_vector( test ) << '\n';
    BOOST_LOG(rdErrorLog) << '\t' << print_vector( ref ) << '\n';
  }
  TEST_ASSERT( test.size() == ref.size() );
  for ( unsigned int ii=0; ii < test.size(); ++ii ) {
    if ( test[ii] < ref[ii] - delta || test[ii] > ref[ii] + delta ) {
      BOOST_LOG(rdErrorLog) << "For " << tag << " vector element " << ii << " does not match" << "\n";
      BOOST_LOG(rdErrorLog) << '\t' << print_vector( test ) << '\n';
      BOOST_LOG(rdErrorLog) << '\t' << print_vector( ref ) << '\n';
      TEST_ASSERT( test[ii] > ref[ii] - delta && test[ii] < ref[ii] + delta );
    }
  }
}

void testNAMSPairing() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test NAMS pairing calculation"
                        << std::endl;

  ROMOL_SPTR mol1( SmilesToMol(SMILES[0]) );
  ROMOL_SPTR mol11( SmilesToMol(SMILES[10]) );
  ROMOL_SPTR mol12( SmilesToMol(SMILES[11]) );

  NAMSMolInfo mi1( *mol1 );
  NAMSMolInfo mi11( *mol11 );
  NAMSMolInfo mi12( *mol12 );

  NAMSParameters params;

  std::unique_ptr< NAMSResult > result1_11( getNAMSResult(mi1, mi11, params) );
  std::unique_ptr< NAMSResult > result11_1( getNAMSResult(mi11, mi1, params) );
  std::unique_ptr< NAMSResult > result11_12( getNAMSResult(mi11, mi12, params) );

  std::vector< int > ref_map1_11 = {1, 2, 3, 5, 4, 10, 9, 12, 7, 8};
  std::vector< int > ref_map11_1 = {-1, 0, 1, 2, 4, 3, -1, 8, 5, 6, 9, -1, 7, -1};
  std::vector< int > ref_map11_12 = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1};
  std::vector< double > ref_scores1_11 = {4.24, 4.71, 5.23, 4.89, 5.79, 5.29, 5.68, 5.26, 5.38, 5.86};
  std::vector< double > ref_scores11_1 = {0, 4.24, 4.71, 5.23, 5.79, 4.89, 0, 5.38, 5.86, 5.68, 5.29, 0, 5.26, 0};
  std::vector< double > ref_scores11_12 = {0, 8.89, 9.76, 10.51, 11.18, 10.66, 10.66, 11.11, 11.35, 10.74, 10.21, 10.02, 10.53, 0 };

  compare_vectors( ref_map1_11, result1_11->mapping1to2, "Cmpd 1 to 11 map" );
  compare_vectors( ref_map11_1, result11_1->mapping1to2, "Cmpd 11 to 1 map" );
  compare_vectors( ref_map11_12, result11_12->mapping1to2, "Cmpd 11 to 12 map" );

  compare_vectors( ref_scores1_11, result1_11->atom_scores, "Cmpd 1 to 11 scores", 0.01 );
  compare_vectors( ref_scores11_1, result11_1->atom_scores, "Cmpd 11 to 1 scores", 0.01 );
  compare_vectors( ref_scores11_12, result11_12->atom_scores, "Cmpd 11 to 12 scores", 0.01 );

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  testNAMSMolInfo();
  testNAMSSelfSimilarity();
  testNAMSCrossSimiliarity();
  testNAMSPairing();

  return 0;
}
