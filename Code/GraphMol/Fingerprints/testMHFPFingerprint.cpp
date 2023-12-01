//
//  2019, Daniel Probst, Reymond Group @ University of Bern
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

#include <GraphMol/Fingerprints/MHFP.h>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <GraphMol/Fingerprints/MHFP.h>

using namespace RDKit;

void testMHFPInit() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test MHFP fingerprint encoder initialization"
                        << std::endl;

  std::string s = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
  std::string t = "Cn1cnc2c1c(=O)[nH]c(=O)n2C";

  ROMol *mol_s = SmilesToMol(s);
  ROMol *mol_t = SmilesToMol(t);

  MHFPFingerprints::MHFPEncoder enc(128, 42);
  auto fp_s = enc.Encode(s);
  auto fp_t = enc.Encode(t);

  auto fp_mol_s = enc.Encode(*mol_s);
  auto fp_mol_t = enc.Encode(*mol_t);

  TEST_ASSERT(fp_s.size() == 128);
  TEST_ASSERT(fp_s[0] == fp_mol_s[0]);
  TEST_ASSERT(fp_t[127] == fp_mol_t[127]);

  delete mol_s;
  delete mol_t;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMHFPHashing() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test MHFP hashing of string and uint arrays"
                        << std::endl;

  MHFPFingerprints::MHFPEncoder enc(8);

  std::vector<uint32_t> input_a = {1, 2, 4, 5, 6, 7, 8, 9};

  std::vector<uint32_t> output_a = {188049437,  364485576, 737251017,
                                    810894466,  300249621, 154369992,
                                    2221926165, 283729444};

  TEST_ASSERT(enc.FromArray(input_a) == output_a);

  std::vector<std::string> input_b = {"a", "b", "c", "d", "e", "f"};

  std::vector<uint32_t> output_b = {631555539, 835857365, 445245415, 4162827301,
                                    955545975, 943207071, 712975995, 363547692};

  TEST_ASSERT(enc.FromStringArray(input_b) == output_b);

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMHFPShingling() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test MHFP shingling creation" << std::endl;

  std::string s = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
  std::string t = "Cn1cnc2c1c(=O)[nH]c(=O)n2C";

  ROMol *mol_s = SmilesToMol(s);

  MHFPFingerprints::MHFPEncoder enc;

  auto fp_s = enc.Encode(s);
  auto fp_t = enc.Encode(t);

  auto sh_a = enc.CreateShingling(s);
  auto sh_b = enc.CreateShingling(*mol_s);

  TEST_ASSERT(sh_a.size() == 44);
  TEST_ASSERT(sh_b.size() == 44);

  TEST_ASSERT(enc.CreateShingling(s, 3, false).size() == 42);
  TEST_ASSERT(enc.CreateShingling(s, 3, true, false, true, 0).size() == 58);

  delete mol_s;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMHFPSECFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test SECFP fingerprint functionality" << std::endl;

  std::string s = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";

  ROMol *mol_s = SmilesToMol(s);

  MHFPFingerprints::MHFPEncoder enc;
  auto fp_s = enc.EncodeSECFP(s, 3, true, false, true, 1, 16);

  auto fp_mol_s = enc.EncodeSECFP(*mol_s, 3, true, false, true, 1, 16);

  TEST_ASSERT(fp_s.size() == 16);
  TEST_ASSERT(fp_s[10]);
  TEST_ASSERT(fp_s[15]);

  delete mol_s;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMHFPDistance() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test MHFP shingling creation" << std::endl;

  std::string s = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
  std::string t = "Cn1cnc2c1c(=O)[nH]c(=O)n2C";

  MHFPFingerprints::MHFPEncoder enc;

  auto fp_s = enc.Encode(s);
  auto fp_t = enc.Encode(t);

  TEST_ASSERT(feq(MHFPFingerprints::MHFPEncoder::Distance(fp_s, fp_s), 0.0));

  TEST_ASSERT(
      feq(MHFPFingerprints::MHFPEncoder::Distance(fp_s, fp_t), 0.7109375));

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  testMHFPInit();
  testMHFPShingling();
  testMHFPHashing();
  testMHFPSECFP();

  return 0;
}
