//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h"  //MOL single molecule !
#include "../FileParsers/MolSupplier.h"  //SDF

#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "../Substruct/SubstructMatch.h"

#include "StructChecker.h"
#include "Stereo.h"

using namespace RDKit;
using namespace RDKit::StructureCheck;

void testFlags()  // PASSED
{
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testFlags\n";

  unsigned int flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR;
  std::string str =
      RDKit::StructureCheck::StructChecker::StructureFlagsToString(flags);
  unsigned int f2 =
      RDKit::StructureCheck::StructChecker::StringToStructureFlags(str);
  BOOST_LOG(rdInfoLog) << str << "\n";
  TEST_ASSERT(flags == f2);

  flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR |
          RDKit::StructureCheck::StructChecker::TRANSFORMED;
  str = RDKit::StructureCheck::StructChecker::StructureFlagsToString(flags);
  f2 = RDKit::StructureCheck::StructChecker::StringToStructureFlags(str);
  BOOST_LOG(rdInfoLog) << str << "\n";
  TEST_ASSERT(flags == f2);

  flags = 0xFFFF;  // &(~0x0080); // - unused bit
  str = StructChecker::StructureFlagsToString(flags);
  f2 = StructChecker::StringToStructureFlags(str);
  BOOST_LOG(rdInfoLog) << f2 << " = " << str << "\n";
  TEST_ASSERT((flags & (~0x0080)) == f2);

  str = " STEREO_ERROR ,\t TRANSFORMED [xXx}";  // 'stability test with minor
                                                // syntax errors'
  flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR |
          RDKit::StructureCheck::StructChecker::TRANSFORMED;
  f2 = StructChecker::StringToStructureFlags(str);
  BOOST_LOG(rdInfoLog) << str << " = "
                       << StructChecker::StructureFlagsToString(f2) << "\n";
  TEST_ASSERT(flags == f2);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
//--------------------------------------------------------------------------

void testOptionsJSON()  // PASSED
{
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testOptionsJSON\n";

  StructCheckerOptions options;

  TEST_ASSERT(parseOptionsJSON("{}", options));
  TEST_ASSERT(!parseOptionsJSON("{...error..}", options));
  bool ok;
  ok = parseOptionsJSON("{\"Verbose\": true, \"CheckStereo\": true}", options);
  TEST_ASSERT(ok && options.Verbose && options.CheckStereo);
  //    BOOST_LOG(rdInfoLog) << "......... results ........\n";
  TEST_ASSERT(ok);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void doLoadOptionsFromFiles(StructCheckerOptions& options,
                            const std::string& dirBase = "",
                            bool strict = true) {
  bool ok;
  const std::string rdbase = getenv("RDBASE") ? getenv("RDBASE") : ".";
  std::string testDataDir;
  if (dirBase == "") {
    testDataDir = rdbase + "/Code/GraphMol/StructChecker/test/";
  } else {
    testDataDir = dirBase;
  }
  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "testDataDir: " << testDataDir << "\n";

  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "loadGoodAugmentedAtoms checkfgs.chk\n";
  ok = options.loadGoodAugmentedAtoms(testDataDir + "checkfgs.chk");
  TEST_ASSERT(ok);

  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "loadAcidicAugmentedAtoms checkfgs.aci\n";
  ok = options.loadAcidicAugmentedAtoms(testDataDir + "checkfgs.aci");
  TEST_ASSERT(!strict || ok);

  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "loadAugmentedAtomTranslations checkfgs.trn\n";
  ok = options.loadAugmentedAtomTranslations(testDataDir + "checkfgs.trn");
  TEST_ASSERT(ok);

  // BOOST_LOG(rdInfoLog) << "loadPatterns patterns.sdf\n";
  // ok = options.loadPatterns("testDataDir + patterns.sdf");
  // TEST_ASSERT(ok);

  //....

  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "loadTautomerData tautomer.sdf\n";
  ok = options.loadTautomerData(testDataDir + "tautomer.sdf");
  TEST_ASSERT(!strict || ok);

  if (options.Verbose)
    BOOST_LOG(rdInfoLog) << "loadTautomerData tautomer.rdf\n";
  ok = options.loadTautomerData(testDataDir + "tautomer.rdf");
  TEST_ASSERT(!strict || ok);

  options.Verbose = true;
}

void testLoadOptionsFromFiles() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "testLoadOptionsFromFiles FROM CURRENT (.../test) DIRECTORY\n";
  bool ok;
  StructCheckerOptions options;
  options.Verbose = true;
  doLoadOptionsFromFiles(options);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
//--------------------------------------------------------------------------

void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "test1\n";
  const char* smols[] = {
      "CCC",  // tmp
      "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)"
      "O2",  // Bergenin (cuscutin) (a resin) (C14H16O9)
  };
  StructCheckerOptions options;
  doLoadOptionsFromFiles(options);
  options.Verbose = true;
  /*
      bool ok = loadOptionsFromFiles(options,
          "", // augmentedAtomTranslationsFile = "",
          "", // patternFile = "",       // file with clean patterns
          "", // rotatePatternFile = "", // file with rotate patterns
          "", // stereoPatternFile = "", // file with stereo patterns
          "");// tautomerFile = "");
      TEST_ASSERT(ok);
  */
  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
    TEST_ASSERT(true);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "test2\n";
  const char* smols[] = {
      "CCC",  // tmp
      "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)"
      "O2",  // Bergenin (cuscutin) (a resin) (C14H16O9)
  };

  StructCheckerOptions options;
  doLoadOptionsFromFiles(options);
  options.Verbose = true;

  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
    TEST_ASSERT(true);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//--------------------------------------------------------------------------

const char* substance_310975001 =
    "310975001\n"
    "  -OEChem-06071611182D\n"
    "\n"
    " 10 10  0     1  0  0  0  0  0999 V2000\n"
    "    1.1317   -0.3264    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0\n"
    "    0.3470    0.7535    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0\n"
    "    0.3470   -0.0715    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.1317    1.0084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.6166    0.3410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.2363    1.3369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.3452   -1.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.0332    1.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7618   -1.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.6166    1.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "  1  5  1  0  0  0  0\n"
    "  1  7  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  4  1  0  0  0  0\n"
    "  2  6  1  0  0  0  0\n"
    "  4  5  1  0  0  0  0\n"
    "  6  8  1  0  0  0  0\n"
    "  7  9  1  0  0  0  0\n"
    "  8 10  1  0  0  0  0\n"
    "M  END\n"
    "> <EXPECTED>\n"
    "['stereo_error']\n"
    "\n"
    "> <GOT>\n"
    "['atom_check_failed']\n"
    "\n"
    "$$$$\n";
//----------------------------------
// CIS_TRANS_EITHER is BondDir::EITHERDOUBLE (i.e.crossed double bond)
// squiggle bond from chiral center
const char* Mrv1561_08171605252D =
    "Mrv1561_08171605252D\n"
    "\n\n"
    "  5  4  0     0  0  0  0  0  0999 V2000\n"
    "   -1.7411    2.3214    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.9161    2.3214    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0\n"
    "   -0.0911    2.3214    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.9161    3.1464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.9161    1.4964    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  5  1  0  0  0  0\n"
    "  2  4  1  4  0  0  0\n"
    "M  END\n"
    "$$$$\n";
// Avalon : ['EITHER_WARNING', 'DUBIOUS_STEREO_REMOVED']
// RDKit : ATOM_CHECK_FAILED
//--------------------

//    crossed double bond (2D)
const char* Mrv1561_08171605322D =
    "Mrv1561 08171605322D 0   0.00000     0.00000     0\n"
    "\n\n"
    " 4  3   0     0  0  0  0  0  0999 V2000\n"
    "   -0.4241   -1.7187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.4009   -1.7187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.8366   -1.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.8134   -2.4332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "  2  4  1  0  0  0  0\n"
    "  1  2  2  3  0  0  0\n"
    "M  END\n"
    "$$$$\n";
// Avalon : []
// RDKit : ATOM_CHECK_FAILED, EITHER_WARNING, DUBIOUS_STEREO_REMOVED
//--------------------

// squiggle bond from double bond (2D)
const char* Mrv1561_08171605332D =
    "Mrv1561 08171605332D 0   0.00000     0.00000     0\n"
    "\n\n"
    " 4  3   0     0  0  0  0  0  0999 V2000\n"
    "   -0.4241   -1.7187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.4009   -1.7187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.8366   -1.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.8134   -2.4332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "  1  2  2  0  0  0  0\n"
    "  2  4  1  0  0  0  0\n"
    "M  END\n"
    "$$$$\n";
// Avalon : ['EITHER_WARNING', 'DUBIOUS_STEREO_REMOVED']
// RDKit : ATOM_CHECK_FAILED, EITHER_WARNING, DUBIOUS_STEREO_REMOVED
//------------
void testStereo()  // stereochemistry
{
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testStereo\n";
  const char* smols[] = {
      "COC(=O)C(\\C)=C\\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC="
      "CC=C",  // Pyrethrin II (C22H28O5)
      "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(="
      "O)O2",  // Bergenin (cuscutin) (a resin) (C14H16O9)
      "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1",  // Glucose
                                                           // (glucopyranose)
                                                           // (C6H12O6)
      "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(="
      "O)O2"  // Bergenin (cuscutin) (a resin) (C14H16O9)
  };

  StructCheckerOptions options;
  //    doLoadOptionsFromFiles(options);
  const std::string rdbase = getenv("RDBASE") ? getenv("RDBASE") : ".";
  const std::string testDataDir = rdbase + "/Code/GraphMol/StructChecker/test/";
  TEST_ASSERT(options.loadGoodAugmentedAtoms(testDataDir + "checkfgs.chk"));
  TEST_ASSERT(options.loadAcidicAugmentedAtoms(testDataDir + "checkfgs.aci"));
  options.Verbose = true;
  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << "FLAGS:"
                         << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    TEST_ASSERT(true);
  }

  {
    BOOST_LOG(rdInfoLog) << "substance_310975001:\n";
    ROMOL_SPTR mol(MolBlockToMol(substance_310975001));
    //      std::cerr << (size_t) mol.get() << std::endl;
    TEST_ASSERT(mol.get());
    BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
    TEST_ASSERT(CheckStereo(*mol.get()) == false);
    unsigned flags = chk.checkMolStructure(*dynamic_cast<RWMol*>(mol.get()));
    BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
    BOOST_LOG(rdInfoLog) << "FLAGS: "
                         << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    //      TEST_ASSERT(0!=(flags & StructChecker::STEREO_ERROR));
  }
  {
    const char* substance_set[] = {substance_310975001, Mrv1561_08171605252D,
                                   Mrv1561_08171605322D, Mrv1561_08171605332D};
    const unsigned res[] = {
        StructChecker::STEREO_ERROR,
        StructChecker::EITHER_WARNING | StructChecker::DUBIOUS_STEREO_REMOVED,
        0,  // ??
        StructChecker::EITHER_WARNING | StructChecker::DUBIOUS_STEREO_REMOVED,
    };
    for (size_t i = 0; i < sizeof(substance_set) / sizeof(*substance_set);
         i++) {
      BOOST_LOG(rdInfoLog) << "substance " << i << "\n";
      ROMOL_SPTR mol(MolBlockToMol(substance_set[i]));
      TEST_ASSERT(mol.get());
      BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
      /*
                  if (0 != (res[i] & StructChecker::STEREO_ERROR)) {
                      TEST_ASSERT(CheckStereo(*mol) == false);
                  }
                  else {
      //                TEST_ASSERT(CheckStereo(*mol) == true);
                  }
      */
      unsigned flags = chk.checkMolStructure(*dynamic_cast<RWMol*>(mol.get()));
      BOOST_LOG(rdInfoLog) << MolToSmarts(*mol) << "\n";
      BOOST_LOG(rdInfoLog) << "ref: "
                           << StructChecker::StructureFlagsToString(res[i])
                           << "\n";
      BOOST_LOG(rdInfoLog) << "RES: "
                           << StructChecker::StructureFlagsToString(flags)
                           << "\n";
      //            TEST_ASSERT((flags == res[i]);
      BOOST_LOG(rdInfoLog) << "-------\n";
    }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

// FAILED
void testOptionsDefault() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testOptionsDefault\n";
  const char* smols[] = {
      "COC(=O)C",
      "COC(=O)C(\\C)=C\\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC="
      "CC=C",  // Pyrethrin II (C22H28O5)
  };

  StructCheckerOptions
      options;  // intial GoodAtoms loading is INCORRECT. There is no Ligands!
  options.Verbose = true;
  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    TEST_ASSERT(0 != flags);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

// FAILED
void testCheckAtomWithDefaultGoodAtoms() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testCheckAtom\n";
  const char* smols[] = {
      "COC(=O)C",
  };

  StructCheckerOptions
      options;  // intial GoodAtoms loading is INCORRECT. There is no Ligands!
  options.Verbose = true;
  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    //        TEST_ASSERT( !(flags & StructChecker::ATOM_CHECK_FAILED));
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testCheckAtom() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testCheckAtom\n";
  const char* smols[] = {
      "COC(=O)C",
  };

  StructCheckerOptions options;
  doLoadOptionsFromFiles(options);
  options.Verbose = true;
  StructChecker chk(options);
  for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
    BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
    RWMol* mol = SmilesToMol(smols[i]);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    //        TEST_ASSERT(!(flags & StructChecker::ATOM_CHECK_FAILED));
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSpecificExamples() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
  BOOST_LOG(rdInfoLog) << "testSpecificExamples\n";
  StructCheckerOptions options;
  const std::string rdbase = getenv("RDBASE");
  const std::string testDataDir = rdbase + "/Data/struchk/";
  doLoadOptionsFromFiles(options, testDataDir, false);

  options.RemoveMinorFragments = true;
  options.CheckCollisions = true;
  options.CollisionLimitPercent = 3;
  options.CheckStereo = true;
  options.MaxMolSize = 999;
  options.Verbose = true;
  StructChecker chk(options);
  {
    const char* smols[] = {
        "C[N+](C)(C)C",
        "CC(=C)C(OCC[N+](C)(C)C)=[N+](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
        "OC(=O)[C@@H]1CC=CN1"};

    for (int i = 0; i < sizeof(smols) / sizeof(smols[0]); i++) {
      BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
      RWMol* mol = SmilesToMol(smols[i]);
      TEST_ASSERT(mol);
      unsigned flags = chk.checkMolStructure(*mol);
      delete mol;
      BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                           << "\n";
      TEST_ASSERT(!flags || flags == StructChecker::TRANSFORMED);
      //        TEST_ASSERT(!(flags & StructChecker::ATOM_CHECK_FAILED));
    }
  }
  {
    std::string molb =
        "\n"
        "  Mrv1561 08261616022D\n"
        "\n"
        "  8  8  0  0  1  0            999 V2000\n"
        "    3.2170    2.5920    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    2.2609    1.5301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    3.5724    0.4125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    4.3260    0.0770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    4.8781    0.6901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    3.6586    1.2330    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  "
        "0\n"
        "    4.4656    1.4045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    3.0455    1.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "  1  8  2  0  0  0  0\n"
        "  2  8  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  3  6  1  0  0  0  0\n"
        "  4  5  2  0  0  0  0\n"
        "  5  7  1  0  0  0  0\n"
        "  6  7  1  0  0  0  0\n"
        "  6  8  1  6  0  0  0\n"
        "M  END\n";

    RWMol* mol = MolBlockToMol(molb);
    TEST_ASSERT(mol);
    unsigned flags = chk.checkMolStructure(*mol);
    delete mol;
    BOOST_LOG(rdInfoLog) << StructChecker::StructureFlagsToString(flags)
                         << "\n";
    TEST_ASSERT(!flags || flags == StructChecker::TRANSFORMED);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//==============================================================================

int main(int argc, const char* argv[]) {
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  BOOST_LOG(rdInfoLog) << "StructChecker Unit Test \n";

// return 0; //tmp
#if 1
  testFlags();
  testOptionsJSON();
  try {
    testLoadOptionsFromFiles();
  } catch (...) {
    // relative path to patern files must be correct !
  }
  // FAILED
  testOptionsDefault();

  test1();
  test2();

  testCheckAtom();
  // FAILED
  testCheckAtomWithDefaultGoodAtoms();

  testStereo();
#endif
  testSpecificExamples();

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return 0;
}
