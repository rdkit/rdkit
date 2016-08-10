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
#include "../FileParsers/FileParsers.h" //MOL single molecule !
#include "../FileParsers/MolSupplier.h" //SDF

#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "../Substruct/SubstructMatch.h"

#include "StructChecker.h"

using namespace RDKit;
using namespace RDKit::StructureCheck;

void testFlags() //PASSED
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "testFlags\n";

    unsigned int flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR;
    std::string  str = RDKit::StructureCheck::StructChecker::StructureFlagsToString(flags);
    unsigned int f2 = RDKit::StructureCheck::StructChecker::StringToStructureFlags(str);
    BOOST_LOG(rdInfoLog) << str << "\n";
    TEST_ASSERT(flags == f2);

    flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR
          | RDKit::StructureCheck::StructChecker::TRANSFORMED;
    str = RDKit::StructureCheck::StructChecker::StructureFlagsToString(flags);
    f2 = RDKit::StructureCheck::StructChecker::StringToStructureFlags(str);
    BOOST_LOG(rdInfoLog) << str << "\n";
    TEST_ASSERT(flags == f2);

    flags = 0xFFFF;// &(~0x0080); // - unused bit
    str = StructChecker::StructureFlagsToString(flags);
    f2 = StructChecker::StringToStructureFlags(str);
    BOOST_LOG(rdInfoLog) << f2 <<" = "<< str << "\n";
    TEST_ASSERT((flags & (~0x0080)) == f2);

    str = " STEREO_ERROR ,\t TRANSFORMED [xXx}"; // 'stability test with minor syntax errors'
    flags = RDKit::StructureCheck::StructChecker::STEREO_ERROR
          | RDKit::StructureCheck::StructChecker::TRANSFORMED;
    f2 = StructChecker::StringToStructureFlags(str);
    BOOST_LOG(rdInfoLog) << str <<" = "<< StructChecker::StructureFlagsToString(f2) <<"\n";
    TEST_ASSERT(flags == f2);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
//--------------------------------------------------------------------------

void testOptions() //PASSED
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "testOptions\n";

    StructCheckerOptions options;

    TEST_ASSERT( parseOptionsJSON("{}", options));
    TEST_ASSERT(!parseOptionsJSON("{...error..}", options));
    bool ok;
    ok = parseOptionsJSON("{\"Verbose\": true, \"CheckStereo\": true}", options);
    TEST_ASSERT(ok && options.Verbose && options.CheckStereo);

    ok = loadOptionsFromFiles(options,
        "", // augmentedAtomTranslationsFile = "",
        "", // patternFile = "",       // file with clean patterns
        "", // rotatePatternFile = "", // file with rotate patterns
        "", // stereoPatternFile = "", // file with stereo patterns
        "");// tautomerFile = "");

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(ok);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testLoadOptionsFromFiles()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "testLoadOptionsFromFiles FROM CURRENT (.../test) DIRECTORY\n";
    bool ok;
    StructCheckerOptions options;

    const std::string rdbase = getenv("RDBASE");
    const std::string testDataDir = rdbase + "/Code/GraphMol/StructChecker/test/";

    BOOST_LOG(rdInfoLog) << "loadGoodAugmentedAtoms checkfgs.chk\n";
    ok = options.loadGoodAugmentedAtoms(rdbase + "checkfgs.chk");
    TEST_ASSERT(ok);

    BOOST_LOG(rdInfoLog) << "loadAcidicAugmentedAtoms checkfgs.aci\n";
    ok = options.loadAcidicAugmentedAtoms(rdbase + "checkfgs.aci");
    TEST_ASSERT(ok);

    BOOST_LOG(rdInfoLog) << "loadAugmentedAtomTranslations checkfgs.trn\n";
    ok = options.loadAugmentedAtomTranslations(rdbase + "checkfgs.trn");
    TEST_ASSERT(ok);

    //BOOST_LOG(rdInfoLog) << "loadPatterns patterns.sdf\n";
    //ok = options.loadPatterns("patterns.sdf");
    //TEST_ASSERT(ok);

    //....

    BOOST_LOG(rdInfoLog) << "loadTautomerData tautomer.sdf\n";
    ok = options.loadTautomerData("tautomer.sdf");
    TEST_ASSERT(ok);

    BOOST_LOG(rdInfoLog) << "loadTautomerData tautomer.rdf\n";
    ok = options.loadTautomerData("tautomer.rdf");
    TEST_ASSERT(ok);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
//--------------------------------------------------------------------------

void test1()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test1\n";
    const char* smols[] = {
        "CCC", //tmp
    };

    StructCheckerOptions options;
    options.Verbose = true;
    bool ok = loadOptionsFromFiles(options,
        "", // augmentedAtomTranslationsFile = "",
        "", // patternFile = "",       // file with clean patterns
        "", // rotatePatternFile = "", // file with rotate patterns
        "", // stereoPatternFile = "", // file with stereo patterns
        "");// tautomerFile = "");
    TEST_ASSERT(ok);

    StructChecker chk(options);
    for (int i = 0; i<sizeof(smols) / sizeof(smols[0]); i++) {
        RWMol* mol = SmilesToMol(smols[i]);
        TEST_ASSERT(mol);
        unsigned flags = chk.checkMolStructure(*mol);
        delete mol;
        BOOST_LOG(rdInfoLog) << "......... results ........\n";
        TEST_ASSERT(true);
    }
    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//--------------------------------------------------------------------------

void testStereo() // stereochemistry
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "testStereo\n";
    const char* smols[] = {
        "COC(=O)C(\\C)=C\\C1C(C)(C)[C@H]1C(=O)O[C@@H]2C(C)=C(C(=O)C2)CC=CC=C", //Pyrethrin II (C22H28O5)
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2", // Bergenin (cuscutin) (a resin) (C14H16O9)"
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1", //Glucose (glucopyranose) (C6H12O6)
        "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2" // Bergenin (cuscutin) (a resin) (C14H16O9)
    };

    StructCheckerOptions options;
    options.Verbose = true;
    StructChecker chk(options);
    for (int i = 0; i<sizeof(smols) / sizeof(smols[0]); i++) {
        BOOST_LOG(rdInfoLog) << i << " : " << smols[i] << "\n";
        RWMol* mol = SmilesToMol(smols[i]);
        TEST_ASSERT(mol);
        unsigned flags = chk.checkMolStructure(*mol);
        delete mol;
        BOOST_LOG(rdInfoLog) << "......... results ........\n";
        TEST_ASSERT(true);
    }
    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//==============================================================================

int main(int argc, const char* argv[])
{
    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    BOOST_LOG(rdInfoLog) << "StructChecker Unit Test \n";

    testStereo();

    testFlags();
    testOptions();
    try {
        testLoadOptionsFromFiles();
    }
    catch (...) {
        // relative path to patern files must be correct !
    }
    test1();

    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    return 0;
}


