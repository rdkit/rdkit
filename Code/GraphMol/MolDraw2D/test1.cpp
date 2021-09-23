//
//  Copyright (C) 2015-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <algorithm>
#include <cstdio>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <map>
#include <fstream>
#include <sstream>

#ifdef RDKIT_USE_BOOST_REGEX
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_search;
#else
#include <regex>
using std::regex;
using std::regex_search;
#endif

#ifdef RDK_BUILD_FREETYPE_SUPPORT
static const std::map<std::string, std::hash_result_t> SVG_HASHES = {
    {"test1_1.svg", 660888398U},
    {"test1_2.svg", 3756518769U},
    {"test1_3.svg", 2466979491U},
    {"test1_4.svg", 3312001249U},
    {"test1_5.svg", 3680433096U},
    {"test1_6.svg", 2302426373U},
    {"test4_1.svg", 3135047047U},
    {"test5_1.svg", 4098074962U},
    {"test5_2.svg", 786626445U},
    {"test5_3.svg", 4145248145U},
    {"test6_1.svg", 2605479453U},
    {"test6_2.svg", 23881357U},
    {"test7_1.svg", 602351070U},
    {"test7_2.svg", 1022417818U},
    {"testGithub781_1.svg", 1687494516U},
    {"testGithub781_2.svg", 3948011009U},
    {"testGithub781_3.svg", 4068544025U},
    {"testGithub781_4.svg", 558275653U},
    {"testGithub781_5.svg", 2840696942U},
    {"test3_1.svg", 3931228883U},
    {"test3_2.svg", 1838614840U},
    {"test3_3.svg", 3265666846U},
    {"test3_4.svg", 2797781880U},
    {"test3_5.svg", 2305295191U},
    {"test3_6.svg", 2987593096U},
    {"test3_7.svg", 1173444046U},
    {"test774_1.svg", 2168226972U},
    {"test774_2.svg", 291877070U},
    {"test9_1.svg", 4293107790U},
    {"test852_1.svg", 2160661054U},
    {"test852_2.svg", 4028481516U},
    {"test852_2a.svg", 2656413440U},
    {"test852_2b.svg", 3528483680U},
    {"test852_2c.svg", 3770633776U},
    {"test852_2d.svg", 403158481U},
    {"test860_1.svg", 3379787385U},
    {"test860_2.svg", 1696713930U},
    {"test860_3.svg", 100405646U},
    {"test910_1.svg", 3222118951U},
    {"test910_2.svg", 106389875U},
    {"test983_1.svg", 3375032549U},
    {"test983_2.svg", 4030659164U},
    {"testNoDeuterium.svg", 2058910502U},
    {"testNoTritium.svg", 2838762713U},
    {"testDeuterium.svg", 1753646811U},
    {"testTritium.svg", 3479638945U},
    {"crossed_bonds.svg", 2936035516U},
    {"test10_1.svg", 966111909U},
    {"test10_2.svg", 949107593U},
    {"test10_3.svg", 4122471991U},
    {"test10_4.svg", 1647099728U},
    {"test10_5.svg", 308826031U},
    {"test10_6.svg", 2601983922U},
    {"test11_1.svg", 2982846572U},
    {"test11_2.svg", 717239576U},
    {"test12_1.svg", 4001334642U},
    {"test12_3.svg", 2718439435U},
    {"test12_4.svg", 2718439435U},
    {"test12_2.svg", 4250647988U},
    {"test13_1.svg", 3533854997U},
    {"testGithub1090_1.svg", 2535099142U},
    {"test1271_1.svg", 3250186735U},
    {"test1271_2.svg", 429631964U},
    {"test1271_3.svg", 3324480069U},
    {"test1271_4.svg", 3324480069U},
    {"test1271_5.svg", 2755333160U},
    {"test1322_1.svg", 1595409146U},
    {"test1322_2.svg", 3317974613U},
    {"test14_1.svg", 1611937867U},
    {"test14_2.svg", 2405055622U},
    {"test15_1.svg", 2891524973U},
    {"test15_2.svg", 238458373U},
    {"test17_1.svg", 867199271U},
    {"test17_2.svg", 3551160267U},
    {"test17_3.svg", 3931471413U},
    {"test17_4.svg", 2364427473U},
    {"test18_1.svg", 2622450140U},
    {"test18_2.svg", 2988606034U},
    {"test18_3.svg", 393339448U},
    {"test18_4.svg", 1084560479U},
    {"test18_5.svg", 393339448U},
    {"test19_1.svg", 2622450140U},
    {"test19_2.svg", 226882495U},
    {"test16_1.svg", 247758111U},
    {"test16_2.svg", 252943288U},
    {"testGithub2063_1.svg", 2158365418U},
    {"testGithub2063_2.svg", 2158365418U},
    {"testGithub2151_1.svg", 739229913U},
    {"testGithub2151_2.svg", 3095783108U},
    {"testGithub2762.svg", 320023621U},
    {"testGithub2931_1.svg", 2708252551U},
    {"testGithub2931_2.svg", 2520170665U},
    {"test20_1.svg", 4167740575U},
    {"test20_2.svg", 3940601415U},
    {"test20_3.svg", 256044325U},
    {"test20_4.svg", 3636815083U},
    {"test21_1.svg", 3896401240U},
    {"test21_2.svg", 2346906091U},
    {"test22_1.svg", 3352529761U},
    {"test22_2.svg", 3623180832U},
    {"test23_1.svg", 1474697688U},
    {"testGithub3112_1.svg", 4164182512U},
    {"testGithub3112_2.svg", 1727510065U},
    {"testGithub3112_3.svg", 1355264573U},
    {"testGithub3112_4.svg", 1958409504U},
    {"testGithub3305_1.svg", 3352529761U},
    {"testGithub3305_2.svg", 1309416130U},
    {"testGithub3305_3.svg", 2451552479U},
    {"testGithub3305_4.svg", 2305295191U},
    {"testGithub3305_5.svg", 3635251234U},
    {"testGithub3305_6.svg", 1996529499U},
    {"testGithub3305_7.svg", 2003696554U},
    {"testGithub3391_1.svg", 3216272110U},
    {"testGithub3391_2.svg", 4241749225U},
    {"testGithub3391_3.svg", 2634130947U},
    {"testGithub3391_4.svg", 3347048145U},
    {"testGithub4156_1.svg", 3587507438U},
    {"testGithub4156_2.svg", 617944934U},
    {"testGithub4496_1.svg", 3884377499U}};
#else
static const std::map<std::string, std::hash_result_t> SVG_HASHES = {
    {"test1_1.svg", 2810365302U},
    {"test1_2.svg", 436586098U},
    {"test1_3.svg", 2057814158U},
    {"test1_4.svg", 3558596694U},
    {"test1_5.svg", 3888582395U},
    {"test1_6.svg", 3118807497U},
    {"test4_1.svg", 2036014054U},
    {"test5_1.svg", 1774913296U},
    {"test5_2.svg", 786626445U},
    {"test5_3.svg", 3167056946U},
    {"test6_1.svg", 2863113368U},
    {"test6_2.svg", 23881357U},
    {"test7_1.svg", 602351070U},
    {"test7_2.svg", 1022417818U},
    {"testGithub781_1.svg", 1868801797U},
    {"testGithub781_2.svg", 645137970U},
    {"testGithub781_3.svg", 3557563258U},
    {"testGithub781_4.svg", 3922698686U},
    {"testGithub781_5.svg", 2840696942U},
    {"test3_1.svg", 1887265010U},
    {"test3_2.svg", 199107565U},
    {"test3_3.svg", 823104548U},
    {"test3_4.svg", 3019005132U},
    {"test3_5.svg", 2891359380U},
    {"test3_6.svg", 3999883443U},
    {"test3_7.svg", 3410424210U},
    {"test774_1.svg", 4249918177U},
    {"test774_2.svg", 3726741547U},
    {"test9_1.svg", 3881616369U},
    {"test852_1.svg", 1747936364U},
    {"test852_2.svg", 1995763141U},
    {"test852_2a.svg", 781368964U},
    {"test852_2b.svg", 855240733U},
    {"test852_2c.svg", 617421834U},
    {"test852_2d.svg", 239450370U},
    {"test860_1.svg", 3684477751U},
    {"test860_2.svg", 2255345324U},
    {"test860_3.svg", 2317297396U},
    {"test910_1.svg", 1390009092U},
    {"test910_2.svg", 358157210U},
    {"test983_1.svg", 3613290859U},
    {"test983_2.svg", 1239543852U},
    {"testNoDeuterium.svg", 1974538289U},
    {"testNoTritium.svg", 3403650090U},
    {"testDeuterium.svg", 1861778686U},
    {"testTritium.svg", 1327537979U},
    {"crossed_bonds.svg", 2936035516U},
    {"test10_1.svg", 1915067848U},
    {"test10_2.svg", 1607134571U},
    {"test10_3.svg", 1339712793U},
    {"test10_4.svg", 230053470U},
    {"test10_5.svg", 3338491447U},
    {"test10_6.svg", 3464662207U},
    {"test11_1.svg", 2664209367U},
    {"test11_2.svg", 2916457382U},
    {"test12_1.svg", 2960364537U},
    {"test12_3.svg", 753038675U},
    {"test12_4.svg", 753038675U},
    {"test12_2.svg", 936927131U},
    {"test13_1.svg", 1540336771U},
    {"testGithub1090_1.svg", 1292284867U},
    {"test1271_1.svg", 3250186735U},
    {"test1271_2.svg", 429631964U},
    {"test1271_3.svg", 3886623470U},
    {"test1271_4.svg", 3886623470U},
    {"test1271_5.svg", 2238686066U},
    {"test1322_1.svg", 2001812365U},
    {"test1322_2.svg", 1025658798U},
    {"test14_1.svg", 972851850U},
    {"test14_2.svg", 3049734321U},
    {"test15_1.svg", 1665312890U},
    {"test15_2.svg", 2368159640U},
    {"test17_1.svg", 3799185367U},
    {"test17_2.svg", 554175906U},
    {"test17_3.svg", 4132524395U},
    {"test17_4.svg", 2032383981U},
    {"test18_1.svg", 1937998685U},
    {"test18_2.svg", 3006977225U},
    {"test18_3.svg", 1655129816U},
    {"test18_4.svg", 1526792634U},
    {"test18_5.svg", 1655129816U},
    {"test19_1.svg", 1937998685U},
    {"test19_2.svg", 2158406784U},
    {"test16_1.svg", 3385311765U},
    {"test16_2.svg", 2735861810U},
    {"testGithub2063_1.svg", 4196564529U},
    {"testGithub2063_2.svg", 4196564529U},
    {"testGithub2151_1.svg", 1816894413U},
    {"testGithub2151_2.svg", 2185139506U},
    {"testGithub2762.svg", 1991842363U},
    {"testGithub2931_1.svg", 1459334551U},
    {"testGithub2931_2.svg", 2613847891U},
    {"test20_1.svg", 962316680U},
    {"test20_2.svg", 2230712067U},
    {"test20_3.svg", 195533391U},
    {"test20_4.svg", 1948583366U},
    {"test22_1.svg", 3352529761U},
    {"test22_2.svg", 2183584384U},
    {"test23_1.svg", 1655940990U},
    {"testGithub3112_1.svg", 2694789798U},
    {"testGithub3112_2.svg", 1912964589U},
    {"testGithub3112_3.svg", 1042493059U},
    {"testGithub3112_4.svg", 463691362U},
    {"testGithub3305_1.svg", 3352529761U},
    {"testGithub3305_2.svg", 1309416130U},
    {"testGithub3305_3.svg", 2451552479U},
    {"testGithub3305_4.svg", 2891359380U},
    {"testGithub3305_5.svg", 2648802668U},
    {"testGithub3305_6.svg", 1394712294U},
    {"testGithub3305_7.svg", 643160167U},
    {"testGithub3391_1.svg", 3043856820U},
    {"testGithub3391_2.svg", 899401090U},
    {"testGithub3391_3.svg", 1321562570U},
    {"testGithub3391_4.svg", 4052907129U},
    {"testGithub4496_1.svg", 2353591199U}};
#endif

// These PNG hashes aren't completely reliable due to floating point cruft,
// but they can still reduce the number of drawings that need visual
// inspection.  At present, the files
// test20_2.png  test3_5.png  test3_7.png  test774_2.png testGithub3305_5.png
// testGithub3305_7.png test2_2.png   test3_6.png  test5_1.png
// testGithub3305_4.png  testGithub3305_6.png
// give different results on my MBP and Ubuntu 20.04 VM.  The SVGs work
// better because the floats are all output to only 1 decimal place so there
// is a much smaller chance of different systems producing different files.
static const std::map<std::string, std::hash_result_t> PNG_HASHES = {
    {"test2_1.png", 3297082290U},
    {"test2_2.png", 103359060U},
    {"test2_3.png", 907340919U},
    {"test4_1.png", 2108306671U},
    {"test5_1.png", 2577909359U},
    {"test5_2.png", 3406225054U},
    {"test5_3.png", 2184033022U},
    {"test7_1.png", 2553239923U},
    {"test7_2.png", 4103514158U},
    {"test3_1.png", 3083663242U},
    {"test3_2.png", 1562500334U},
    {"test3_3.png", 69009007U},
    {"test3_4.png", 1894736182U},
    {"test3_5.png", 2575151141U},
    {"test3_6.png", 1485731281U},
    {"test3_7.png", 3257189464U},
    {"test774_1.png", 1790608443U},
    {"test774_2.png", 2217836790U},
    {"test852_1.png", 1671563761U},
    {"test852_2.png", 3595190908U},
    {"test860_1.png", 3676395309U},
    {"test860_2.png", 1193304483U},
    {"test860_3.png", 3150966270U},
    {"test20_1.png", 1403848138U},
    {"test20_2.png", 2342474031U},
    {"test20_3.png", 2287710945U},
    {"test20_4.png", 3550202242U},
    {"testGithub3305_1.png", 2793539378U},
    {"testGithub3305_2.png", 4013123469U},
    {"testGithub3305_3.png", 2854242566U},
    {"testGithub3305_4.png", 2575151141U},
    {"testGithub3305_5.png", 4034543887U},
    {"testGithub3305_6.png", 4150166195U},
    {"testGithub3305_7.png", 2478532666U}};

using namespace RDKit;

// if the generated SVG hashes to the value we're expecting, delete
// the file.  That way, only the files that need inspection will be
// left at the end of the run.
static const bool DELETE_WITH_GOOD_HASH = true;

std::hash_result_t hash_file(const std::string &filename) {
  std::ifstream ifs(filename, std::ios_base::binary);
  std::string file_contents(std::istreambuf_iterator<char>{ifs}, {});
  if (filename.substr(filename.length() - 4) == ".svg") {
    // deal with MSDOS newlines.
    file_contents.erase(
        remove(file_contents.begin(), file_contents.end(), '\r'),
        file_contents.end());
  }
  return gboost::hash_range(file_contents.begin(), file_contents.end());
}

void check_file_hash(const std::string &filename,
                     std::hash_result_t exp_hash = 0U) {
  std::cout << filename << " : " << hash_file(filename) << "U" << std::endl;

  std::map<std::string, std::hash_result_t>::const_iterator it;
  if (filename.substr(filename.length() - 4) == ".svg") {
    it = SVG_HASHES.find(filename);
  } else {
    it = PNG_HASHES.find(filename);
  }
  std::hash_result_t file_hash = hash_file(filename);
  if (exp_hash == 0U) {
    exp_hash = it == SVG_HASHES.end() ? 0U : it->second;
  }
  if (it != SVG_HASHES.end() && file_hash == exp_hash) {
    if (DELETE_WITH_GOOD_HASH) {
      std::remove(filename.c_str());
    }
  } else {
    std::cout << "file " << filename << " gave hash " << file_hash
              << "U not the expected " << exp_hash << "U" << std::endl;
  }
}

void test1() {
  std::cout << " ----------------- Test 1" << std::endl;
  {
    std::string smiles = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test1_1.svg");
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test1_1.svg");
    delete m;
  }
  {
    // make sure this works with the stringstream too:
    std::string smiles = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    TEST_ASSERT(text.find("<svg") != std::string::npos);
    TEST_ASSERT(text.find("</svg>") != std::string::npos);
    delete m;
  }
  {
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test1_2.svg");
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test1_2.svg");
    delete m;
  }
  {
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test1_3.svg");
    MolDraw2DSVG drawer(300, 300, outs);
    std::vector<int> highlights;
    highlights.push_back(0);
    highlights.push_back(4);
    highlights.push_back(5);
    drawer.drawMolecule(*m, &highlights);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test1_3.svg");
    delete m;
  }
  {
    std::string smiles = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test1_4.svg");
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawOptions().additionalAtomLabelPadding = 0.25;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test1_4.svg");
    delete m;
  }
  {
    // in this one, all three double bonds in the phenyl ring need to be inside
    // the aromatic ring.  There was a time when one of them strayed into the
    // aliphatic ring.
    std::string smiles =
        "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::unique_ptr<ROMol> romol(MolOps::removeAllHs(*m));
    RDDepict::compute2DCoords(*romol);
    WedgeMolBonds(*romol, &(romol->getConformer()));
    std::ofstream outs("test1_5.svg");
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawMolecule(*romol);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    // note that this hash check is likely to fail at the moment, due
    // to issue #4205.
    check_file_hash("test1_5.svg");
    delete m;
  }
  {
    // Here, the H should be between the two bonds off the N, not
    // on top of the vertical one.
    std::string smiles = "C[NH+](C)CCC";
    std::string nameBase = "test1_6";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs("test1_6.svg");
    outs << txt;
    outs.flush();
    outs.close();
    check_file_hash("test1_6.svg");
    delete m;
  }

  std::cout << " Done" << std::endl;
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <cairo.h>
#include "MolDraw2DCairo.h"
void test2() {
  std::cout << " ----------------- Test 2" << std::endl;
  {
    std::string smiles = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    MolDraw2DCairo drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    drawer.writeDrawingText("test2_1.png");
    check_file_hash("test2_1.png");
    delete m;
  }
  {
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    MolDraw2DCairo drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    std::string drawing = drawer.getDrawingText();
    TEST_ASSERT(drawing.size() > 0);
    std::ofstream ofs("test2_2.png");
    ofs.write(drawing.c_str(), drawing.size());
    check_file_hash("test2_2.png");
    delete m;
  }
  {
    // ensure we still work with a client-provided drawing context
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create(surface);

    MolDraw2DCairo drawer(300, 300, cr);
    std::vector<int> highlights;
    highlights.push_back(0);
    highlights.push_back(4);
    highlights.push_back(5);
    drawer.drawMolecule(*m, &highlights);
    drawer.finishDrawing();

    cairo_destroy(cr);
    cairo_surface_write_to_png(surface, "test2_3.png");
    cairo_surface_destroy(surface);
    check_file_hash("test2_3.png");

    delete m;
  }
  std::cout << " Done" << std::endl;
}
#else  // RDK_BUILD_CAIRO_SUPPORT
void test2() {}
#endif

void test3() {
  std::cout << " ----------------- Test 3" << std::endl;
  {
    std::string smiles = "C1CC1CC1ON1";
    std::string nameBase = "test3_1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    static const int ha[] = {0, 3, 4, 5};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, std::string> atomLabels;
    atomLabels[2] = "C1";
    atomLabels[1] = "a<sub>3</sub><sup>4</sup>";
    atomLabels[0] = "[CH2;X2:4]";
    atomLabels[6] = "[NH2+:7]";

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions().atomLabels = atomLabels;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();

      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().atomLabels = atomLabels;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }

  {
    std::string smiles = "C1CC1CC1ON1";
    std::string nameBase = "test3_2";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    static const int ha[] = {0, 3, 4, 5};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions().circleAtoms = false;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().circleAtoms = false;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    std::string nameBase = "test3_3";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    static const int ha[] = {11, 12, 13, 14, 15, 16};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, DrawColour> highlight_colors;
    highlight_colors[12] = DrawColour(0, 0, 1);
    highlight_colors[13] = DrawColour(0, 1, 0);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions().circleAtoms = true;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().circleAtoms = true;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    std::string nameBase = "test3_4";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    static const int ha[] = {11, 12, 13, 14, 15, 16, 3};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, DrawColour> highlight_colors;
    highlight_colors[12] = DrawColour(.5, .5, 1);
    highlight_colors[13] = DrawColour(.5, 1, .5);
    MolDrawOptions options;
    options.circleAtoms = true;
    options.highlightColour = DrawColour(1, .5, .5);
    options.continuousHighlight = true;

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles =
        "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1";
    std::string nameBase = "test3_5";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    static const int ha[] = {17, 18, 19, 20, 21, 6, 7, 8, 9, 31, 32};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, DrawColour> highlight_colors;
    MolDrawOptions options;
    options.circleAtoms = true;
    options.highlightColour = DrawColour(1, .5, .5);
    options.continuousHighlight = true;

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(200, 200, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }

  {
    std::string smiles =
        "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1";
    std::string nameBase = "test3_6";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    MolDrawOptions options;
    static const int ha1[] = {17, 18, 19, 20, 21};
    std::vector<int> highlight_atoms1(ha1, ha1 + sizeof(ha1) / sizeof(int));
    options.atomRegions.push_back(highlight_atoms1);
    static const int ha2[] = {6, 7, 8, 9, 31, 32};
    std::vector<int> highlight_atoms2(ha2, ha2 + sizeof(ha2) / sizeof(int));
    options.atomRegions.push_back(highlight_atoms2);
    options.includeAtomTags = true;

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles =
        "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1";
    std::string nameBase = "test3_7";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    MolDrawOptions options;
    options.continuousHighlight = true;
    static const int ha[] = {17, 20, 25};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, double> highlight_radii;
    highlight_radii[17] = 0.5;
    highlight_radii[20] = 1.0;
    std::map<int, DrawColour> highlight_colors;
    highlight_colors[17] = DrawColour(.5, .5, 1.);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors,
                          &highlight_radii);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors,
                          &highlight_radii);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }

  std::cout << " Done" << std::endl;
}

void test4() {
  std::cout << " ----------------- Test 4" << std::endl;
  {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/MolDraw2D/test_dir";
    fName += "/clash.mol";
    ROMol *m = MolFileToMol(fName);
    std::string nameBase = "test4_1";
    TEST_ASSERT(m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }

  std::cout << " Done" << std::endl;
}

void test5() {
  std::cout << " ----------------- Test 5" << std::endl;
  {
    std::string smiles = "*c1cc(*)cc(*)c1";
    std::string nameBase = "test5_1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDrawOptions options;
    options.dummiesAreAttachments = true;
    options.atomLabels[0] = "R1";
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles = "*C";
    std::string nameBase = "test5_2";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m, nullptr, true);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDrawOptions options;
    options.dummiesAreAttachments = true;
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles = "CC(F)(Cl)Br";
    std::string nameBase = "test5_3";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    m->getBondBetweenAtoms(1, 2)->setBondDir(Bond::UNKNOWN);
    RDDepict::compute2DCoords(*m, nullptr, true);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDrawOptions options;
    options.dummiesAreAttachments = true;
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  std::cout << " Done" << std::endl;
}

#ifdef RDK_TEST_MULTITHREADED
#include <thread>
#include <future>
namespace {
void runblock(const std::vector<ROMol *> &mols,
              const std::vector<std::string> &refData, unsigned int count,
              unsigned int idx) {
  for (unsigned int j = 0; j < 200; j++) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      ROMol *mol = mols[i];
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*mol);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      TEST_ASSERT(text == refData[i]);
    }
  }
}
}  // namespace

void testMultiThreaded() {
  std::cout << " ----------------- Test multi-threaded drawing" << std::endl;
  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  RDKit::SDMolSupplier suppl(fName);
  std::cerr << "reading molecules" << std::endl;
  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }

  std::cerr << "generating reference drawings" << std::endl;
  std::vector<std::string> refData(mols.size());
  for (unsigned int i = 0; i < mols.size(); ++i) {
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*(mols[i]));
    drawer.finishDrawing();
    refData[i] = drawer.getDrawingText();
    TEST_ASSERT(refData[i].find("<svg") != std::string::npos);
    TEST_ASSERT(refData[i].find("</svg>") != std::string::npos);
  }

  std::vector<std::future<void>> tg;
  unsigned int count = 4;
  std::cerr << "processing" << std::endl;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(
        std::async(std::launch::async, runblock, mols, refData, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  for (auto &&mol : mols) {
    delete mol;
  }
  std::cerr << " Done" << std::endl;
}
#else
void testMultiThreaded() {}
#endif

void test6() {
  std::cout << " ----------------- Test 6 (atom labels)" << std::endl;
  {
    std::string smiles = "CC[13CH2][CH2:7][CH-]C[15NH2+]C";
    std::string nameBase = "test6_1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs(nameBase + ".svg");
    outs << txt;
    // TEST_ASSERT(txt.find("<svg")!=std::string::npos);
    outs.flush();
    outs.close();
    outs.close();
    check_file_hash(nameBase + ".svg");
    delete m;
  }
  {
    auto m = "[C]1[C][C][CH][CH][CH]1"_smiles;
    std::string nameBase = "test6_2";
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    std::ofstream outs(nameBase + ".svg");
    outs << txt;
    outs.close();
    // start of bond-0
    TEST_ASSERT(
        txt.find("<path class='bond-0 atom-0 atom-1' d='M 273.6,147.5") !=
        std::string::npos);
    // start of first radical spot
    TEST_ASSERT(txt.find("<path d='M 286.5,143.5 L 286.5,143.4") !=
                std::string::npos);
    check_file_hash(nameBase + ".svg");
  }

  std::cerr << " Done" << std::endl;
}

void test7() {
  std::cout << " ----------------- Test 7 (backgrounds)" << std::endl;
  std::string smiles = "CCC";
  ROMol *m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  RDDepict::compute2DCoords(*m);
  {
    std::string nameBase = "test7_1";
    MolDraw2DSVG drawer(300, 300);
    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs((nameBase + ".svg").c_str());
    outs << txt;
    outs.flush();
    outs.close();
    check_file_hash(nameBase + ".svg");
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<rect") == std::string::npos);
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  {
    std::string nameBase = "test7_1";
    MolDraw2DCairo drawer(300, 300);
    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + ".png");
    check_file_hash(nameBase + ".png");
  }
#endif
  {
    std::string nameBase = "test7_2";
    MolDraw2DSVG drawer(300, 300);
    drawer.drawOptions().backgroundColour = DrawColour(.8, .8, .8);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs((nameBase + ".svg").c_str());
    outs << txt;
    outs.flush();
    outs.close();
    check_file_hash(nameBase + ".svg");
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<rect") != std::string::npos);
    TEST_ASSERT(txt.find("fill:#CCCCCC") != std::string::npos);
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  {
    std::string nameBase = "test7_2";
    MolDraw2DCairo drawer(300, 300);
    drawer.drawOptions().backgroundColour = DrawColour(.8, .8, .8);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + ".png");
    check_file_hash(nameBase + ".png");
  }
#endif
  delete m;
  std::cerr << " Done" << std::endl;
}

void test8PrepareMolForDrawing() {
  std::cout << " ----------------- Test8: PrepareMolDrawing" << std::endl;
  {
    std::string smiles = "c1ccccc1[C@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm);
      TEST_ASSERT(nm.getNumAtoms() == 9);  // this is a test for github #982
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(!nm.getConformer().is3D());
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() !=
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
      TEST_ASSERT(nm.getBondBetweenAtoms(6, 7)->getBondType() == Bond::SINGLE);
      TEST_ASSERT(nm.getBondBetweenAtoms(6, 7)->getBondDir() ==
                  Bond::BEGINWEDGE);

      // make sure we can do it again:
      MolDraw2DUtils::prepareMolForDrawing(nm);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() !=
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
    }
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm, false);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() ==
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
    }
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm, false, false);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() ==
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
    }
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm, false, true);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() ==
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
    }

    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm, true, true, false);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getBondType() !=
                  Bond::AROMATIC);
      TEST_ASSERT(nm.getBondBetweenAtoms(0, 1)->getIsAromatic());
      TEST_ASSERT(nm.getBondBetweenAtoms(6, 7)->getBondType() == Bond::SINGLE);
      TEST_ASSERT(nm.getBondBetweenAtoms(6, 7)->getBondDir() == Bond::NONE);
    }

    {
      // by default we don't force conformer generation
      RWMol nm(*m);
      RDDepict::compute2DCoords(nm);
      nm.getConformer().set3D(true);  // it's not really, we're cheating
      TEST_ASSERT(nm.getNumAtoms() == 9)
      MolDraw2DUtils::prepareMolForDrawing(nm);
      TEST_ASSERT(nm.getNumAtoms() == 9);
      TEST_ASSERT(nm.getNumConformers() == 1);  // we have a conformer anyway
      TEST_ASSERT(nm.getConformer().is3D());

      // but if we do force, it blows out that conformer:
      MolDraw2DUtils::prepareMolForDrawing(nm, true, true, true, true);
      TEST_ASSERT(!nm.getConformer().is3D());
    }

    delete m;
  }
  {
    std::string smiles = "C1CC[C@H]2NCCCC2C1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 10)
      MolDraw2DUtils::prepareMolForDrawing(nm);
      TEST_ASSERT(nm.getNumAtoms() == 11);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(!nm.getConformer().is3D());
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 10)->getBondType() == Bond::SINGLE);
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 10)->getBondDir() ==
                  Bond::BEGINDASH);

      // make sure we can do it again:
      MolDraw2DUtils::prepareMolForDrawing(nm);
      TEST_ASSERT(nm.getNumAtoms() == 11);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 10)->getBondType() == Bond::SINGLE);
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 10)->getBondDir() ==
                  Bond::BEGINDASH);
    }
    {
      RWMol nm(*m);
      TEST_ASSERT(nm.getNumAtoms() == 10)
      MolDraw2DUtils::prepareMolForDrawing(nm, false, false);
      TEST_ASSERT(nm.getNumAtoms() == 10);
      TEST_ASSERT(nm.getNumConformers() == 1);
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 2)->getBondType() == Bond::SINGLE);
      TEST_ASSERT(nm.getBondBetweenAtoms(3, 2)->getBondDir() ==
                  Bond::BEGINWEDGE);
    }
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub781() {
  std::cout << " ----------------- Test Github #781: Rendering single-atom "
               "molecules"
            << std::endl;

  {
    auto m = "C"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    // write the file so we can update the coords below more easily
    // if the font changes, for example.
    std::ofstream outs("testGithub781_1.svg");
    outs << txt;
    outs.close();
    check_file_hash("testGithub781_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // the start of the C
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 139.1 150.1") !=
                std::string::npos)
    // the start of the H
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 164.3 136.0") !=
                std::string::npos)
#else
    TEST_ASSERT(txt.find(">C</text>") != std::string::npos);
    TEST_ASSERT(txt.find(">H</text>") != std::string::npos);
    TEST_ASSERT(txt.find(">4</text>") != std::string::npos);
    TEST_ASSERT(txt.find("font-size:40px") != std::string::npos);
    TEST_ASSERT(txt.find("font-size:26px") != std::string::npos);
#endif
  }
  {
    auto m = "O"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    std::ofstream outs("testGithub781_2.svg");
    outs << txt;
    outs.close();
    check_file_hash("testGithub781_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // start of the H
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 98.9 136.0") !=
                std::string::npos);
    // start of the O
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 137.0 150.1") !=
                std::string::npos);
#else
    TEST_ASSERT(txt.find("<tspan>OH</tspan>") == std::string::npos);
#endif
  }
  {
    auto m = "[C]"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(600, 600);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    std::ofstream outs("testGithub781_3.svg");
    outs << txt;
    outs.close();
    check_file_hash("testGithub781_3.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // The C
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 289.1 300.1") !=
                std::string::npos);
    // the first radical marker
    TEST_ASSERT(txt.find("<path d='M 317.4,288.0 L 317.4,287.8") !=
                std::string::npos);
#else
    TEST_ASSERT(txt.find(">C</text>") != std::string::npos);
    // the first radical marker
    TEST_ASSERT(
        txt.find("<path d='M 318.4,288.0 L 318.4,287.8 L 318.3,287.7") !=
        std::string::npos);
#endif
  }
  {
    auto m = "C.CC.[Cl-]"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    std::ofstream outs("testGithub781_4.svg");
    outs << txt;
    outs.close();
    check_file_hash("testGithub781_4.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // start of C
    TEST_ASSERT(txt.find("<path  class='atom-0' d='M 27.4 198.0") !=
                std::string::npos);
    // start of l
    TEST_ASSERT(txt.find("<path  class='atom-3' d='M 36.4 74.8") !=
                std::string::npos);
#else
    TEST_ASSERT(txt.find(">C</text>") != std::string::npos);
    TEST_ASSERT(txt.find(">H</text>") != std::string::npos);
    TEST_ASSERT(txt.find(">l</text>") != std::string::npos);
#endif
  }
  {  // empty molecule
    auto *m = new ROMol();
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    // the Freetype version just draws the white canvas.
    std::ofstream outs("testGithub781_5.svg");
    outs << txt;
    outs.close();
    check_file_hash("testGithub781_5.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#else
    TEST_ASSERT(txt.find("<tspan>") == std::string::npos);
#endif
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub774() {
  std::cout << " ----------------- Test Github774" << std::endl;
  {
    std::string smiles =
        "Cc1c(C(=O)NCC[NH3+])[n+](=O)c2cc(CC[C@](F)(Cl)Br)ccc2n1[O-]";
    std::string nameBase = "test774_1";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolOps::Kekulize(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
      Point2D ocoords(1.0, 2.0);
      Point2D dcoords =
          drawer.getAtomCoords(std::make_pair(ocoords.x, ocoords.y));
      Point2D acoords = drawer.getDrawCoords(dcoords);
      TEST_ASSERT(feq(acoords.x, 1.0));
      TEST_ASSERT(feq(acoords.y, 2.0));
    }
    // m->setProp("_Name","mol");
    // std::cerr<<MolToMolBlock(*m)<<std::endl;
    delete m;
  }
  {
    std::string smiles =
        "CC(=O)\\C=C\\CC1[C@H]2N([C@@H](C(=O)O)C(C)(C)S2(=O)=O)C1=O";
    std::string nameBase = "test774_2";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolOps::Kekulize(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    // m->setProp("_Name","mol");
    // std::cerr<<MolToMolBlock(*m)<<std::endl;
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void test9MolLegends() {
  std::cout << " ----------------- Test 9 (molecule legends)" << std::endl;
  {
    auto m = "CC[13CH2][CH2:7][CH-]C[15NH2+]C"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m, "mol legend");
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs("test9_1.svg");
    outs << txt;
    // There's a bizarre thing whereby this file comes up as identical on my
    // MBP and Ubuntu 20.04 systems, but the hash codes are different.
    check_file_hash("test9_1.svg");
    // TEST_ASSERT(txt.find("<svg")!=std::string::npos);
  }
  std::cerr << " Done" << std::endl;
}

void testGithub852() {
  std::cout << " ----------------- Test Github852: Lines used to wedge bonds "
               "are too thick"
            << std::endl;
  {
    std::string smiles =
        "COc1cccc(NC(=O)[C@H](Cl)Sc2nc(ns2)c3ccccc3Cl)c1";  // made up
    std::string nameBase = "test852_1";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles =
        "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O";  // estradiol
    std::string nameBase = "test852_2";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    delete m;
  }
  {
    std::string smiles =
        "C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O";  // estradiol
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    {
      std::cerr << "----------------" << std::endl;
      std::string nameBase = "test852_2a";
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(200, 200, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    {
      std::cerr << "----------------" << std::endl;
      std::string nameBase = "test852_2b";
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(250, 250, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    {
      std::cerr << "----------------" << std::endl;
      std::string nameBase = "test852_2c";
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(400, 400, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
    {
      std::cerr << "----------------" << std::endl;
      std::string nameBase = "test852_2d";
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(500, 500, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }

    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub860() {
  std::cout << " ----------------- Test Github860: Atom symbols in wrong order "
               "if bond comes from right"
            << std::endl;
  {
    auto m = "[15NH3+:1]-C#C-[15NH3+:2]"_smiles;
    std::string nameBase = "test860_1";
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
  }
  {
    auto m = "[15NH3+:1]-C#C-C([15NH3+:2])([15NH3+:3])-C#C-[15NH3+:4]"_smiles;
    std::string nameBase = "test860_2";
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
  }
  {
    auto m = "[15NH3+:1]-CCCCCCCC-[15NH3+:4]"_smiles;
    std::string nameBase = "test860_3";
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
      check_file_hash(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
  }
  std::cerr << " Done" << std::endl;
}

void testGithub910() {
  std::cout << " ----------------- Test Github #910: ugly coordinates "
               "generated for peptide chains"
            << std::endl;
  // this really isn't much of a test, but it does help visually confirm that
  // things are actually ok
  {
    // this is a ChEMBL molecule
    std::string smiles =
        "CSCC[C@H](NC(=O)[C@@H](CCC(N)=O)NC(=O)[C@@H](N)Cc1c[nH]c2ccccc12)C(="
        "O)"
        "NCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CO)C(=O)O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test910_1.svg");
    MolDraw2DSVG drawer(600, 300, outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test910_1.svg");
    delete m;
  }
  {  // now with Hs
    // this is a ChEMBL molecule
    std::string smiles =
        "CSCC[C@H](NC(=O)[C@@H](CCC(N)=O)NC(=O)[C@@H](N)Cc1c[nH]c2ccccc12)C(="
        "O)"
        "NCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CO)C(=O)O";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    std::ofstream outs("test910_2.svg");
    MolDraw2DSVG drawer(600, 300, outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test910_2.svg");
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub932() {
  std::cout << " ----------------- Test Github #932: mistake in SVG for "
               "wedged bonds"
            << std::endl;
  {
    std::string smiles = "CC[C@](F)(Cl)Br";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    TEST_ASSERT(text.find("evenoddstroke") == std::string::npos);
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub953() {
  std::cout << " ----------------- Test Github #953: default color should "
               "not be cyan"
            << std::endl;
  {
    std::string smiles = "[Nb]";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    TEST_ASSERT(text.find("#00FFFF") == std::string::npos);
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub983() {
  std::cout << " ----------------- Test Github #983: wedged bonds between "
               "chiral centers drawn improperly"
            << std::endl;
  {
    // this has an ugly drawing (wedged bond between chiral centers) but we
    // force it to be drawn that way just to check.
    std::string mb =
        "\n\
  Mrv1561 07241608122D\n\
\n\
  6  5  0  0  0  0            999 V2000\n\
    8.6830   -9.5982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    9.3975   -9.1857    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n\
   10.1120   -9.5982    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n\
    9.3975   -8.3607    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n\
   10.8264   -9.1857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   10.1120  -10.4232    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
  3  5  1  0  0  0  0\n\
  3  2  1  1  0  0  0\n\
  2  4  1  1  0  0  0\n\
  3  6  1  0  0  0  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb, false, false);
    TEST_ASSERT(m);
    MolOps::sanitizeMol(*m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test983_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test983_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    TEST_ASSERT(text.find("<path class='bond-1 atom-2 atom-4' d='M 126.9,116.0"
                          " L 184.0,90.8 L 177.2,79.1 Z'"
                          " style='fill:#000000;") != std::string::npos);
#else
    TEST_ASSERT(text.find("<path class='bond-1 atom-2 atom-4' d='M 126.5,111.6"
                          " L 182.7,86.9 L 176.0,75.3 Z'"
                          " style='fill:#000000;fill-rule:evenodd;"
                          "fill-opacity:1;stroke:#000000;") !=
                std::string::npos);
#endif
    delete m;
  }
  {
    std::string mb =
        "\n\
  Mrv1561 07241616282D\n\
\n\
 12 12  0  0  1  0            999 V2000\n\
   10.4656   -7.9623    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
    9.7496   -8.3748    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n\
    8.9075   -9.4746    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n\
    7.5671   -9.4746    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    8.2373   -8.9934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    8.6497  -10.2651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    9.0392   -7.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    7.8249  -10.2651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    7.1547  -10.1792    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
    6.8567   -9.0622    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n\
   10.3338   -8.9591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
    8.6841   -8.6669    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n\
  2  1  1  0  0  0  0\n\
  3  2  1  0  0  0  0\n\
  4  5  1  0  0  0  0\n\
  5  3  1  0  0  0  0\n\
  6  3  1  0  0  0  0\n\
  7  2  1  0  0  0  0\n\
  8  6  1  0  0  0  0\n\
  9  4  1  0  0  0  0\n\
 10  4  1  0  0  0  0\n\
  2 11  1  6  0  0  0\n\
  3 12  1  6  0  0  0\n\
  8  4  1  0  0  0  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 1)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 4)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 4)->getBondDir() == Bond::BEGINWEDGE);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test983_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test983_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    TEST_ASSERT(text.find("<path class='bond-3 atom-2 atom-4' d='M 103.7,117.6"
                          " L 76.9,93.5 L 72.3,99.8 Z'"
                          " style='fill:#000000;") != std::string::npos);
#else
    TEST_ASSERT(text.find("<path class='bond-3 atom-2 atom-4' d='M 105.1,114.8"
                          " L 78.5,90.9 L 74.0,97.2 Z'"
                          " style='fill:#000000;fill-rule:evenodd;"
                          "fill-opacity:1;stroke:#000000;") !=
                std::string::npos);
#endif

    MolDraw2DUtils::prepareMolForDrawing(*m);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 1)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 4)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 4)->getBondDir() == Bond::BEGINWEDGE);

    RWMol nm(*m);
    MolDraw2DUtils::prepareMolForDrawing(nm);
    TEST_ASSERT(nm.getBondBetweenAtoms(2, 1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(nm.getBondBetweenAtoms(2, 1)->getBondDir() == Bond::NONE);
    TEST_ASSERT(nm.getBondBetweenAtoms(2, 4)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(nm.getBondBetweenAtoms(2, 4)->getBondDir() == Bond::BEGINWEDGE);

    delete m;
  }

  std::cerr << " Done" << std::endl;
}

void testDeuteriumTritium() {
  std::cout << " ----------------- Test Deuterium, Tritium" << std::endl;
  {
    auto m = "C([2H])([2H])([2H])[2H]"_smiles;
    RDDepict::compute2DCoords(*m);
    std::string nameBase = "testNoDeuterium";
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawOptions().atomLabelDeuteriumTritium = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.close();
    std::ifstream ins((nameBase + ".svg").c_str());
    bool ok = true;
    unsigned int count = 0;
    while (ok) {
      std::string line;
      std::getline(ins, line);
      ok = (ins.good() && !ins.eof());
      if (!ok) {
        continue;
      }
        // there are no characters to look for, but each atom should
        // be made of 2 glyphs, the superscript 2 and the H.
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      if ((line.find("atom-") != std::string::npos)) {
        if ((line.find("bond-") == std::string::npos)) {
          ++count;
        }
      }
#else
      // a bit kludgy, but...
      if (line.find("<text x='250.5' y='152.7' class='atom-1'"
                    " style='font-size:26px;font-style:normal;"
                    "font-weight:normal;fill-opacity:1;stroke:none;"
                    "font-family:sans-serif;text-anchor:start;"
                    "fill:#000000' >2</text>") != std::string::npos) {
        ++count;
      }
#endif
    }
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    TEST_ASSERT(count == 8);
#else
    // the first superscript 2
    TEST_ASSERT(count == 1);
#endif
    check_file_hash(nameBase + ".svg");
  }
  {
    auto m = "C([3H])([3H])([3H])[3H]"_smiles;
    RDDepict::compute2DCoords(*m);
    std::string nameBase = "testNoTritium";
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawOptions().atomLabelDeuteriumTritium = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.close();
    std::ifstream ins((nameBase + ".svg").c_str());
    bool ok = true;
    unsigned int count = 0;
    while (ok) {
      std::string line;
      std::getline(ins, line);
      ok = (ins.good() && !ins.eof());
      if (!ok) {
        continue;
      }
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // there are no characters to look for, but each atom should
      // be made of 2 glyphs, the superscript 3 and the H.
      if ((line.find("atom-") != std::string::npos)) {
        if ((line.find("bond-") == std::string::npos)) {
          ++count;
        }
      }
#else
      if (line.find("<text x='250.5' y='152.7' class='atom-1'"
                    " style='font-size:26px;font-style:normal;"
                    "font-weight:normal;fill-opacity:1;stroke:none;"
                    "font-family:sans-serif;text-anchor:start;"
                    "fill:#000000' >3</text>") != std::string::npos) {
        ++count;
      }
#endif
    }
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    TEST_ASSERT(count == 8);
#else
    TEST_ASSERT(count == 1);
#endif
    check_file_hash(nameBase + ".svg");
  }
  {
    auto m = "C([2H])([2H])([2H])[2H]"_smiles;
    RDDepict::compute2DCoords(*m);
    std::string nameBase = "testDeuterium";
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawOptions().atomLabelDeuteriumTritium = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.close();
    std::ifstream ins((nameBase + ".svg").c_str());
    bool ok = true;
    unsigned int count = 0;
    while (ok) {
      std::string line;
      std::getline(ins, line);
      ok = (ins.good() && !ins.eof());
      if (!ok) {
        continue;
      }
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // there should be just 1 glyph per atom - a D
      if ((line.find("atom-") != std::string::npos)) {
        if ((line.find("bond-") == std::string::npos)) {
          ++count;
        }
      }
#else
      if ((line.find("baseline-shift:super") == std::string::npos) &&
          (line.find(">2<") == std::string::npos) &&
          (line.find(">D<") != std::string::npos)) {
        ++count;
      }
#endif
    }
    TEST_ASSERT(count == 4);
    check_file_hash(nameBase + ".svg");
  }
  {
    auto m = "C([3H])([3H])([3H])[3H]"_smiles;
    RDDepict::compute2DCoords(*m);
    std::string nameBase = "testTritium";
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawOptions().atomLabelDeuteriumTritium = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.close();
    std::ifstream ins((nameBase + ".svg").c_str());
    bool ok = true;
    unsigned int count = 0;
    while (ok) {
      std::string line;
      std::getline(ins, line);
      ok = (ins.good() && !ins.eof());
      if (!ok) {
        continue;
      }
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // there should be just 1 glyph per atom - a T
      if ((line.find("atom-") != std::string::npos)) {
        if ((line.find("bond-") == std::string::npos)) {
          ++count;
        }
      }
#else
      if ((line.find("baseline-shift:super") == std::string::npos) &&
          (line.find(">3<") == std::string::npos) &&
          (line.find(">T<") != std::string::npos)) {
        ++count;
      }
#endif
    }
    TEST_ASSERT(count == 4);
    check_file_hash(nameBase + ".svg");
  }
  std::cerr << " Done" << std::endl;
}

void testCrossedBonds() {
  std::cerr << " ----------------- Test crossed bonds" << std::endl;
  {
    std::string smiles = "CC=CC";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    m->getBondWithIdx(1)->setBondDir(Bond::EITHERDOUBLE);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string nameBase = "crossed_bonds";
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(300, 300, outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.close();
    check_file_hash(nameBase + ".svg");
    delete m;
  }
  std::cerr << " Done" << std::endl;
}
void test10DrawSecondMol() {
  std::cout << " ----------------- Testing drawing a second molecule"
            << std::endl;
  std::string mb1 =
      "\n\
  Mrv1561 08301611102D\n\
\n\
  3  2  0  0  0  0            999 V2000\n\
   -2.5670    1.3616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.8525    1.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.1380    1.3616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
M  END";

  RWMol *m1 = MolBlockToMol(mb1);
  TEST_ASSERT(m1);
  MolOps::sanitizeMol(*m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  RDGeom::Point3D c1 = MolTransforms::computeCentroid(m1->getConformer());
  for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
    RDGeom::Point3D &p = m1->getConformer().getAtomPos(i);
    p -= c1;
  }
  std::string mb2 =
      "\n\
  Mrv1561 08301611122D\n\
\n\
  3  2  0  0  0  0            999 V2000\n\
   -1.9900    2.2136    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.5775    1.4991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -1.9900    0.7846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0  0  0  0\n\
  2  3  1  0  0  0  0\n\
M  END";
  RWMol *m2 = MolBlockToMol(mb2);
  TEST_ASSERT(m2);
  MolOps::sanitizeMol(*m2);
  MolDraw2DUtils::prepareMolForDrawing(*m2);
  RDGeom::Point3D c2 = MolTransforms::computeCentroid(m2->getConformer());
  for (unsigned int i = 0; i < m2->getNumAtoms(); ++i) {
    RDGeom::Point3D &p = m2->getConformer().getAtomPos(i);
    p -= c2;
  }

  {
    MolDraw2DSVG drawer(200, 200);
    drawer.drawOptions().padding = 0.2;
    drawer.drawMolecule(*m1);
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_1.svg");
  }
  {
    MolDraw2DSVG drawer(200, 200);
    drawer.drawOptions().padding = 0.2;
    drawer.drawMolecule(*m2);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_2.svg");
  }
  {
    MolDraw2DSVG drawer(400, 200, 200, 200);
    drawer.drawOptions().padding = 0.2;
    drawer.drawMolecule(*m1);
    drawer.setOffset(200, 0);
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_3.svg");
  }
  {
    MolDraw2DSVG drawer(200, 400, 200, 200);
    drawer.drawOptions().padding = 0.2;
    drawer.drawMolecule(*m1);
    drawer.setOffset(0, 200);
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_4.svg");
  }
  {
    MolDraw2DSVG drawer(200, 400, 200, 200);
    Point2D minv(1000, 1000);
    Point2D maxv(-1000, -1000);
    for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
      const RDGeom::Point3D &pti = m1->getConformer().getAtomPos(i);
      minv.x = std::min(minv.x, pti.x);
      minv.y = std::min(minv.y, pti.y);
      maxv.x = std::max(maxv.x, pti.x);
      maxv.y = std::max(maxv.y, pti.y);
    }
    drawer.setScale(200, 200, minv, maxv);
    drawer.drawMolecule(*m1);
    drawer.setOffset(0, 200);
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_5.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_5.svg");
  }
  {
    MolDraw2DSVG drawer(200, 400, 200, 200);
    Point2D minv(1000, 1000);
    Point2D maxv(-1000, -1000);
    for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
      const RDGeom::Point3D &pti = m1->getConformer().getAtomPos(i);
      minv.x = std::min(minv.x, pti.x);
      minv.y = std::min(minv.y, pti.y);
      maxv.x = std::max(maxv.x, pti.x);
      maxv.y = std::max(maxv.y, pti.y);
    }
    drawer.drawOptions().padding = 0.2;
    drawer.setScale(200, 200, minv, maxv);
    drawer.drawMolecule(*m1);
    drawer.setOffset(0, 200);
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test10_6.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test10_6.svg");
  }

  delete m1;
  delete m2;
  std::cerr << " Done" << std::endl;
}

void test11DrawMolGrid() {
  std::cout << " ----------------- Testing drawing a grid of molecules"
            << std::endl;

  std::string smiles =
      "COc1cccc(NC(=O)[C@H](Cl)Sc2nc(ns2)c3ccccc3Cl)c1";  // made up
  RWMol *m1 = SmilesToMol(smiles);
  TEST_ASSERT(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  RDGeom::Point3D c1 = MolTransforms::computeCentroid(m1->getConformer());
  for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
    RDGeom::Point3D &p = m1->getConformer().getAtomPos(i);
    p -= c1;
  }
  smiles = "NC(=O)[C@H](Cl)Sc1ncns1";  // made up
  RWMol *m2 = SmilesToMol(smiles);
  TEST_ASSERT(m2);
  MolDraw2DUtils::prepareMolForDrawing(*m2);
  RDGeom::Point3D c2 = MolTransforms::computeCentroid(m2->getConformer());
  for (unsigned int i = 0; i < m2->getNumAtoms(); ++i) {
    RDGeom::Point3D &p = m2->getConformer().getAtomPos(i);
    p -= c2;
  }
  smiles = "BrCNC(=O)[C@H](Cl)Sc1ncns1";  // made up
  RWMol *m3 = SmilesToMol(smiles);
  TEST_ASSERT(m3);
  MolDraw2DUtils::prepareMolForDrawing(*m3);
  RDGeom::Point3D c3 = MolTransforms::computeCentroid(m3->getConformer());
  for (unsigned int i = 0; i < m3->getNumAtoms(); ++i) {
    RDGeom::Point3D &p = m3->getConformer().getAtomPos(i);
    p -= c3;
  }

  {
    MolDraw2DSVG drawer(500, 400, 250, 200);
    drawer.drawMolecule(*m1, "m1");
    drawer.setOffset(250, 0);
    drawer.drawMolecule(*m2, "m2");
    drawer.setOffset(0, 200);
    drawer.drawMolecule(*m3, "m3");
    drawer.setOffset(250, 200);
    drawer.drawMolecule(*m1, "m4");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test11_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test11_1.svg");
  }
  {  // drawing "out of order"
    MolDraw2DSVG drawer(500, 400, 250, 200);
    drawer.setOffset(250, 0);
    drawer.drawMolecule(*m1, "m1");
    drawer.setOffset(0, 0);
    drawer.drawMolecule(*m2, "m2");
    drawer.setOffset(0, 200);
    drawer.drawMolecule(*m1, "m3");
    drawer.setOffset(250, 200);
    drawer.drawMolecule(*m2, "m4");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test11_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test11_2.svg");
  }
  delete m1;
  delete m2;
  delete m3;
  std::cerr << " Done" << std::endl;
}

void test12DrawMols() {
  std::cout << " ----------------- Testing drawMolecules" << std::endl;

  auto setup_mol = [](const std::string &smi, const std::string leg,
                      std::vector<ROMol *> &mols,
                      std::vector<std::string> &legends) {
    mols.push_back(SmilesToMol(smi));
    TEST_ASSERT(mols.back());
    legends.push_back(leg);
  };
  std::vector<ROMol *> mols;
  std::unique_ptr<std::vector<std::string>> legends(
      new std::vector<std::string>());
  // made up SMILES, each with sequence F, Cl, Br so we can see which
  // ones are drawn, which ones are missing.
  setup_mol("COc1cccc(NC(=O)[C@H](F)Sc2nc(ns2)c3ccccc3F)c1", "m1", mols,
            *legends);
  setup_mol("NC(=O)[C@H](F)Sc1ncns1", "m2", mols, *legends);
  setup_mol("COc1cccc(NC(=O)[C@H](Cl)Sc2nc(ns2)c3ccccc3F)c1", "m3", mols,
            *legends);
  setup_mol("NC(=O)[C@H](Cl)Sc1ncns1", "m4", mols, *legends);
  setup_mol("COc1cccc(NC(=O)[C@H](Br)Sc2nc(ns2)c3ccccc3F)c1", "m5", mols,
            *legends);
  setup_mol("NC(=O)[C@H](Br)Sc1ncns1", "m6", mols, *legends);

  {
    MolDraw2DSVG drawer(750, 400, 250, 200);
    drawer.drawMolecules(mols, legends.get());
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test12_1.svg");
  }

  {  // github #1325: multiple molecules in one pane
    MolDraw2DSVG drawer(300, 300, 300, 300);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test12_3.svg");
  }

  {  // github #1325: multiple molecules in one pane
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test12_4.svg");
  }
  {
    delete mols[2];
    delete mols[4];
    mols[2] = nullptr;
    mols[4] = nullptr;
    MolDraw2DSVG drawer(750, 400, 250, 200);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test12_2.svg");
  }
  for (auto m : mols) {
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void test13JSONConfig() {
  std::cerr << " ----------------- Test JSON Configuration" << std::endl;
  {
    auto m = "CCO"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    const char *json =
        "{\"legendColour\":[1.0,0.5,1.0], \"rotate\": 90, "
        "\"bondLineWidth\": 5}";
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, json);
    drawer.drawMolecule(*m, "foo");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test13_1.svg");
    outs << text;
    outs.close();
    check_file_hash("test13_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // we'll just have to assume that this pink is for the legend
    TEST_ASSERT(text.find("' fill='#FF7FFF") != std::string::npos);
    TEST_ASSERT(text.find("<path class='bond-0 atom-0 atom-1' d='M 119.4,8.2"
                          " L 162.9,83.6'") != std::string::npos);
#else
    TEST_ASSERT(text.find("sans-serif;text-anchor:start;fill:#FF7FFF") !=
                std::string::npos);
    TEST_ASSERT(text.find("<path class='bond-0 atom-0 atom-1' d='M 119.8,8.2"
                          " L 162.1,81.5'") != std::string::npos);
#endif
    // these days the bond line width scales with the rest of the
    // drawing, and at this size this comes out as 6px.
    TEST_ASSERT(text.find("stroke-width:5.0px") != std::string::npos);
  }
  std::cerr << " Done" << std::endl;
}

void testGithub1090() {
  std::cout << " ----------------- Testing github 1090: escape html characters "
               "in SVG output"
            << std::endl;

  std::string smiles = "CCOC";  // made up
  RWMol *m1 = SmilesToMol(smiles);
  TEST_ASSERT(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  {
    ROMol lm(*m1);
    MolDraw2DSVG drawer(250, 200);
    drawer.drawOptions().atomLabels[0] = "C&1";
    drawer.drawOptions().atomLabels[1] = "[CH2<1]";
    drawer.drawOptions().atomLabels[3] = "[C>1H3]";
    drawer.drawMolecule(lm, "legend&legend>1");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub1090_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub1090_1.svg");
    TEST_ASSERT(text.find("C&1") == std::string::npos);
    TEST_ASSERT(text.find("<<") == std::string::npos);
    TEST_ASSERT(text.find(">>") == std::string::npos);
    TEST_ASSERT(text.find("d&l") == std::string::npos);
  }
  delete m1;
  std::cerr << " Done" << std::endl;
}

void testGithub1035() {
  std::cout << " ----------------- Testing github 1035: overflow bug in SVG "
               "color generation"
            << std::endl;

  std::string smiles = "CCOC";  // made up
  RWMol *m1 = SmilesToMol(smiles);
  TEST_ASSERT(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  std::vector<int> highlights;
  highlights.push_back(0);
  highlights.push_back(1);
  {
    MolDraw2DSVG drawer(250, 200);
    drawer.drawOptions().highlightColour = DrawColour(1.1, .5, .5);
    bool ok = false;
    try {
      drawer.drawMolecule(*m1, &highlights);
    } catch (const ValueErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    MolDraw2DSVG drawer(250, 200);
    drawer.drawOptions().highlightColour = DrawColour(.1, -.5, .5);
    bool ok = false;
    try {
      drawer.drawMolecule(*m1, &highlights);
    } catch (const ValueErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    MolDraw2DSVG drawer(250, 200);
    drawer.drawOptions().highlightColour = DrawColour(1., .5, 1.5);
    bool ok = false;
    try {
      drawer.drawMolecule(*m1, &highlights);
    } catch (const ValueErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  delete m1;
  std::cerr << " Done" << std::endl;
}

void testGithub1271() {
  std::cout << " ----------------- Testing github 1271: MolDraw2D not drawing "
               "anything for molecules aligned with the X or Y axes"
            << std::endl;
  {
    std::string mb =
        "ethane\n\
     RDKit          2D\n\
\n\
  2  1  0  0  0  0  0  0  0  0999 V2000\n\
   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test1271_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test1271_1.svg");
    TEST_ASSERT(text.find("d='M 0,200 0,200") == std::string::npos);
    delete m;
  }
  {
    std::string mb =
        "ethane\n\
     RDKit          2D\n\
\n\
  2  1  0  0  0  0  0  0  0  0999 V2000\n\
   -0.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.0000   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test1271_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test1271_2.svg");
    TEST_ASSERT(text.find("d='M 0,200 0,200") == std::string::npos);
    delete m;
  }
  {
    std::string mb =
        "water\n\
     RDKit          2D\n\
\n\
  1  0  0  0  0  0  0  0  0  0999 V2000\n\
   -0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test1271_3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test1271_3.svg");
    TEST_ASSERT(text.find("d='M 0,200 0,200") == std::string::npos);
    delete m;
  }
  {
    std::string mb =
        "water\n\
     RDKit          2D\n\
\n\
  1  0  0  0  0  0  0  0  0  0999 V2000\n\
   -0.0000    0.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
M  END";
    RWMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test1271_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test1271_4.svg");
    TEST_ASSERT(text.find("d='M 0,200 0,200") == std::string::npos);
    delete m;
  }

  {
    std::string smiles = "C=C(O)C(O)";  // made up
    RWMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    smiles = "O";
    RWMol *m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    MolDraw2DUtils::prepareMolForDrawing(*m2);

    MolDraw2DSVG drawer(500, 200, 250, 200);
    drawer.drawMolecule(*m1, "m1");
    drawer.setOffset(250, 0);
    drawer.drawMolecule(*m2, "m2");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test1271_5.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test1271_5.svg");
    delete m1;
    delete m2;
  }

  std::cerr << " Done" << std::endl;
}

void testGithub1322() {
  std::cout << " ----------------- Testing github 1322: add custom atom labels"
            << std::endl;
  {
    auto m1 = "CCC[Se]"_smiles;  // made up
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test1322_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test1322_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // there should be 2 paths of class atom-3, one for the S,
      // the other for the e + plus one bond = 3.
      size_t start_pos = 0;
      int count = 0;
      while (true) {
        start_pos = text.find("atom-3", start_pos);
        if (start_pos == std::string::npos) {
          break;
        }
        ++count;
        ++start_pos;
      }
      TEST_ASSERT(count == 3);
#else
      TEST_ASSERT(text.find(">S</text>") != std::string::npos);
      TEST_ASSERT(text.find(">e</text>") != std::string::npos);
#endif
    }
    {
      m1->getAtomWithIdx(3)->setProp(common_properties::atomLabel,
                                     "customlabel");
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test1322_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test1322_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // there should be 11 paths of class atom-3, one for each letter
      // of customlabel + plus one bond = 12.
      size_t start_pos = 0;
      int count = 0;
      while (true) {
        start_pos = text.find("atom-3", start_pos);
        if (start_pos == std::string::npos) {
          break;
        }
        ++count;
        ++start_pos;
      }
      TEST_ASSERT(count == 12);
#else
      TEST_ASSERT(text.find(">S</text>") == std::string::npos);
      TEST_ASSERT(text.find(">s</text>") != std::string::npos);
      TEST_ASSERT(text.find(">b</text>") != std::string::npos);
#endif
    }
  }
  std::cerr << " Done" << std::endl;
}

void test14BWPalette() {
  std::cout << " ----------------- Testing use of a black & white palette"
            << std::endl;
  {
    std::string smiles = "CNC(Cl)C(=O)O";
    RWMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    {  // start with color
      MolDraw2DSVG drawer(200, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();

      TEST_ASSERT(text.find("stroke:#00CC00") != std::string::npos);

      std::ofstream outs("test14_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test14_1.svg");
    }
    {  // now B&W
      MolDraw2DSVG drawer(200, 200);
      assignBWPalette(drawer.drawOptions().atomColourPalette);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();

      TEST_ASSERT(text.find("stroke:#00CC00") == std::string::npos);

      std::ofstream outs("test14_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test14_2.svg");
    }
    delete m1;
  }
  std::cerr << " Done" << std::endl;
}

void test15ContinuousHighlightingWithGrid() {
  std::cerr << " ----------------- Testing use of continuous highlighting with "
               "drawMolecules"
            << std::endl;

  {
    std::string smiles =
        "COc1cccc(NC(=O)[C@H](Cl)Sc2nc(ns2)c3ccccc3Cl)c1";  // made up
    RWMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    smiles = "NC(=O)[C@H](Cl)Sc1ncns1";  // made up
    RWMol *m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    std::vector<ROMol *> mols;
    mols.push_back(m1);
    mols.push_back(m2);
    std::vector<std::vector<int>> atHighlights(2);
    atHighlights[0].push_back(0);
    atHighlights[0].push_back(1);
    atHighlights[0].push_back(2);
    atHighlights[0].push_back(6);

    atHighlights[1].push_back(0);
    atHighlights[1].push_back(1);
    atHighlights[1].push_back(2);
    atHighlights[1].push_back(6);

    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().continuousHighlight = false;
      drawer.drawMolecules(mols, nullptr, &atHighlights);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test15_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test15_1.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4.0px;") ==
                  std::string::npos);
    }

    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().continuousHighlight = true;
      drawer.drawMolecules(mols, nullptr, &atHighlights);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test15_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test15_2.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4.3px") !=
                  std::string::npos);
    }
    for (auto &&mol : mols) {
      delete mol;
    }
  }

  std::cerr << " Done" << std::endl;
}

void testGithub1829() {
  std::cerr << " ----------------- Testing github 1829: crash when "
               "drawMolecules() is called with an empty list"
            << std::endl;
  {
    std::vector<ROMol *> mols;
    MolDraw2DSVG drawer(750, 400, 250, 200);
    // this should run quietly without complaining
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
  }
  std::cerr << " Done" << std::endl;
}

void test16MoleculeMetadata() {
  std::cout << " ----------------- Testing inclusion of molecule metadata"
            << std::endl;
  {
    std::string smiles = "CN[C@H](Cl)C(=O)O";
    std::unique_ptr<RWMol> m1(SmilesToMol(smiles));
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    {  // one molecule
      MolDraw2DSVG drawer(200, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.addMoleculeMetadata(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test16_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test16_1.svg");
      TEST_ASSERT(text.find("idx=\"2\" atom-smiles=\"[NH]\" drawing-x=\"55.") !=
                  std::string::npos);
      TEST_ASSERT(text.find("idx=\"2\" begin-atom-idx=\"2\" end-atom-idx=\"3\" "
                            "bond-smiles=\"-\"") != std::string::npos);
    }

    {  // multiple molecules
      MolDraw2DSVG drawer(400, 400, 200, 200);
      auto *rom = rdcast<ROMol *>(m1.get());
      std::vector<ROMol *> ms = {new ROMol(*rom), new ROMol(*rom),
                                 new ROMol(*rom), new ROMol(*rom)};
      drawer.drawMolecules(ms);
      drawer.addMoleculeMetadata(ms);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test16_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("test16_2.svg");

      TEST_ASSERT(text.find("idx=\"2\" atom-smiles=\"[NH]\" drawing-x=\"55.") !=
                  std::string::npos);
      TEST_ASSERT(
          text.find("idx=\"2\" atom-smiles=\"[NH]\" drawing-x=\"255.") !=
          std::string::npos);

      for (auto ptr : ms) {
        delete ptr;
      }
    }
  }

  std::cerr << " Done" << std::endl;
}

void test17MaxMinFontSize() {
  std::cout << " ----------------- Test 17 - Testing maximum font size"
            << std::endl;
  {
    auto m = R"CTAB(
  Mrv2014 03142110062D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.4137 4.3863 0 0
M  V30 2 C -5.7475 3.6163 0 0
M  V30 3 C -5.7475 2.0763 0 0
M  V30 4 C -4.4137 1.3063 0 0
M  V30 5 N -3.08 2.0763 0 0
M  V30 6 C -3.0801 3.6163 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 6
M  V30 3 2 2 3
M  V30 4 1 3 4
M  V30 5 2 4 5
M  V30 6 1 5 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    TEST_ASSERT(m);
    std::string nameBase = "test17_";

    {
      std::ofstream outs((nameBase + "1.svg").c_str());
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // where it starts drawing the N is a poor surrogate for checking
      // the font size, but all we have.
      TEST_ASSERT(text.find("<path  class='atom-4' d='M 257.0 204.0") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:40px") != std::string::npos);
#endif
    }
    {
      std::ofstream outs((nameBase + "2.svg").c_str());
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().maxFontSize = -1;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // where it starts drawing the N
      TEST_ASSERT(text.find("<path  class='atom-4' d='M 252.7 199.4") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:53px") != std::string::npos);
#endif
    }
    {
      std::ofstream outs((nameBase + "3.svg").c_str());
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().maxFontSize = 20;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "3.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // where it starts drawing the N
      TEST_ASSERT(text.find("<path  class='atom-4' d='M 262.5 211.1") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:20px") != std::string::npos);
#endif
    }
    {
      auto m1 =
          "C[C@H](C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3"
          "=CN(N=C3)C4CCNCC4)N"_smiles;
      std::ofstream outs((nameBase + "4.svg").c_str());
      MolDraw2DSVG drawer(200, 200);
      // this is currently the default min font size.  Repeated for
      // documentation of test.
      drawer.drawOptions().minFontSize = 12;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "4.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-8' d='M 164.3 92.8") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:12px") != std::string::npos);
#endif
    }
  }

  std::cerr << " Done" << std::endl;
}

void test18FixedScales() {
  std::cout << " ----------------- Testing use of fixed scales for drawing."
            << std::endl;
  std::string nameBase = "test18_";
  {
    auto m = "Clc1ccccc1"_smiles;
    TEST_ASSERT(m);
    {
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "1.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // where it starts drawing the l is a poor surrogate for checking
      // the font size, but all we have.
      TEST_ASSERT(text.find("<path  class='atom-0' d='M 283.2 137.0") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:33px") != std::string::npos);
#endif
    }
    {
      MolDraw2DSVG drawer(300, 300);
      // fix scale so bond is 5% if window width.
      drawer.drawOptions().fixedScale = 0.05;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "2.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // where it starts drawing the l.
      TEST_ASSERT(text.find("<path  class='atom-0' d='M 185.0 146.6") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:9px") != std::string::npos);
#endif
    }
  }
  {
    std::string smi =
        "C[C@@H](N[C@@H]1CC[C@@H](C(=O)N2CCC(C(=O)N3CCCC3)"
        "(c3ccccc3)CC2)C(C)(C)C1)c1ccc(Cl)cc1";
    std::unique_ptr<ROMol> m(SmilesToMol(smi));
    TEST_ASSERT(m);

    {
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "3.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "3.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-2' d='M 72.1 191.7") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:10px") != std::string::npos);
#endif
    }
    {
      // fix bond length to 10 pixels.
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().fixedBondLength = 10;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "4.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "4.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-2' d='M 104.0 170.3") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:6px") != std::string::npos);
#endif
    }
    {
      // this one should be the same size as the first (_3), as it won't scale
      // up if the picture won't fit.
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().fixedBondLength = 30;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "5.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "5.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-2' d='M 73.3 183.9") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("font-size:10px") != std::string::npos);
#endif
    }
  }
  std::cerr << " Done" << std::endl;
}

void test19RotateDrawing() {
  std::cout << " ----------------- Testing rotation of 2D drawing."
            << std::endl;
  std::string nameBase = "test19_";
  {
    std::string smi = "Clc1ccccc1";
    std::unique_ptr<ROMol> m(SmilesToMol(smi));
    TEST_ASSERT(m);
    {
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "1.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-0' d='M 283.2 137.0") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("<text x='256.8' y='166.9' class='atom-0'") !=
                  std::string::npos);
#endif
    }
    {
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().rotate = 90.0;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "2.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("<path  class='atom-0' d='M 140.6 273.8") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("<text x='139.8' y='286.4' class='atom-0'") !=
                  std::string::npos);
#endif
    }
  }
  std::cerr << " Done" << std::endl;
}

void testGithub2063() {
  std::cout << " ----------------- Testing Github2063: Drawing racemic bond "
               "stereo as crossed bonds should be the default"
            << std::endl;
  {
    std::string molb = R"molb(squiggle bond
  Mrv1810 09301816112D          

  4  3  0  0  0  0            999 V2000
    0.5804   -0.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2948    0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1341    0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0093   -0.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  2  4  1  4  0  0  0
M  END)molb";
    std::unique_ptr<RWMol> m1(MolBlockToMol(molb));
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m1, "m1");
    drawer.addMoleculeMetadata(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub2063_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub2063_1.svg");
    TEST_ASSERT(text.find("<path class='bond-0 atom-0 atom-1' d='M "
                          "65.9,100.9 L 134.1,79.1'") != std::string::npos);
    TEST_ASSERT(text.find("<path class='bond-1 atom-0 atom-2' d='M "
                          "69.7,107.5 L 9.1,72.5'") != std::string::npos);
  }
  {
    std::string molb = R"molb(crossed bond
  Mrv1810 09301816112D          

  4  3  0  0  0  0            999 V2000
    0.5804   -0.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2948    0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1341    0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0093   -0.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  3  0  0  0
  1  3  1  0  0  0  0
  2  4  1  0  0  0  0
M  END)molb";
    std::unique_ptr<RWMol> m1(MolBlockToMol(molb));
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    MolDraw2DSVG drawer(200, 200);
    drawer.drawMolecule(*m1, "m1");
    drawer.addMoleculeMetadata(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub2063_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub2063_2.svg");
    TEST_ASSERT(text.find("<path class='bond-0 atom-0 atom-1' d='M "
                          "65.9,100.9 L 134.1,79.1'") != std::string::npos);
    TEST_ASSERT(text.find("<path class='bond-1 atom-0 atom-2' d='M "
                          "69.7,107.5 L 9.1,72.5'") != std::string::npos);
  }
  std::cerr << " Done" << std::endl;
}

void testGithub2151() {
  std::cout << " ----------------- Testing Github2151: MolDraw2D: line width "
               "should be controlled by MolDrawOptions"
            << std::endl;
  {
    auto m1 = "C[C@H](F)c1ccc(C#N)cc1"_smiles;
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    {
      MolDraw2DSVG drawer(200, 200);
      drawer.drawMolecule(*m1);
      drawer.addMoleculeMetadata(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub2151_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("testGithub2151_1.svg");
      TEST_ASSERT(text.find("stroke-width:2.0px") != std::string::npos);
      TEST_ASSERT(text.find("stroke-width:3.0px") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(200, 200);
      drawer.drawOptions().bondLineWidth = 8;
      drawer.drawMolecule(*m1);
      drawer.addMoleculeMetadata(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub2151_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("testGithub2151_2.svg");
      TEST_ASSERT(text.find("stroke-width:2.0px") == std::string::npos);
      TEST_ASSERT(text.find("stroke-width:8.0px") != std::string::npos);
    }
  }
  std::cerr << " Done" << std::endl;
}

void testGithub2762() {
  std::cout << " ----------------- Testing testGithub2762: MolDraw2D: HCl "
               "and ethane should be drawn"
            << std::endl;
  {
    auto m1 = "Cl"_smiles;
    TEST_ASSERT(m1);
    auto m2 = "CC"_smiles;
    TEST_ASSERT(m2);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    MolDraw2DUtils::prepareMolForDrawing(*m2);
    std::vector<ROMol *> mols;
    mols.push_back(m1.get());
    mols.push_back(m2.get());
    MolDraw2DSVG drawer(500, 250, 250, 250);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub2762.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub2762.svg");
    TEST_ASSERT(text.find("font-size:0px") == std::string::npos);
    TEST_ASSERT(text.find("'bond-0' d='M 0.0,200.0 L 0.0,200.0'") ==
                std::string::npos);
  }
  std::cerr << " Done" << std::endl;
}

void testGithub2931() {
  std::cout << " ----------------- Testing testGithub2931: multi-coloured"
               " molecule highlights."
            << std::endl;

  auto get_all_hit_atoms = [](ROMol &mol,
                              const std::string &smt) -> std::vector<int> {
    std::vector<int> hit_atoms;
    RWMol *query = SmartsToMol(smt);
    std::vector<MatchVectType> hits_vect;
    SubstructMatch(mol, *query, hits_vect);
    for (size_t i = 0; i < hits_vect.size(); ++i) {
      for (size_t j = 0; j < hits_vect[i].size(); ++j) {
        hit_atoms.emplace_back(hits_vect[i][j].second);
      }
    }
    delete query;
    return hit_atoms;
  };

  auto get_all_hit_bonds =
      [](ROMol &mol, const std::vector<int> &hit_atoms) -> std::vector<int> {
    std::vector<int> hit_bonds;
    for (int i : hit_atoms) {
      for (int j : hit_atoms) {
        if (i > j) {
          Bond *bnd = mol.getBondBetweenAtoms(i, j);
          if (bnd) {
            hit_bonds.emplace_back(bnd->getIdx());
          }
        }
      }
    }
    return hit_bonds;
  };

  auto update_colour_map = [](const std::vector<int> &ats, DrawColour col,
                              std::map<int, std::vector<DrawColour>> &ha_map) {
    for (auto h : ats) {
      auto ex = ha_map.find(h);
      if (ex == ha_map.end()) {
        std::vector<DrawColour> cvec(1, col);
        ha_map.insert(make_pair(h, cvec));
      } else {
        if (ex->second.end() ==
            find(ex->second.begin(), ex->second.end(), col)) {
          ex->second.emplace_back(col);
        }
      }
    }
  };

  {
    std::string smiles = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]";
    std::unique_ptr<ROMol> m(SmilesToMol(smiles));
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    std::vector<std::string> smarts = {"CONN", "N#CC~CO", "C=CON", "CONNCN"};
    std::vector<DrawColour> colours = {
        DrawColour(1.0, 0.0, 0.0), DrawColour(0.0, 1.0, 0.0),
        DrawColour(0.0, 0.0, 1.0), DrawColour(1.0, 0.55, 0.0)};
    std::map<int, std::vector<DrawColour>> ha_map;
    std::map<int, std::vector<DrawColour>> hb_map;
    for (size_t i = 0; i < smarts.size(); ++i) {
      std::vector<int> hit_atoms = get_all_hit_atoms(*m, smarts[i]);
      std::vector<int> hit_bonds = get_all_hit_bonds(*m, hit_atoms);
      update_colour_map(hit_atoms, colours[i], ha_map);
      update_colour_map(hit_bonds, colours[i], hb_map);
    }
    std::map<int, double> h_rads;
    std::map<int, int> h_lw_mult;
    {
      MolDraw2DSVG drawer(500, 500);
      drawer.drawOptions().fillHighlights = false;
      drawer.drawOptions().continuousHighlight = true;
      drawer.drawMoleculeWithHighlights(*m, "Test 1", ha_map, hb_map, h_rads,
                                        h_lw_mult);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub2931_1.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("testGithub2931_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:8.0px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='242.2' cy='367.5'"
                            " rx='10.4' ry='10.7' "
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:8.0px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='242.2' cy='313.0'"
                            " rx='10.4' ry='10.4' "
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
#endif
    }
    {
      MolDraw2DSVG drawer(500, 500);
      drawer.drawOptions().fillHighlights = false;
      drawer.drawOptions().continuousHighlight = true;
      drawer.drawOptions().atomHighlightsAreCircles = true;
      drawer.drawMoleculeWithHighlights(*m, "Test 2", ha_map, hb_map, h_rads,
                                        h_lw_mult);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub2931_2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash("testGithub2931_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:8.0px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='242.2' cy='367.0'"
                            " rx='10.5' ry='10.5' "
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
#else
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:8.0px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='242.2' cy='312.7'"
                            " rx='10.4' ry='10.4' "
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
#endif
    }
  }
  std::cerr << " Done" << std::endl;
}

void testGithub3112() {
  std::cout << " ----------------- Testing drawing of legends." << std::endl;
  {
    auto m = "CCCC"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    drawer.drawMolecule(*m, "foobar");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3112_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3112_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // this is the b
    TEST_ASSERT(text.find("<path  d='M 126.9 183.1") != std::string::npos);
#else
    TEST_ASSERT(text.find("<text x='121.0' y='195.2'"
                          " style='font-size:15px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >b</text>") != std::string::npos);
#endif
  }
  {
    auto m = "CCCC"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    drawer.drawMolecule(*m, "foo\nbar");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3112_2.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3112_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // this is the b on the 2nd line.
    TEST_ASSERT(text.find("<path  d='M 117.4 189.2") != std::string::npos);
#else
    TEST_ASSERT(text.find("<text x='110.1' y='196.2'"
                          " style='font-size:12px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >b</text>") != std::string::npos);
#endif
  }
  {
    auto m = "CCCC"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    drawer.drawMolecule(
        *m,
        "No one in their right mind would have a legend this long, surely.");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3112_3.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3112_3.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // The first letter, N.
    TEST_ASSERT(text.find("<path  d='M -1.4 187.2") != std::string::npos);
#else
    TEST_ASSERT(text.find("<text x='-2.5' y='196.8'"
                          " style='font-size:10px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >N</text>") != std::string::npos);
#endif
  }
  {
    auto m = "CCCC"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    drawer.drawMolecule(
        *m,
        "No one in their right mind would\nhave a legend this long, surely.");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3112_4.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3112_4.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // The first letter, N
    TEST_ASSERT(text.find("<path  d='M 58.3 176.7") != std::string::npos);
#else
    TEST_ASSERT(text.find("<text x='50.2' y='186.2'"
                          " style='font-size:12px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >N</text>") != std::string::npos);
#endif
  }
  std::cerr << " Done" << std::endl;
}

void test20Annotate() {
  std::cout << " ----------------- Testing annotation of 2D Drawing."
            << std::endl;

  // add serial numbers to the atoms in the molecule.
  // There's an option for this, but for debugging it's sometimes
  // useful to be able to put notes on just a few atoms.
  auto addAtomSerialNumbers = [](ROMol &mol) {
    for (auto atom : mol.atoms()) {
      atom->setProp(common_properties::atomNote, atom->getIdx());
    }
  };
  auto addBondSerialNumbers = [](ROMol &mol) {
    for (auto bond : mol.bonds()) {
      bond->setProp(common_properties::bondNote, bond->getIdx());
    }
  };
  {
    auto m1 = "S=C1N=C(NC(CC#N)(C)C=C=C)NC2=NNN=C21"_smiles;
    addAtomSerialNumbers(*m1);
    addBondSerialNumbers(*m1);
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(500, 500);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("test20_1.png");
      check_file_hash("test20_1.png");
    }
#endif

    MolDraw2DSVG drawer(500, 500);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test20_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test20_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // first note (atom 0)
    TEST_ASSERT(text.find("<path  class='note' d='M 44.9 115.7") !=
                std::string::npos);
#else
    // first one of atom note 11
    TEST_ASSERT(text.find("<text x='414.1' y='253.5' class='note'"
                          " style='font-size:12px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >1</text>") != std::string::npos);
#endif
  }

  {
    auto m1 = "C[C@@H](F)/C=C/[C@H](O)C"_smiles;
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();

      drawer.writeDrawingText("test20_2.png");
      check_file_hash("test20_2.png");
    }
#endif
    MolDraw2DSVG drawer(500, 500);
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test20_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test20_2.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // last note
    TEST_ASSERT(text.find("<path  class='note' d='M 278.1 227.5") !=
                std::string::npos);
#else
    // this is the (E)
    TEST_ASSERT(text.find("<text x='261.0' y='231.6' class='note'"
                          " style='font-size:20px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;"
                          "stroke:none;font-family:sans-serif;"
                          "text-anchor:start;fill:#000000' >E</text>") !=
                std::string::npos);
#endif
  }
  {
    auto m1 = "S=C1N=C(NC(CC#N)(C)C=C=C)NC2=NNN=C21"_smiles;
    auto atom = m1->getAtomWithIdx(3);
    atom->setProp("atomNote", "foolish annotation");
    auto bond = m1->getBondWithIdx(5);
    bond->setProp("bondNote", "way too long to be useful");
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("test20_3.png");
      check_file_hash("test20_3.png");
    }
#endif
    MolDraw2DSVG drawer(500, 500);
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test20_3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test20_3.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // first note
    TEST_ASSERT(text.find("<path  class='note' d='M 157.6 176.7") !=
                std::string::npos);
#else
    // f of foolish
    TEST_ASSERT(text.find("<text x='145.7' y='181.5' class='note'"
                          " style='font-size:12px;font-style:normal;"
                          "font-weight:normal;fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;text-anchor:start;"
                          "fill:#000000' >f</text>") != std::string::npos);
#endif
  }
  {
    auto m1 = "S=C1N=C(NC(CC#N)(C)C=C=C)NC2=NNN=C21"_smiles;
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("test20_4.png");
      check_file_hash("test20_4.png");
    }
#endif

    MolDraw2DSVG drawer(200, 200);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test20_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test20_4.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    // first note (atom 0)
    TEST_ASSERT(text.find("<path  class='note' d='M 18.0 46.3") !=
                std::string::npos);
#else
    // first one of atom note 11
    TEST_ASSERT(text.find("<text x='165.6' y='101.4' class='note' "
                          "style='font-size:4px;font-style:normal;font-weight:"
                          "normal;fill-opacity:1;stroke:none;font-family:sans-"
                          "serif;text-anchor:start;fill:#000000' >1</text>") !=
                std::string::npos);
#endif
  }
  std::cerr << " Done" << std::endl;
}

void test21FontFile() {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
  std::cout << " ----------------- Test 21 - different font" << std::endl;
  // You have to look at this one, there's no algorithmic check.
  {
    auto m = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test21_1.svg");
    MolDraw2DSVG drawer(500, 500, outs);
    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/Amadeus.ttf";
    drawer.drawOptions().fontFile = fName;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test21_1.svg");
  }
  {
    auto m = "CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    std::ofstream outs("test21_2.svg");
    MolDraw2DSVG drawer(500, 500, outs);
    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/No_Such_Font_File.ttf";
    drawer.drawOptions().fontFile = fName;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
    outs.close();
    check_file_hash("test21_2.svg");
  }
  std::cerr << "Done" << std::endl;
#endif
}

void test22ExplicitMethyl() {
  std::cout << " ----------------- Test 22 - draw explicit methyls."
            << std::endl;
  auto m = "CCC(C#C)C=C"_smiles;
  TEST_ASSERT(m);
  RDDepict::compute2DCoords(*m);
  {
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test22_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test22_1.svg");
    TEST_ASSERT(text.find("class='atom-") == std::string::npos);
  }
  {
    MolDraw2DSVG drawer(300, 300);
    drawer.drawOptions().explicitMethyl = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test22_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("test22_2.svg");
    TEST_ASSERT(text.find("class='atom-") != std::string::npos);
  }
  std::cerr << "Done" << std::endl;
}

void testGithub3305() {
  std::cout
      << " ----------------- Test Github 3305 - change and scale line widths."
      << std::endl;
  auto m = "CCC(C#C)C=C"_smiles;
  TEST_ASSERT(m);
  RDDepict::compute2DCoords(*m);
  std::string nameBase = "testGithub3305_";

  {
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(nameBase + "1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash(nameBase + "1.svg");
    TEST_ASSERT(text.find("stroke-width:2.0px") != std::string::npos);
  }
  {
    MolDraw2DSVG drawer(600, 600);
    drawer.drawOptions().bondLineWidth = 2;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(nameBase + "2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash(nameBase + "2.svg");
    TEST_ASSERT(text.find("stroke-width:2.0px") != std::string::npos);
  }
  {
    MolDraw2DSVG drawer(600, 600);
    drawer.drawOptions().bondLineWidth = 2;
    drawer.drawOptions().scaleBondWidth = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(nameBase + "3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash(nameBase + "3.svg");
    TEST_ASSERT(text.find("stroke-width:4.2px") != std::string::npos);
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  {
    MolDraw2DCairo drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + "1.png");
    check_file_hash(nameBase + "1.png");
  }
  {
    MolDraw2DCairo drawer(600, 600);
    drawer.drawOptions().bondLineWidth = 2;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + "2.png");
    check_file_hash(nameBase + "2.png");
  }
  {
    MolDraw2DCairo drawer(600, 600);
    drawer.drawOptions().bondLineWidth = 2;
    drawer.drawOptions().scaleBondWidth = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + "3.png");
    check_file_hash(nameBase + "3.png");
  }
#endif
  {
    auto m = "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"_smiles;
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));

    static const int ha[] = {17, 18, 19, 20, 21, 6, 7, 8, 9, 31, 32};
    std::vector<int> highlight_atoms(ha, ha + sizeof(ha) / sizeof(int));
    std::map<int, DrawColour> highlight_colors;
    MolDrawOptions options;
    options.circleAtoms = true;
    options.highlightColour = DrawColour(1, .5, .5);
    options.continuousHighlight = true;

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      options.scaleHighlightBondWidth = true;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + "4.png");
      check_file_hash(nameBase + "4.png");
    }
#endif
    {
      MolDraw2DSVG drawer(200, 200);
      options.scaleHighlightBondWidth = true;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "4.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "4.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:2.7") !=
                  std::string::npos);
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:16.0px") ==
                  std::string::npos);
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      options.scaleHighlightBondWidth = false;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + "5.png");
      check_file_hash(nameBase + "5.png");
    }
#endif
    {
      // This picture has very wide bond highlights as a test - it
      // looks pretty unsavoury.  I mention it so that when you flick
      // through the test images you don't panic and start searching
      // for the bug.  Been there, done that!
      MolDraw2DSVG drawer(200, 200);
      options.scaleHighlightBondWidth = false;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "5.svg").c_str());
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "5.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:2.7") ==
                  std::string::npos);
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:16.0px") !=
                  std::string::npos);
    }
    options.continuousHighlight = false;
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      options.scaleHighlightBondWidth = true;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + "6.png");
      check_file_hash(nameBase + "6.png");
    }
#endif
    {
      MolDraw2DSVG drawer(200, 200);
      options.scaleHighlightBondWidth = true;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      std::ofstream outs((nameBase + "6.svg").c_str());
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "6.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:0.7") !=
                  std::string::npos);
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4.0px") ==
                  std::string::npos);
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 200);
      options.scaleHighlightBondWidth = false;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + "7.png");
      check_file_hash(nameBase + "7.png");
    }
#endif
    {
      MolDraw2DSVG drawer(200, 200);
      options.scaleHighlightBondWidth = false;
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      std::ofstream outs((nameBase + "7.svg").c_str());
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "7.svg");
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:0.7") ==
                  std::string::npos);
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4.0px") !=
                  std::string::npos);
    }
  }
  std::cerr << "Done" << std::endl;
}

void testGithub3391() {
  std::cout
      << " ----------------- Test Github 3391 - maxFontSize interacting badly"
         " with DrawMolecules."
      << std::endl;
  auto m = "C"_smiles;
  auto m2 = "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"_smiles;
  auto m3 = "CCl"_smiles;
  {
    MolDraw2DSVG drawer(400, 200, 200, 200);
    drawer.drawOptions().maxFontSize = 14;
    std::vector<ROMol *> mols;
    mols.push_back(m.get());
    mols.push_back(m.get());
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3391_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub3391_1.svg");
  }
  {
    MolDraw2DSVG drawer(400, 200, 200, 200);
    drawer.drawOptions().maxFontSize = 14;
    std::vector<ROMol *> mols;
    mols.push_back(m.get());
    mols.push_back(m2.get());
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3391_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub3391_2.svg");
  }
  {
    MolDraw2DSVG drawer(400, 200, 200, 200);
    drawer.drawOptions().maxFontSize = 14;
    std::vector<ROMol *> mols;
    mols.push_back(m2.get());
    mols.push_back(m.get());
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3391_3.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub3391_3.svg");
  }
  {
    MolDraw2DSVG drawer(600, 200, 200, 200);
    drawer.drawOptions().maxFontSize = 14;
    drawer.drawOptions().minFontSize = 8;
    std::vector<ROMol *> mols;
    auto m1 = "CO"_smiles;
    auto m2 = "CCCCCCCCCCO"_smiles;
    auto m3 = "CCCCCCCCCCCCCCCCCCCCCCO"_smiles;
    mols.push_back(m3.get());
    mols.push_back(m2.get());
    mols.push_back(m1.get());
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3391_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub3391_4.svg");
  }

  std::cerr << "Done" << std::endl;
}

void testGithub4156() {
  std::cout
      << " ----------------- Test Github 4156 - bad scale for radicals in grid"
      << std::endl;
  auto m1 = "C1[CH]C1[C@H](F)C1CCC1"_smiles;
  auto m2 = "F[C@H]1CC[C@H](O)CC1"_smiles;
#ifdef RDK_BUILD_FREETYPE_SUPPORT
  {
    std::vector<ROMol *> mols;
    mols.push_back(m1.get());
    mols.push_back(m2.get());
    MolDraw2DSVG drawer(500, 200, 250, 200);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4156_1.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub4156_1.svg");
    // this is the start of the radical spot.
    regex qry("<path d='M 58.[0-9]*,79.[0-9]* L 58.[0-9]*,79.[0-9]*");
    TEST_ASSERT(regex_search(text, qry));
  }
  {
    std::vector<ROMol *> mols;
    mols.push_back(m2.get());
    mols.push_back(m1.get());
    MolDraw2DSVG drawer(500, 200, 250, 200);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4156_2.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub4156_2.svg");
    // this is the start of the radical spot.
    regex qry("<path d='M 308.[0-9]*,79.[0-9]* L 308.[0-9]*,79.[0-9]*");
    TEST_ASSERT(regex_search(text, qry));
  }
#endif
  std::cerr << " Done" << std::endl;
}

void test23JSONAtomColourPalette() {
  std::cerr << " ----------------- Test JSON atomColourPalette" << std::endl;
  {
    auto m = "c1cncs1"_smiles;
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    const char *json =
        R"JSON({"atomColourPalette": {"7": [0.2, 0.4, 0.9], "16": [0.9, 0.6, 0.0]},
 "rotate": 90, "bondLineWidth": 5})JSON";
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, json);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test23_1.svg");
    outs << text;
    outs.close();
    check_file_hash("test23_1.svg");
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    TEST_ASSERT(text.find("' fill='#3366E5") != std::string::npos);
    TEST_ASSERT(text.find("' fill='#E59900") != std::string::npos);
#else
    TEST_ASSERT(text.find("fill:#3366E5") != std::string::npos);
    TEST_ASSERT(text.find("fill:#E59900") != std::string::npos);
#endif
  }
  std::cerr << " Done" << std::endl;
}

void testGithub4496() {
  std::cerr << " ----------------- Test draw aromatic OR queries" << std::endl;
  auto m = "[c,n]1[c,n][c,n][c,n][c,n][c,n]1"_smarts;
  MolDraw2DSVG drawer(200, 200);
  MolDrawOptions options;
  options.prepareMolsBeforeDrawing = false;
  drawer.drawOptions() = options;
  MolDraw2DUtils::prepareMolForDrawing(*m, false);
  drawer.drawMolecule(*m);
  drawer.finishDrawing();
  std::ofstream outs("testGithub4496_1.svg");
  std::string text = drawer.getDrawingText();
  outs << text;
  outs.flush();
  outs.close();
  check_file_hash("testGithub4496_1.svg");
  std::cerr << " Done" << std::endl;
}

int main() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = false;
#endif

  RDLog::InitLogs();

#if 1
  test1();
  test2();
  test4();
  test5();
  test6();
  test7();
  test8PrepareMolForDrawing();
  testMultiThreaded();
  testGithub781();
  test3();
  testGithub774();
  test9MolLegends();
  testGithub852();
  testGithub860();
  testGithub910();
  testGithub932();
  testGithub953();
  testGithub983();
  testDeuteriumTritium();
  testCrossedBonds();
  test10DrawSecondMol();
  test11DrawMolGrid();
  test12DrawMols();
  test13JSONConfig();
  testGithub1090();
  testGithub1035();
  testGithub1271();
  testGithub1322();
  test14BWPalette();
  test15ContinuousHighlightingWithGrid();
  test17MaxMinFontSize();
  testGithub1829();
  test18FixedScales();
  test19RotateDrawing();
  test16MoleculeMetadata();
  testGithub2063();
  testGithub2151();
  testGithub2762();
  testGithub2931();
  test20Annotate();
  test21FontFile();
  test22ExplicitMethyl();
  testGithub3112();
  testGithub3305();
  testGithub3391();
  testGithub4156();
  test23JSONAtomColourPalette();
  testGithub4496();
#endif
}
