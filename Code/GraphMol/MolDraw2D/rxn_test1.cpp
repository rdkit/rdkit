//
//  Copyright (C) 2017 Greg Landrum
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
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>

using namespace RDKit;

#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <cairo.h>
#include "MolDraw2DCairo.h"
#endif

namespace {

// if the generated SVG hashes to the value we're expecting, delete
// the file.  That way, only the files that need inspection will be
// left at the end of the run.
static const bool DELETE_WITH_GOOD_HASH = true;
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    static const std::map<std::string, std::hash_result_t> SVG_HASHES = {
    {"rxn_test1_1.svg", 4142488780U},
    {"rxn_test1_2.svg", 2391235449U},
    {"rxn_test1_3.svg", 1596835824U},
    {"rxn_test1_4.svg", 2658830567U},
    {"rxn_test1_5.svg", 3734154337U},
    {"rxn_test1_6.svg", 1382455147U},
    {"rxn_test1_7.svg", 850006879U},
    {"rxn_test2_1.svg", 3558612822U},
    {"rxn_test2_2_1.svg", 2748718550U},
    {"rxn_test2_2_2.svg", 3634147033U},
    {"rxn_test2_2_3.svg", 1683554313U},
    {"rxn_test2_2_4.svg", 3008925605U},
    {"rxn_test3_1.svg", 409355256U},
    {"rxn_test4_1.svg", 1446031830U},
    {"rxn_test4_2.svg", 3254692708U}
};
#else
static const std::map<std::string, std::hash_result_t> SVG_HASHES = {
    {"rxn_test1_1.svg", 3455349472U},
    {"rxn_test1_2.svg", 1549953602U},
    {"rxn_test1_3.svg", 3887971655U},
    {"rxn_test1_4.svg", 4154131138U},
    {"rxn_test1_5.svg", 817859479U},
    {"rxn_test1_6.svg", 3148235567U},
    {"rxn_test1_7.svg", 1494204470U},
    {"rxn_test2_1.svg", 3794156425U},
    {"rxn_test2_2_1.svg", 905466032U},
    {"rxn_test2_2_2.svg", 937555927U},
    {"rxn_test2_2_3.svg", 571938327U},
    {"rxn_test2_2_4.svg", 2056889618U},
    {"rxn_test3_1.svg", 2623086145U},
    {"rxn_test4_1.svg", 3737522161U},
    {"rxn_test4_2.svg", 873707471U}
};
#endif

// These PNG hashes aren't completely reliable due to floating point cruft,
// but they can still reduce the number of drawings that need visual
// inspection.  At present, the files
// rxn_test1_2.png  rxn_test1_5.png  rxn_test2_2_1.png  rxn_test2_2_4.png
// rxn_test4_2.png rxn_test1_3.png  rxn_test1_7.png  rxn_test2_2_2.png
// rxn_test3_1.png rxn_test1_4.png  rxn_test2_1.png  rxn_test2_2_3.png
// rxn_test4_1.png
// give different results on my MBP and Ubuntu 20.04 VM.  The SVGs work
// better because the floats are all output to only 1 decimal place so there
// is a much smaller chance of different systems producing different files.
static const std::map<std::string, std::hash_result_t> PNG_HASHES = {
    {"rxn_test1_1.png", 1360104321U},
    {"rxn_test1_2.png", 2193860546U},
    {"rxn_test1_3.png", 3637813871U},
    {"rxn_test1_4.png", 2689827672U},
    {"rxn_test1_5.png", 1081681982U},
    {"rxn_test1_6.png", 1652373318U},
    {"rxn_test1_7.png", 1129854251U},
    {"rxn_test2_1.png", 2184440504U},
    {"rxn_test2_2_1.png", 2843973557U},
    {"rxn_test2_2_2.png", 3681916546U},
    {"rxn_test2_2_3.png", 2917064875U},
    {"rxn_test2_2_4.png", 1045873872U},
    {"rxn_test3_1.png", 1460167503U},
    {"rxn_test4_1.png", 2326178421U},
    {"rxn_test4_2.png", 380646467U}
};

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
                     std::hash_result_t exp_hash=0U) {
//    std::cout << filename << " : " << hash_file(filename) << "U" << std::endl;

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

void drawit(ChemicalReaction *rxn, std::string nameBase,
            bool highlight_map = false,
            const std::vector<DrawColour> *highlight_colors = nullptr) {
  double panex = 200, paney = 150;
  double width = panex * (rxn->getNumReactantTemplates() +
                          rxn->getNumProductTemplates() + 1);
  double height = paney;
#ifdef RDK_BUILD_CAIRO_SUPPORT
  {
    MolDraw2DCairo drawer(width, height);
    drawer.drawReaction(*rxn, highlight_map, highlight_colors);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + ".png");
    check_file_hash(nameBase + ".png");
  }
#endif
  {
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(width, height, outs);
    drawer.drawReaction(*rxn, highlight_map, highlight_colors);
    drawer.finishDrawing();
    outs.close();
    check_file_hash(nameBase + ".svg");
  }
}
}  // namespace

void test1() {
  std::cout << " ----------------- Test 1" << std::endl;
  {
    std::string smiles =
        "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=["
        "O:3])[NH:6][CH3:5].[OH2:4]";
    std::string nameBase = "rxn_test1_1";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }

  {
    std::string smiles =
        "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:"
        "2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
    std::string nameBase = "rxn_test1_2";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }

  {
    std::string smiles =
        ">>[N:1]1[C:"
        "2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
    std::string nameBase = "rxn_test1_3";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }
  {
    std::string smiles =
        "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>";
    std::string nameBase = "rxn_test1_4";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }

  {
    std::string smiles =
        "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>O.ClCl>";
    std::string nameBase = "rxn_test1_5";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }

  {
    std::string smiles =
        "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=["
        "O:3])[NH:6][CH3:5].[OH2:4]";
    std::string nameBase = "rxn_test1_6";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase, true);
    delete rxn;
  }

  {
    std::string smiles =
        "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:"
        "2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
    std::string nameBase = "rxn_test1_7";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase, true);
    delete rxn;
  }
  std::cout << " Done" << std::endl;
}

void test2() {
  std::cout << " ----------------- Test 2: some examples that caused problems "
               "in testing"
            << std::endl;
  {  // from the reaction role assignment paper
    std::string smiles =
        "[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2]"
        "(=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1["
        "cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])"
        "[N:13]=[N+:14]=[N-:15]";
    std::string nameBase = "rxn_test2_1";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase);
    delete rxn;
  }
  {  // from Ed Griffen
    //     order is: SMIRKS	substrate	product
    std::string indata =
        "[cH:1]1[cH:2][cH:3][c:4][cH:5][cH:6]1>>CC(C)[c:1]1[cH:2][cH:3][c:4][cH:5][cH:6]1	CN1CCC(CC1)N(C)S(=O)(=O)c2ccccc2	CC(C)c1ccc(cc1)S(=O)(=O)N(C)C2CCN(CC2)C\n\
[cH:1]1[cH:2][cH:3][c:4][cH:5][cH:6]1>>[cH:3]1[cH:2][c:1]([cH:6][cH:5][c:4]1)Cl	Cc1cccc(c1C)N2CCN(CC2)CCCNC(=O)c3cc([nH]c3C)c4ccccc4	Cc1cccc(c1C)N2CCN(CC2)CCCNC(=O)c3cc([nH]c3C)c4ccc(cc4)Cl\n\
Cc1c(cc([nH]1)[c:1]2[cH:2][cH:3][cH:4][cH:5][cH:6]2)[C:7](=[O:8])[NH:9][CH2:10][CH2:11]>>CCc1nc(c(n1[c:1]2[cH:2][cH:3][cH:4][cH:5][cH:6]2)C)[C:7](=[O:8])[NH:9][CH2:10][CH2:11]	Cc1c(cc([nH]1)c2ccccc2)C(=O)NCCN3CCN(C\n\
C3)c4cccc(c4Cl)Cl	CCc1nc(c(n1c2ccccc2)C)C(=O)NCCN3CCN(CC3)c4cccc(c4Cl)Cl\n\
[c:1][C:2](=[O:3])[NH:4]CCC[CH2:5][N:6]([CH2:7])[CH2:8]>>[c:1][C:2](=[O:3])[NH:4]C[CH2:5][N:6]([CH2:7])[CH2:8]	CCCn1c(c(cc1c2ccccc2)C(=O)NCCCCN3CCN(CC3)c4cccc(c4Cl)Cl)C	CCCn1c(c(cc1c2ccccc2)C(=O)NCCN3CCN(CC3)c4cccc(c4Cl)Cl)C\n\
[CH3:1]c1c(cc(n1[CH3:2])C(C)(C)C)[C:3]>>[CH3:2]c1c(c(nc(n1)[CH3:1])[C:3])OC	Cc1cccc(c1C)N2CCN(CC2)CCNC(=O)c3cc(n(c3C)C)C(C)(C)C	Cc1cccc(c1C)N2CCN(CC2)CCNC(=O)c3c(c(nc(n3)C)C)OC";
    std::vector<std::string> lines;
    boost::split(lines, indata, boost::is_any_of("\n"));
    unsigned int idx = 0;
    for (auto &line : lines) {
      std::vector<std::string> tokens;
      boost::split(tokens, line, boost::is_any_of("\t "));
      std::cerr << tokens.size() << " " << line << std::endl;
      if (tokens.size() > 2) {
        ++idx;

        std::string nameBase = "rxn_test2_2_" + std::to_string(idx);
        bool useSmiles = true;
        ChemicalReaction *rxn =
            RxnSmartsToChemicalReaction(tokens[0], nullptr, useSmiles);
        TEST_ASSERT(rxn);
        drawit(rxn, nameBase);
        delete rxn;
      }
    }
  }
  std::cout << " Done" << std::endl;
}

void test3() {
  std::cout << " ----------------- Test 3: test reactant highlighting"
            << std::endl;
  {  // from the reaction role assignment paper
    std::string smiles =
        "[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2]"
        "(=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1["
        "cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])"
        "[N:13]=[N+:14]=[N-:15]";
    std::string nameBase = "rxn_test3_1";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    std::vector<DrawColour> highlight_colors;
    highlight_colors.push_back(DrawColour(1., 1., .67));
    highlight_colors.push_back(DrawColour(1., .71, .76));
    highlight_colors.push_back(DrawColour(.8, 1., .8));
    highlight_colors.push_back(DrawColour(.67, .67, 1.));
    drawit(rxn, nameBase, true, &highlight_colors);
    delete rxn;
  }
  std::cout << " Done" << std::endl;
}

void test4() {
  std::cout << " ----------------- Test 4: some examples from SMARTS "
               "in testing"
            << std::endl;
  {  // from the reaction role assignment paper
    std::string smiles =
        "[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2]"
        "(=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1["
        "cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])"
        "[N:13]=[N+:14]=[N-:15]";
    std::string nameBase = "rxn_test4_1";
    bool useSmiles = false;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    std::cerr << "draw!" << std::endl;
    drawit(rxn, nameBase, true);
    std::cerr << " done" << std::endl;
    delete rxn;
  }
  {  // from the reaction role assignment paper
    std::string smiles =
        "[CH3:16][CH2:15][c:14]1[c:7]([cH:6][c:5]([cH:22][c:17]1[O:18][CH2:19]["
        "CH:20]=[CH2:21])[O:4][CH2:3][CH:2]=[CH2:1])[CH2:8][CH2:9][O:10][CH2:"
        "11][CH2:12][O:13]C2CCCCO2.CO.C1COCCO1.C(=O)(O)[O-].O.[Na+].Cl>>[CH3:"
        "16][CH2:15][c:14]1[c:7]([cH:6][c:5]([cH:22][c:17]1[O:18][CH2:19][CH:"
        "20]=[CH2:21])[O:4][CH2:3][CH:2]=[CH2:1])[CH2:8][CH2:9][O:10][CH2:11]["
        "CH2:12][OH:13]";
    std::string nameBase = "rxn_test4_2";
    bool useSmiles = false;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, nullptr, useSmiles);
    TEST_ASSERT(rxn);
    drawit(rxn, nameBase, true);
    delete rxn;
  }
  std::cout << " Done" << std::endl;
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
#endif
}
