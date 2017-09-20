//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
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

using namespace RDKit;

#ifdef RDK_CAIRO_BUILD
#include <cairo.h>
#include "MolDraw2DCairo.h"
#endif

namespace {
void drawit(ChemicalReaction *rxn, std::string nameBase,
            bool highlight_map = false,
            const std::vector<DrawColour> *highlight_colors = NULL) {
  double panex = 200, paney = 150;
  double width = panex * (rxn->getNumReactantTemplates() +
                          rxn->getNumProductTemplates() + 1);
  double height = paney;
#ifdef RDK_CAIRO_BUILD
  {
    MolDraw2DCairo drawer(width, height);
    drawer.drawReaction(*rxn, highlight_map, highlight_colors);
    drawer.finishDrawing();
    drawer.writeDrawingText(nameBase + ".png");
  }
#endif
  {
    std::ofstream outs((nameBase + ".svg").c_str());
    MolDraw2DSVG drawer(width, height, outs);
    drawer.drawReaction(*rxn, highlight_map, highlight_colors);
    drawer.finishDrawing();
    outs.flush();
  }
}
}

void test1() {
  std::cout << " ----------------- Test 1" << std::endl;
  {
    std::string smiles =
        "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][NH2:6]>CC(O)C.[Pt]>[CH3:1][C:2](=["
        "O:3])[NH:6][CH3:5].[OH2:4]";
    std::string nameBase = "rxn_test1_1";
    bool useSmiles = true;
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
    BOOST_FOREACH (std::string &line, lines) {
      std::vector<std::string> tokens;
      boost::split(tokens, line, boost::is_any_of("\t "));
      std::cerr << tokens.size() << " " << line << std::endl;
      if (tokens.size() > 2) {
        ++idx;

        std::string nameBase =
            "rxn_test2_2_" + boost::lexical_cast<std::string>(idx);
        bool useSmiles = true;
        ChemicalReaction *rxn =
            RxnSmartsToChemicalReaction(tokens[0], NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
        RxnSmartsToChemicalReaction(smiles, NULL, useSmiles);
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
