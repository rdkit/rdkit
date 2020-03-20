//
//  Copyright (C) 2015-2017 Greg Landrum
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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace RDKit;

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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().atomLabels = atomLabels;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().circleAtoms = false;
      drawer.drawMolecule(*m, &highlight_atoms);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions().circleAtoms = true;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(200, 200, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m, &highlight_atoms, &highlight_colors);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawOptions() = options;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    std::string nameBase = "test5_1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs("test6_1.svg");
    outs << txt;
    // TEST_ASSERT(txt.find("<svg")!=std::string::npos);
    delete m;
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
    std::string smiles = "C";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<tspan>CH</tspan>") != std::string::npos);
    delete m;
  }
  {
    std::string smiles = "O";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<tspan>OH</tspan>") == std::string::npos);
    delete m;
  }
  {
    std::string smiles = "[C]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<tspan>C</tspan>") != std::string::npos);
    delete m;
  }
  {
    std::string smiles = "C.CC.[Cl-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    TEST_ASSERT(txt.find("<svg") != std::string::npos);
    TEST_ASSERT(txt.find("<tspan>CH</tspan>") != std::string::npos);
    TEST_ASSERT(txt.find("<tspan>Cl</tspan>") != std::string::npos);
    delete m;
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
    TEST_ASSERT(txt.find("<tspan>") == std::string::npos);
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    std::string smiles = "CC[13CH2][CH2:7][CH-]C[15NH2+]C";
    std::string nameBase = "test5_1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m, &(m->getConformer()));
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m, "mol legend");
    drawer.finishDrawing();
    std::string txt = drawer.getDrawingText();
    std::ofstream outs("test9_1.svg");
    outs << txt;
    // TEST_ASSERT(txt.find("<svg")!=std::string::npos);
    delete m;
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
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
    std::string smiles = "[15NH3+:1]-C#C-[15NH3+:2]";
    std::string nameBase = "test860_1";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
    }
    delete m;
  }
  {
    std::string smiles =
        "[15NH3+:1]-C#C-C([15NH3+:2])([15NH3+:3])-C#C-[15NH3+:4]";
    std::string nameBase = "test860_2";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
    }
    delete m;
  }
  {
    std::string smiles = "[15NH3+:1]-CCCCCCCC-[15NH3+:4]";
    std::string nameBase = "test860_3";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);

#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText(nameBase + ".png");
    }
#endif
    {
      std::ofstream outs((nameBase + ".svg").c_str());
      MolDraw2DSVG drawer(300, 300, outs);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      outs.flush();
    }
    delete m;
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
    TEST_ASSERT(
        text.find("<path class='bond-1' d='M 125.352,114.634"
                  " L 179.234,90.896 L 172.848,79.8357 Z'"
                  " style='fill:#000000") != std::string::npos);
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
    TEST_ASSERT(
        text.find("<path class='bond-3' d='M 102.962,116.968 L"
                  " 77.0085,93.6784 L 72.5976,99.8217 Z'"
                  " style='fill:#000000;") != std::string::npos);

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
    std::string deuterium = "C([2H])([2H])([2H])[2H]";
    ROMol *m = SmilesToMol(deuterium);
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
      if ((line.find("baseline-shift:super") != std::string::npos) &&
          (line.find(">2<") != std::string::npos) &&
          (line.find(">H<") != std::string::npos)) {
        ++count;
      }
    }
    TEST_ASSERT(count == 4);
    delete m;
  }
  {
    std::string tritium = "C([3H])([3H])([3H])[3H]";
    ROMol *m = SmilesToMol(tritium);
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
      if ((line.find("baseline-shift:super") != std::string::npos) &&
          (line.find(">3<") != std::string::npos) &&
          (line.find(">H<") != std::string::npos)) {
        ++count;
      }
    }
    TEST_ASSERT(count == 4);
    delete m;
  }
  {
    std::string deuterium = "C([2H])([2H])([2H])[2H]";
    ROMol *m = SmilesToMol(deuterium);
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
      if ((line.find("baseline-shift:super") == std::string::npos) &&
          (line.find(">2<") == std::string::npos) &&
          (line.find(">D<") != std::string::npos)) {
        ++count;
      }
    }
    TEST_ASSERT(count == 4);
    delete m;
  }
  {
    std::string tritium = "C([3H])([3H])([3H])[3H]";
    ROMol *m = SmilesToMol(tritium);
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
      if ((line.find("baseline-shift:super") == std::string::npos) &&
          (line.find(">3<") == std::string::npos) &&
          (line.find(">T<") != std::string::npos)) {
        ++count;
      }
    }
    TEST_ASSERT(count == 4);
    delete m;
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
  }
  delete m1;
  delete m2;
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
  std::unique_ptr<std::vector<std::string>> legends(new std::vector<std::string>());
  // made up SMILES, each with sequence F, Cl, Br so we can see which
  // ones are drawn, which ones are missing.
  setup_mol("COc1cccc(NC(=O)[C@H](F)Sc2nc(ns2)c3ccccc3F)c1", "m1",
            mols, *legends);
  setup_mol("NC(=O)[C@H](F)Sc1ncns1", "m2", mols, *legends);
  setup_mol("COc1cccc(NC(=O)[C@H](Cl)Sc2nc(ns2)c3ccccc3F)c1", "m3",
            mols, *legends);
  setup_mol("NC(=O)[C@H](Cl)Sc1ncns1", "m4", mols, *legends);
  setup_mol("COc1cccc(NC(=O)[C@H](Br)Sc2nc(ns2)c3ccccc3F)c1", "m5",
            mols, *legends);
  setup_mol("NC(=O)[C@H](Br)Sc1ncns1", "m6", mols, *legends);

  {
    MolDraw2DSVG drawer(750, 400, 250, 200);
    drawer.drawMolecules(mols, legends.get());
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_1.svg");
    outs << text;
    outs.flush();
  }

  {  // github #1325: multiple molecules in one pane
    MolDraw2DSVG drawer(300, 300, 300, 300);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_3.svg");
    outs << text;
    outs.flush();
  }

  {  // github #1325: multiple molecules in one pane
    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_4.svg");
    outs << text;
    outs.flush();
  }
  {
    mols[2] = nullptr;
    mols[4] = nullptr;
    MolDraw2DSVG drawer(750, 400, 250, 200);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test12_2.svg");
    outs << text;
    outs.flush();
  }
  for(auto m: mols) {
    delete m;
  }
  std::cerr << " Done" << std::endl;
}

void test13JSONConfig() {
  std::cerr << " ----------------- Test JSON Configuration" << std::endl;
  {
    std::string smiles = "CCO";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 200);
    const char *json = "{\"legendColour\":[1.0,0.5,1.0], \"rotate\": 90, "
                       "\"bondLineWidth\": 5}";
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, json);
    drawer.drawMolecule(*m, "foo");
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test13_1.svg");
    outs << text;
    outs.close();
    TEST_ASSERT(text.find("sans-serif;fill:#FF7FFF") !=
                std::string::npos);
    TEST_ASSERT(text.find("'bond-0' d='M 122.883,9.09091 L 170.762,92.0201'")
                != std::string::npos);
    // these days the bond line width scales with the rest of the
    // drawing, and at this size this comes out as 6px.
    TEST_ASSERT(text.find("stroke-width:6px") !=
                std::string::npos);
    delete m;
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
    TEST_ASSERT(text.find("C&1") == std::string::npos);
    TEST_ASSERT(text.find("<1") == std::string::npos);
    TEST_ASSERT(text.find(">1") == std::string::npos);
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
    delete m1;
    delete m2;
  }

  std::cerr << " Done" << std::endl;
}

void testGithub1322() {
  std::cout << " ----------------- Testing github 1322: add custom atom labels"
            << std::endl;
  {
    std::string smiles = "CCC[Se]";  // made up
    RWMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      TEST_ASSERT(text.find("Se") != std::string::npos);
    }
    {
      m1->getAtomWithIdx(3)->setProp(common_properties::atomLabel,
                                     "customlabel");
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test1322_1.svg");
      outs << text;
      outs.flush();
      TEST_ASSERT(text.find("Se") == std::string::npos);
      TEST_ASSERT(text.find("customlabel") != std::string::npos);
    }
    delete m1;
  }
  std::cerr << " Done" << std::endl;
}

void testGithub565() {
  std::cout << " ----------------- Testing github 565: support a fixed bond "
               "length in the MolDraw2D code"
            << std::endl;
  {
    std::string smiles = "CCCCC";
    RWMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    Point2D minV, maxV;
    const Conformer &cnf = m1->getConformer();
    minV.x = maxV.x = cnf.getAtomPos(0).x;
    minV.y = maxV.y = cnf.getAtomPos(0).y;
    for (unsigned int i = 1; i < m1->getNumAtoms(); i++) {
      minV.x = std::min(minV.x, cnf.getAtomPos(i).x);
      minV.y = std::min(minV.y, cnf.getAtomPos(i).y);
      maxV.x = std::max(maxV.x, cnf.getAtomPos(i).x);
      maxV.y = std::max(maxV.y, cnf.getAtomPos(i).y);
    }

    {
      unsigned int dpa = 100;
      unsigned int w = dpa * (maxV.x - minV.x);
      unsigned int h = dpa * (maxV.y - minV.y);

      MolDraw2DSVG drawer(w, h);
      drawer.setScale(w, h, minV, maxV);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test565_1.svg");
      outs << text;
      outs.flush();
    }
    {
      unsigned int dpa = 50;
      unsigned int w = dpa * (maxV.x - minV.x);
      unsigned int h = dpa * (maxV.y - minV.y);

      MolDraw2DSVG drawer(w, h);
      drawer.setScale(w, h, minV, maxV);
      drawer.drawMolecule(*m1, "m1");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("test565_2.svg");
      outs << text;
      outs.flush();
    }
    delete m1;
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
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4px;") ==
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
      TEST_ASSERT(text.find("stroke:#FF7F7F;stroke-width:4px;") !=
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
      TEST_ASSERT(text.find("idx=\"2\" atom-smiles=\"[NH]\" drawing-x=\"54.") !=
                  std::string::npos);
      TEST_ASSERT(text.find("idx=\"2\" begin-atom-idx=\"3\" end-atom-idx=\"2\" "
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

      TEST_ASSERT(text.find("atom-smiles=\"[NH]\" drawing-x=\"54.") !=
                  std::string::npos);
      TEST_ASSERT(text.find("atom-smiles=\"[NH]\" drawing-x=\"254.") !=
                  std::string::npos);

      for (auto ptr : ms) {
        delete ptr;
      }
    }
  }

  std::cerr << " Done" << std::endl;
}

void test17MaxFontSize() {
  std::cout << " ----------------- Test 16 - Testing maximum font size"
            << std::endl;
  {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/MolDraw2D/test_dir";
    fName += "/clash.mol";
    std::unique_ptr<ROMol> m(MolFileToMol(fName));
    std::string nameBase = "test16_";
    TEST_ASSERT(m);

    {
      std::ofstream outs((nameBase + "1.svg").c_str());
      MolDraw2DSVG drawer(300, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      outs << text;
      outs.flush();
      TEST_ASSERT(text.find("font-size:40px") != std::string::npos);
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
      TEST_ASSERT(text.find("font-size:47px") != std::string::npos);
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
      TEST_ASSERT(text.find("font-size:20px") != std::string::npos);
    }
  }

  std::cerr << " Done" << std::endl;
}

void test18FixedScales() {
  std::cout << " ----------------- Testing use of fixed scales for drawing."
            << std::endl;
  std::string nameBase = "test18_";
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
      TEST_ASSERT(text.find("font-size:28px") != std::string::npos);
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
      TEST_ASSERT(text.find("font-size:7px") != std::string::npos);
    }
  }
  {
    std::string smi = "C[C@@H](N[C@@H]1CC[C@@H](C(=O)N2CCC(C(=O)N3CCCC3)"
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
      TEST_ASSERT(text.find("font-size:8px") != std::string::npos);
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
      TEST_ASSERT(text.find("font-size:5px") != std::string::npos);
    }
    {
      // this one should be the same size as the first, as it won't scale
      // up if the picture won't fit.
      MolDraw2DSVG drawer(300, 300);
      drawer.drawOptions().fixedBondLength = 30;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs((nameBase + "5.svg").c_str());
      outs << text;
      outs.flush();
      TEST_ASSERT(text.find("font-size:8px") != std::string::npos);
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
      TEST_ASSERT(text.find("text-anchor=\"start\" x='256.907' y='154.276'")
                  != std::string::npos);
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
      TEST_ASSERT(text.find("text-anchor=\"start\" x='136.604' y='276.316'")
                  != std::string::npos);
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
    TEST_ASSERT(
        text.find(
            "<path class='bond-0' d='M 65.8823,110.884 L 134.118,89.1159'") !=
        std::string::npos);
    TEST_ASSERT(
        text.find(
            "<path class='bond-1' d='M 69.6998,117.496 L 9.09091,82.5044'") !=
        std::string::npos);
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
    TEST_ASSERT(
        text.find(
            "<path class='bond-0' d='M 65.8823,110.884 L 134.118,89.1159'") !=
        std::string::npos);
    TEST_ASSERT(
        text.find(
            "<path class='bond-1' d='M 69.6998,117.496 L 9.09091,82.5044'") !=
        std::string::npos);
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
      TEST_ASSERT(text.find("stroke-width:2px") != std::string::npos);
      TEST_ASSERT(text.find("stroke-width:3px") == std::string::npos);
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
      // the bonds are scaled in thickness, so it won't be 8 pixels.
      // Experiment finds 3 on my machine.
      TEST_ASSERT(text.find("stroke-width:2px") == std::string::npos);
      TEST_ASSERT(text.find("stroke-width:3px") != std::string::npos);
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
    TEST_ASSERT(text.find("font-size:0px") == std::string::npos);
    TEST_ASSERT(text.find("'bond-0' d='M 0,200 L 0,200'") == std::string::npos);
  }
  std::cerr << " Done" << std::endl;
}

void testGithub2931() {
  std::cout << " ----------------- Testing testGithub2931: multi-coloured"
               " molecule highlights."
            << std::endl;

  auto get_all_hit_atoms = [](ROMol &mol, const std::string &smt) -> std::vector<int> {
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

  auto get_all_hit_bonds = [](ROMol &mol, const std::vector<int> &hit_atoms) -> std::vector<int> {
    std::vector<int> hit_bonds;
    for (int i: hit_atoms) {
      for (int j: hit_atoms) {
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
                              std::map<int, std::vector<DrawColour> > &ha_map) {
    for (auto h: ats) {
      auto ex = ha_map.find(h);
      if (ex == ha_map.end()) {
        std::vector<DrawColour> cvec(1, col);
        ha_map.insert(make_pair(h, cvec));
      } else {
        if (ex->second.end() == find(ex->second.begin(), ex->second.end(), col)) {
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
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:5px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='241.942' cy='386.438'"
                            " rx='11.9978' ry='12.846'"
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
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
      TEST_ASSERT(text.find("stroke:#FF8C00;stroke-width:5px") !=
                  std::string::npos);
      TEST_ASSERT(text.find("<ellipse cx='241.766' cy='385.788'"
                            " rx='10.9782' ry='10.9782'"
                            " style='fill:none;stroke:#00FF00;") !=
                  std::string::npos);
    }
  }
  std::cerr << " Done" << std::endl;
}

void test20Annotate() {
    std::cout << " ----------------- Testing annotation of 2D Drawing."
            << std::endl;

  // add serial numbers to the atoms in the molecule
  auto addAtomSerialNumbers = [](ROMol &mol) {
      for(auto atom: mol.atoms()) {
        atom->setProp("atomNote", atom->getIdx());
      }
  };
  auto addBondSerialNumbers = [](ROMol &mol) {
    for(auto bond: mol.bonds()) {
      bond->setProp("bondNote", bond->getIdx());
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
    }
#endif

    MolDraw2DSVG drawer(500, 500);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("test20_1.svg");
    outs << text;
    outs.flush();
    TEST_ASSERT(text.find("x='402.344' y='249.673' style='font-size:14px;"
                          "font-style:normal;font-weight:normal;"
                          "fill-opacity:1;stroke:none;font-family:sans-serif;"
                          "fill:#000000' ><tspan>11") != std::string::npos);
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
    TEST_ASSERT(text.find("x='278.629' y='221.042' style='font-size:27px;"
                          "font-style:normal;font-weight:normal;"
                          "fill-opacity:1;stroke:none;font-family:sans-serif;"
                          "fill:#000000") != std::string::npos);
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
    TEST_ASSERT(text.find("x='209.32' y='180.482' style='font-size:15px;"
                          "font-style:normal;font-weight:normal;"
                          "fill-opacity:1;stroke:none;"
                          "font-family:sans-serif;fill:#000000'"
                          " ><tspan>foolish annotation</tspan>")
                != std::string::npos);
  }
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
  testGithub565();
  test14BWPalette();
  test15ContinuousHighlightingWithGrid();
  test17MaxFontSize();
  testGithub1829();
  test18FixedScales();
  test19RotateDrawing();
#endif
  test16MoleculeMetadata();
  testGithub2063();
  testGithub2151();
  testGithub2762();
  testGithub2931();
  test20Annotate();
}
