// $Id$
//
//  Copyright (C) 2009-2012 Greg Landrum
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define USE_CAIRO 1

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>
#ifdef USE_CAIRO
#include <GraphMol/MolDrawing/DrawingToCairo.h>
#endif

#include <RDGeneral/RDLog.h>
#include <vector>
#include <sstream>
#include <fstream>


using namespace RDKit;

std::string MolToSVG(const ROMol &mol){
  std::vector<int> drawing=RDKit::Drawing::MolToDrawing(mol);
  std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
  return svg;
}


#ifdef USE_CAIRO
void MolToCairo(const ROMol &mol,cairo_t *cr,int width,int height){
  PRECONDITION(cr,"no context");
  PRECONDITION(width>0 && height>0,"bad dimensions");
  std::vector<int> drawing=RDKit::Drawing::MolToDrawing(mol);
  RDKit::Drawing::DrawingToCairo(drawing,cr,width,height);
}
#endif

void DrawDemo(){

#if 1
  {
    RWMol *mol=SmilesToMol("c1c(C#N)cccc1C~C2CC2");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah3.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCCCCCC(CCCCCCCCCCC(CCCCCCCCC)(CCCCCCCCCC)CCCCCC)CCCCCCCCCCCCCCCCCCCCCc2sc([NH3+])c([NH3+])c2)cc1");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah4.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("BrO");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah2.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("BrC(O)(Cl)N");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah5.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("[14NH2+]=[14NH2+]");
    std::string svg=MolToSVG(*mol);
    std::ofstream ostr("blah6.svg");
    ostr<<svg<<std::endl;
    delete mol;
  }


#endif
#ifdef USE_CAIRO
  {
    RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 200, 200);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,200,200);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol1.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("c1c[12cH]ccn1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol2.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("Nc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol3.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("Nccc(CCO)n",0,false);
    mol->updatePropertyCache();
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol4.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("N~ccc(CCO)n",0,false);
    mol->updatePropertyCache();
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol5.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("BrO",0,false);
    mol->updatePropertyCache();
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol6.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCCCCCC(CCCCCCCCCCC(CCCCCCCCC)(CCCCCCCCCC)CCCCCC)CCCCCCCCCCCCCCCCCCCCCc2sc([NH3+])c([NH3+])c2)cc1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol7.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("[NH3+][NH3+]",0,false);
    mol->updatePropertyCache();
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol8.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
#endif
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  DrawDemo();
}
