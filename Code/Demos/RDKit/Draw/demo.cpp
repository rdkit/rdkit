// $Id$
//
//  Copyright (C) 2009-2012 Greg Landrum
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>

#include <RDGeneral/RDLog.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cairo.h>

using namespace RDKit;

std::string MolToSVG(const ROMol &mol){
  RWMol cp(mol);
  RDKit::MolOps::Kekulize(cp);
  RDDepict::compute2DCoords(cp);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(cp);
  std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
  return svg;
}


void setColorCairo(int atNum,cairo_t *cr){
  PRECONDITION(cr,"no context");
  switch(atNum){
  case 7:
    cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);
    break;
  case 8:
    cairo_set_source_rgb (cr, 1.0, 0.0, 0.0);
    break;
  case 9:
    cairo_set_source_rgb (cr, 0.2, 0.8, 0.8);
    break;
  case 15:
    cairo_set_source_rgb (cr, 1.0, 0.5, 0.0);
    break;
  case 16:
    cairo_set_source_rgb (cr, 0.8, 0.8, 0.0);
    break;
  case 17:
    cairo_set_source_rgb (cr, 0.0, 0.8, 0.0);
    break;
  case 35:
    cairo_set_source_rgb (cr, 0.5, 0.3, 0.1);
    break;
  case 0:
    cairo_set_source_rgb (cr, 0.5, 0.5, 0.5);
    break;
  default:
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
  }
}
void drawLineCairo(std::vector<int>::const_iterator &pos,
                 cairo_t *cr){
  PRECONDITION(cr,"no context");
  int width=*pos++;
  cairo_set_line_width(cr,width*10);
  int dashed=*pos++;
  double dashes[2];
  switch(dashed){
  case 0:
    cairo_set_dash(cr,0,0,0);
    break;
  case 2:
    dashes[0]=5.0;dashes[1]=20.0;
    cairo_set_dash(cr,dashes,sizeof(dashes)/sizeof(dashes[0]),0);    
    break;
  default:
    dashes[0]=20.0;dashes[1]=20.0;
    cairo_set_dash(cr,dashes,sizeof(dashes)/sizeof(dashes[0]),0);    
  }
  int an1=*pos++;
  int an2=*pos++;
  int xp1 = *pos++;
  int yp1 = *pos++;
  int xp2 = *pos++;
  int yp2 = *pos++;
  if(an1==an2){
    setColorCairo(an1,cr);
    cairo_move_to(cr,xp1,yp1);
    cairo_line_to(cr,xp2,yp2);
    cairo_stroke(cr);
  } else {
    int mx = xp1+(xp2-xp1)/2;
    int my = yp1+(yp2-yp1)/2;

    setColorCairo(an1,cr);
    cairo_move_to(cr,xp1,yp1);
    cairo_line_to(cr,mx,my);
    cairo_stroke(cr);
    setColorCairo(an2,cr);
    cairo_move_to(cr,mx,my);
    cairo_line_to(cr,xp2,yp2);
    cairo_stroke(cr);
  }
}
void drawAtomCairo(std::vector<int>::const_iterator &pos,
                 cairo_t *cr){
  PRECONDITION(cr,"no context");
  int fontSz=50;
  int atNum=*pos++;
  double xp=static_cast<double>(*pos++);
  double yp=static_cast<double>(*pos++);
  int slen=*pos++;
  std::string label="";
  for(unsigned int i=0;i<slen;++i){
    label+=(char)*pos++;
  }
  RDKit::Drawing::OrientType orient=static_cast<RDKit::Drawing::OrientType>(*pos++);
  cairo_text_extents_t extents;
  cairo_text_extents(cr,label.c_str(),&extents);
  double twidth=extents.width,theight=extents.height;

  switch(orient){
  case RDKit::Drawing::W:
    xp -= twidth;
    yp += theight/2;
    break;
  case RDKit::Drawing::E:
    yp += theight/2;
    break;
  case RDKit::Drawing::S:
    xp -= twidth/2;
    yp += theight;
    break;
  case RDKit::Drawing::N:
    xp -= twidth/2;
    yp -= theight/2;
    break;
  default:
    xp -= twidth/2;
    yp += theight/2;
  }
  cairo_set_source_rgb (cr, 1.0, 1., 1.);
  cairo_rectangle(cr,
                  xp-10,yp-theight-10,
                  twidth+20,theight+20);
  cairo_fill(cr);
                  

  cairo_move_to(cr,xp,yp);
  setColorCairo(atNum,cr);
  cairo_show_text(cr,label.c_str());
  cairo_stroke(cr);
}


void MolToCairo(const ROMol &mol,cairo_t *cr,int width,int height,
                int fontSize=14,int maxDotsPerAngstrom=30){
  PRECONDITION(cr,"no context");
  PRECONDITION(width>0 && height>0,"bad dimensions");
  RWMol *cp = new RWMol(mol);
  try{
    RDKit::MolOps::Kekulize(*cp);
  } catch (...) {
    delete cp;
    cp = new RDKit::RWMol(mol);
  }
  RDDepict::compute2DCoords(*cp);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(*cp);
  delete cp;

  cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr,0,0,width,height);
  cairo_fill(cr);

  cairo_select_font_face (cr, "sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);

  std::vector<int>::const_iterator pos=drawing.begin();
  int resolution=0;
  if(*pos!=Drawing::RESOLUTION){
    std::cerr<<"no RESOLUTION token found"<<std::endl;
    return;
  }
  resolution=*(pos+1);
  PRECONDITION(resolution>0,"bad resolution");
  pos+=2;
  if(*pos!=Drawing::BOUNDS){
    std::cerr<<"no bounds token found"<<std::endl;
    return;
  }
  pos+=3;
  int dwidth,dheight;
  dwidth = *pos++;
  dheight = *pos++;
  PRECONDITION(dwidth>0 && dheight>0,"bad dimensions");

  // size of the image in angstroms:
  double xSize,ySize;
  xSize=static_cast<double>(dwidth)/resolution;
  ySize=static_cast<double>(dheight)/resolution;
  double scale=1.0;

  if(dwidth>width || dheight>height){
    if(static_cast<double>(dwidth)/width > static_cast<double>(dheight)/height){
      scale = static_cast<double>(width)/dwidth;
    } else {
      scale = static_cast<double>(height)/dheight;
    }
  } else {
    if(width/xSize > height/ySize){
      if( width/xSize > maxDotsPerAngstrom ){
        scale = maxDotsPerAngstrom*xSize/width;
      }
    } else {
      if( height/ySize > maxDotsPerAngstrom ){
        scale = maxDotsPerAngstrom*ySize/height;
      }
    }
  }
  
  scale*=0.80;
  cairo_translate(cr,
                  .5*(width-dwidth*scale),
                  .5*(height-dheight*scale));
  cairo_scale(cr,scale,scale);

  // scaling factors here are empirically determined
  double dFontSize=2.5*maxDotsPerAngstrom*fontSize/14;
  cairo_set_font_size (cr, dFontSize);


  while(pos!=drawing.end()){
    int token=*pos++;
    switch(token){
    case Drawing::LINE:
      drawLineCairo(pos,cr);
      break;
    case Drawing::ATOM:
      drawAtomCairo(pos,cr);
      break;
    default:
      std::cerr<<"unrecognized token: "<<token<<std::endl;
    }
  }
}

void DrawDemo(){

#if 1
  //RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1~C1CC1");
  RWMol *mol=SmilesToMol("c1c(C#N)cccc1C~C2CC2");
  std::string svg=MolToSVG(*mol);
  std::ofstream ostr("blah.svg");
  ostr<<svg<<std::endl;
  delete mol;
#else
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
#endif
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  DrawDemo();
}
