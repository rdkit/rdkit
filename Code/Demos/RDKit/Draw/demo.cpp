// $Id: sample.cpp 793 2008-08-17 14:33:30Z glandrum $
//
//  Copyright (C) 2009 Greg Landrum
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <Demos/RDKit/Draw/MolDrawing.h>

#include <RDGeneral/RDLog.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cairo.h>

using namespace RDKit;

std::string getColorSVG(int atNum){
  static std::map<int,std::string> colors;
  if(colors.empty()){
    colors[7]="#0000FF";
    colors[8]="#FF0000";
    colors[9]="#33CCCC";
    colors[15]="#FF7F00";
    colors[16]="#CCCC00";
    colors[17]="#00CC00";
    colors[35]="#7F4C19";
    colors[0]="#7F7F7F";
  }
  std::string res="#000000";
  if(colors.find(atNum)!=colors.end()) res= colors[atNum];
  return res;
}
void drawLineSVG(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
  int width=*pos++;
  int dashed=*pos++;
  int an1=*pos++;
  int an2=*pos++;
  std::string c1=getColorSVG(an1);
  std::string c2=getColorSVG(an2);
  if(c1==c2){
    sstr<<"<svg:path ";
    sstr<<"d='M "<<*pos<<","<<*(pos+1)<<" "<<*(pos+2)<<","<<*(pos+3)<<"' ";
    pos+=4;
    sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
    sstr<<" />\n";
  } else {
    int xp1 = *pos++;
    int yp1 = *pos++;
    int xp2 = *pos++;
    int yp2 = *pos++;
    int mx = xp1+(xp2-xp1)/2;
    int my = yp1+(yp2-yp1)/2;
    sstr<<"<svg:path ";
    sstr<<"d='M "<<xp1<<","<<yp1<<" "<<mx<<","<<my<<"' ";
    sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
    sstr<<" />\n";
    sstr<<"<svg:path ";
    sstr<<"d='M "<<mx<<","<<my<<" "<<xp2<<","<<yp2<<"' ";
    sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c2<<";stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
    sstr<<" />\n";
  }
}
void drawAtomSVG(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
  int fontSz=50;
  int atNum=*pos++;
  int xp=*pos++;
  int yp=*pos++;
  int slen=*pos++;
  std::string label="";
  for(unsigned int i=0;i<slen;++i){
    label+=(char)*pos++;
  }
  int orient=*pos++;
  int width=fontSz*label.length();
  int height=fontSz;
  sstr<<"<svg:g transform='translate("<<xp<<","<<yp<<")'><svg:rect ";
  sstr<<"style='opacity:1.0;fill:#FFFFFF;stroke:none'";
  sstr<<" width='"<<width<<"' height='"<<height<<"'";
  sstr<<" x='-"<<width/2<<"' y='-"<<height/2<<"'";
  sstr<<"> </svg:rect>\n";
  sstr<<"<svg:text";
  sstr<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:Sans;text-anchor:middle"<<";fill:"<<getColorSVG(atNum)<<"'";
  sstr<<" y='"<<.75*fontSz/2<<"'>";
  sstr<<"<svg:tspan>";
  sstr<<label<<"</svg:tspan>";
  sstr<<"</svg:text>";
  sstr<<"</svg:g>\n";
}

std::string ToSVG(const std::vector<int> &drawing){
  std::vector<int>::const_iterator pos=drawing.begin()+2;
  if(*pos!=Drawing::BOUNDS){
    std::cerr<<"no bounds token found"<<std::endl;
    return "";
  }
  pos+=3;
  int width=300,height=300;
  width = *pos++;
  height = *pos++;
  std::ostringstream sstr;
  sstr<<"<?xml version='1.0' encoding='iso-8859-1'?>\n";
  sstr << "<svg:svg version='1.1' baseProfile='full'\n      \
        xmlns:svg='http://www.w3.org/2000/svg'\n                \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
  sstr<<"width='"<<width<<"px' height='"<<height<<"px' >\n";
  sstr<<"<svg:g transform='translate("<<width*.05<<","<<height*.05<<") scale(.9,.9)'>";
  while(pos!=drawing.end()){
    int token=*pos++;
    switch(token){
    case Drawing::LINE:
      drawLineSVG(pos,sstr);
      break;
    case Drawing::ATOM:
      drawAtomSVG(pos,sstr);
      break;
    default:
      std::cerr<<"unrecognized token: "<<token<<std::endl;
    }
  }

  sstr<<"</svg:g></svg:svg>";
  return sstr.str();
}

std::string MolToSVG(const ROMol &mol){
  RWMol cp(mol);
  RDKit::MolOps::Kekulize(cp);
  RDDepict::compute2DCoords(cp);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(cp);
  std::string svg=ToSVG(drawing);
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
  std::cerr<<" label: "<<label<<" "<<orient<<std::endl;
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
  RWMol cp(mol);
  RDKit::MolOps::Kekulize(cp);
  RDDepict::compute2DCoords(cp);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(cp);

  cairo_select_font_face (cr, "sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

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
  
  std::cerr<<" sz: "<<width<<" "<<height<<" : "<<dwidth<<" "<<dheight<<std::endl;
  std::cerr<<" scale: "<<scale<<std::endl;
  scale*=0.80;
  cairo_translate(cr,
                  .5*(width-dwidth*scale),
                  .5*(height-dheight*scale));
  cairo_scale(cr,scale,scale);

  cairo_set_font_size (cr, static_cast<double>(fontSize)/scale);


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

#if 0
  RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
  std::string svg=MolToSVG(*mol);
  std::ofstream ostr("blah.svg");
  ostr<<svg<<std::endl;
  delete mol;
#else
  {
    RWMol *mol=SmilesToMol("[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol1.png");
    cairo_surface_destroy (surface);
    delete mol;
  }
  {
    RWMol *mol=SmilesToMol("c1ccccn1");
    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);
    MolToCairo(*mol,cr,300,300);

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "mol2.png");
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
