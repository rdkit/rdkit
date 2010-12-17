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


using namespace RDKit;

std::string getColor(int atNum){
  static std::map<int,std::string> colors;
  if(colors.empty()){
    colors[7]="#0000FF";
    colors[8]="#FF0000";
    colors[9]="#33CCCC";
    colors[15]="#FF7F00";
    colors[16]="#CCCC00";
    colors[16]="#00CC00";
    colors[35]="#7F4C19";
    colors[0]="#7F7F7F";
  }
  std::string res="#000000";
  if(colors.find(atNum)!=colors.end()) res= colors[atNum];
  return res;
}
void drawLine(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
      pos+=4;
      sstr<<"<svg:path ";
      sstr<<"d='M "<<*pos<<","<<*(pos+1)<<" "<<*(pos+2)<<","<<*(pos+3)<<"' ";
      pos+=4;
      sstr<<"style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:3px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1'";
      sstr<<" />\n";
}
void drawAtom(std::vector<int>::const_iterator &pos,std::ostringstream &sstr){
  int atNum=*pos++;
  int xp=*pos++;
  int yp=*pos++;
  int slen=*pos++;
  std::string label="";
  for(unsigned int i=0;i<slen;++i){
    label+=(char)*pos++;
  }
  sstr<<"<svg:text x='"<<xp<<"' y='"<<yp<<"'";
  sstr<<" style='font-size:40px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:Sans;text-anchor:middle;baseline-shift:sub"<<";fill:"<<getColor(atNum)<<"'";
  sstr<<">";
  sstr<<"<svg:tspan x='"<<xp<<"' y='"<<yp<<"'>";
  sstr<<label<<"</svg:tspan>";
  sstr<<"</svg:text>\n";
}

std::string ToSVG(const std::vector<int> &drawing){
  std::ostringstream sstr;
  sstr<<"<?xml version='1.0' encoding='iso-8859-1'?>\n";
  int width=300,height=300;
  std::vector<int>::const_iterator pos=drawing.begin()+2;
  TEST_ASSERT(*pos==Drawing::BOUNDS);
  pos+=3;
  width = *pos++;
  height = *pos++;
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
      drawLine(pos,sstr);
      break;
    case Drawing::ATOM:
      drawAtom(pos,sstr);
      break;
    default:
      std::cerr<<"unrecognized token: "<<token<<std::endl;
    }
  }

  sstr<<"</svg:g></svg:svg>";
  return sstr.str();
}

void DrawDemo(){
  RWMol *mol=SmilesToMol("Clc1c(C#N)cc(C(=O)NCc2sccc2)cc1");
  //RWMol *mol=SmilesToMol("c1ncncn1");
  RDKit::MolOps::Kekulize(*mol);

  // generate the 2D coordinates:
  RDDepict::compute2DCoords(*mol);
  std::vector<int> drawing=RDKit::Drawing::DrawMol(*mol);

  std::cerr<<"var codes=[";
  std::copy(drawing.begin(),drawing.end(),std::ostream_iterator<int>(std::cerr,","));
  std::cerr<<"];"<<std::endl;

  std::string svg=ToSVG(drawing);
  std::cout<<svg<<std::endl;

  delete mol;
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  DrawDemo();
}
