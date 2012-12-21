// $Id$
/* 
*
*  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _RD_DRAWING_TO_SVG_H_
#define _RD_DRAWING_TO_SVG_H_

#include <GraphMol/MolDrawing/MolDrawing.h>

namespace RDKit {
  namespace Drawing {
    namespace {
      std::string getColor(int atNum){
        std::string res="#000000";
        switch(atNum){
        case 7: res="#0000FF";break;
        case 8: res="#FF0000";break;
        case 9: res="#33CCCC";break;
        case 15: res="#FF7F00";break;
        case 16: res="#CCCC00";break;
        case 17: res="#00CC00";break;
        case 35: res="#7F4C19";break;
        case 0: res="#7F7F7F";break;
        default: res="#000000";
        }
        return res;
      }
      void drawLine(std::vector<int>::const_iterator &pos,std::ostringstream &sstr,
                    unsigned int lineWidthMult){
        int width=*pos++;
        width*=lineWidthMult;

        int dashed=*pos++;
        std::string dashString="";
        switch(dashed){
        case 0:
          break;
        case 2:
          dashString=";stroke-dasharray:2, 6";
          break;
        default:
          dashString=";stroke-dasharray:6, 6";
        }
        int an1=*pos++;
        int an2=*pos++;
        std::string c1=getColor(an1);
        std::string c2=getColor(an2);
        if(c1==c2){
          sstr<<"<svg:path ";
          sstr<<"d='M "<<*pos<<","<<*(pos+1)<<" "<<*(pos+2)<<","<<*(pos+3)<<"' ";
          pos+=4;
          sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
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
          sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c1<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
          sstr<<" />\n";
          sstr<<"<svg:path ";
          sstr<<"d='M "<<mx<<","<<my<<" "<<xp2<<","<<yp2<<"' ";
          sstr<<"style='fill:none;fill-rule:evenodd;stroke:"<<c2<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
          sstr<<" />\n";
        }
      }

      void drawAtom(std::vector<int>::const_iterator &pos,std::ostringstream &sstr,unsigned int fontSz){
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
        sstr<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:Sans;text-anchor:middle"<<";fill:"<<getColor(atNum)<<"'";
        sstr<<" y='"<<.75*fontSz/2<<"'>";
        sstr<<"<svg:tspan>";
        sstr<<label<<"</svg:tspan>";
        sstr<<"</svg:text>";
        sstr<<"</svg:g>\n";
      }
    } // end of anonymous namespace 

    std::string DrawingToSVG(const std::vector<int> &drawing,
                             unsigned int lineWidthMult=2,unsigned int fontSize=50){
      std::vector<int>::const_iterator pos=drawing.begin()+2;
      if(*pos!= RDKit::Drawing::BOUNDS){
        std::cerr<<"no bounds token found"<<std::endl;
        return "";
      }
      pos+=3;
      unsigned int width,height;
      width = *pos++;
      height = *pos++;
      std::ostringstream sstr;
      sstr<<"<?xml version='1.0' encoding='iso-8859-1'?>\n";
      sstr << "<svg:svg version='1.1' baseProfile='full'\n      \
        xmlns:svg='http://www.w3.org/2000/svg'\n                \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
      sstr<<"width='"<<width<<"px' height='"<<height<<"px' >\n";
      sstr<<"<svg:rect ";
      sstr<<"style='opacity:1.0;fill:#FFFFFF;stroke:none'";
      sstr<<" width='"<<width<<"' height='"<<height<<"'";
      sstr<<" x='0' y='0'";
      sstr<<"> </svg:rect>\n";
      sstr<<"<svg:g transform='translate("<<width*.05<<","<<height*.05<<") scale(.85,.85)'>";
      while(pos!=drawing.end()){
        int token=*pos++;
        switch(token){
        case  RDKit::Drawing::LINE:
          drawLine(pos,sstr,lineWidthMult);
          break;
        case  RDKit::Drawing::ATOM:
          drawAtom(pos,sstr,fontSize);
          break;
        default:
          std::cerr<<"unrecognized token: "<<token<<std::endl;
        }
      }
      sstr<<"</svg:g></svg:svg>";
      return sstr.str();
    }
  } // end of namespace Drawing
} // end of namespace RDKit
#endif
