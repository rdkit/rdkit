//
//  Copyright (C) 2015 Greg Landrum 
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//

#include "MolDraw2DSVG.h"
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <sstream>

namespace RDKit {
  namespace {
    std::string DrawColourToSVG(const MolDraw2D::DrawColour &col){
      const char *convert="0123456789ABCDEF";
      std::string res(7,' ');
      res[0]='#';
      unsigned int v;
      unsigned int i=1;
      v=int(255 * col.get<0>());
      res[i++] = convert[v/16];
      res[i++] = convert[v%16];
      v=int(255 * col.get<1>());
      res[i++] = convert[v/16];
      res[i++] = convert[v%16];
      v=int(255 * col.get<2>());
      res[i++] = convert[v/16];
      res[i++] = convert[v%16];
      return res;
    }
  }
  
  void MolDraw2DSVG::initDrawing() {
    d_os<<"<?xml version='1.0' encoding='iso-8859-1'?>\n";
    d_os << "<svg:svg version='1.1' baseProfile='full'\n      \
        xmlns:svg='http://www.w3.org/2000/svg'\n                \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
    d_os<<"width='"<<width()<<"px' height='"<<height()<<"px' >\n";
    //d_os<<"<svg:g transform='translate("<<width()*.05<<","<<height()*.05<<") scale(.85,.85)'>";
  }

  // ****************************************************************************
  void MolDraw2DSVG::finishDrawing() {
    //d_os << "</svg:g>";
    d_os << "</svg:svg>\n";
  }

  // ****************************************************************************
  void MolDraw2DSVG::setColour( const DrawColour &col ) {
    MolDraw2D::setColour( col );
  }

  // ****************************************************************************
  void MolDraw2DSVG::drawLine( const std::pair<float,float> &cds1 ,
                               const std::pair<float,float> &cds2 ) {

    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );
    std::string col=DrawColourToSVG(colour());
    unsigned int width=2;
    std::string dashString="";
    const DashPattern &dashes=dash();
    if(dashes.size()){
      std::stringstream dss;
      dss<<";stroke-dasharray:";
      std::copy(dashes.begin(),dashes.end()-1,std::ostream_iterator<unsigned int>(dss,","));
      dss<<dashes.back();
      dashString = dss.str();
    }
    d_os<<"<svg:path ";
    d_os<< "d='M " << c1.first << "," << c1.second << " " << c2.first << "," << c2.second << "' ";
    d_os<<"style='fill:none;fill-rule:evenodd;stroke:"<<col<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
    d_os<<" />\n";
  }

  // ****************************************************************************
  // draw the char, with the bottom left hand corner at cds
  void MolDraw2DSVG::drawChar( char c , const std::pair<float,float> &cds ) {
    unsigned int fontSz=scale()*fontSize();
    std::string col = DrawColourToSVG(colour());

    d_os<<"<svg:text";
    d_os<<" x='" << cds.first;
    d_os<< "' y='" << cds.second + fontSz <<"'"; // doesn't seem like this should be necessary, but vertical text alignment seems impossible
    d_os<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;"<<"fill:"<<col<<"'";
    d_os<<" >";
    d_os<<c;
    d_os<<"</svg:text>";
  }

  // ****************************************************************************
  void MolDraw2DSVG::drawTriangle( const std::pair<float , float> &cds1 ,
                                   const std::pair<float , float> &cds2 ,
                                   const std::pair<float, float> &cds3 ) {
    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );
    std::pair<float,float> c3 = getDrawCoords( cds3 );

    std::string col=DrawColourToSVG(colour());
    unsigned int width=2;
    std::string dashString="";
    d_os<<"<svg:path ";
    d_os<< "d='M " << c1.first << "," << c1.second << " " << c2.first << "," << c2.second << " " << c3.first << "," << c3.second << "' ";
    d_os<<"style='fill:"<<col<<";fill-rule:evenodd;stroke:"<<col<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
    d_os<<" />\n";

  }

  // ****************************************************************************
  void MolDraw2DSVG::clearDrawing() {
    d_os << "<svg:rect";
    d_os << " style='opacity:1.0;fill:#ffffff;stroke:none'";
    d_os << " width='" << width() << "' height='" << height() << "'";
    d_os << " x='0' y='0'";
    d_os << "> </svg:rect>\n";
  }


  
  // ****************************************************************************
  void MolDraw2DSVG::setFontSize( float new_size ) {
    MolDraw2D::setFontSize( new_size );
    float font_size_in_points = fontSize() * scale();
  }

  // ****************************************************************************
  // using the current scale, work out the size of the label in molecule coordinates
  void MolDraw2DSVG::getStringSize( const std::string &label , float &label_width ,
                                    float &label_height ) const {

    label_width = 0.0;
    label_height = 0.0;

    int draw_mode = 0; // 0 for normal, 1 for superscript, 2 for subscript

    bool had_a_super = false;

    for( int i = 0 , is = label.length() ; i < is ; ++i ) {

      // setStringDrawMode moves i along to the end of any <sub> or <sup>
      // markup
      if( '<' == label[i] && setStringDrawMode( label , draw_mode , i ) ) {
        continue;
      }

      label_height = fontSize();
      float char_width = fontSize() * static_cast<float>(MolDraw2D_detail::char_widths[label[i]]) / MolDraw2D_detail::char_widths['M'];
      //char_width *= 0.75; // extremely empirical
      if( 2 == draw_mode ) {
        char_width *= 0.75;
      } else if( 1 == draw_mode ) {
        char_width *= 0.75;
        had_a_super = true;
      }
      label_width += char_width;

    }

    // subscript keeps its bottom in line with the bottom of the bit chars,
    // superscript goes above the original char top by a quarter
    if( had_a_super ) {
      label_height *= 1.25;
    }
  }

  // ****************************************************************************
  // draws the string centred on cds
  void MolDraw2DSVG::drawString( const std::string &str ,
                                 const std::pair<float,float> &cds ) {

    unsigned int fontSz=scale()*fontSize();
    std::string col = DrawColourToSVG(colour());

    float string_width , string_height;
    getStringSize( str , string_width , string_height );

    float draw_x = cds.first - string_width / 2.0;
    float draw_y = cds.second - string_height / 2.0;
    std::pair<float,float> draw_coords = getDrawCoords(std::make_pair(draw_x,draw_y));

    d_os<<"<svg:text";
    d_os<<" x='" << draw_coords.first;
    d_os<< "' y='" << draw_coords.second + fontSz <<"'"; // doesn't seem like this should be necessary, but vertical text alignment seems impossible
    d_os<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;"<<"fill:"<<col<<"'";
    d_os<<" >";

    int draw_mode = 0; // 0 for normal, 1 for superscript, 2 for subscript
    std::string span;
    bool first_span=true;
    for( int i = 0 , is = str.length() ; i < is ; ++i ) {
      // setStringDrawMode moves i along to the end of any <sub> or <sup>
      // markup
      if( '<' == str[i] && setStringDrawMode( str , draw_mode , i ) ) {
        if(!first_span){
          d_os<<span<<"</svg:tspan>";
          span="";
        }
        first_span=false;
        d_os<<"<svg:tspan";
        switch(draw_mode){
        case 1:
          d_os<<" style='baseline-shift:super;font-size:"<<fontSz*0.75<<"px;"<< "'"; break;
        case 2:
          d_os<<" style='baseline-shift:sub;font-size:"<<fontSz*0.75<<"px;"<< "'"; break;
        default:
          break;
        }
        d_os<<">";
        continue;
      }
      if(first_span){
        first_span=false;
        d_os<<"<svg:tspan>";
        span="";
      }
      span += str[i];
    }
    d_os<<span<<"</svg:tspan>";
    d_os<<"</svg:text>\n";
  }
} // EO namespace RDKit
