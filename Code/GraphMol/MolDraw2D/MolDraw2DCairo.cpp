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

#include "MolDraw2DCairo.h"
#include <cairo.h>

namespace RDKit {
  void MolDraw2DCairo::initDrawing() {
    cairo_select_font_face(d_cr, "sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  }

  // ****************************************************************************
  void MolDraw2DCairo::finishDrawing() {
  }

  // ****************************************************************************
  void MolDraw2DCairo::setColour( const DrawColour &col ) {
    MolDraw2D::setColour( col );
    cairo_set_source_rgb(d_cr, col.get<0>(), col.get<1>(), col.get<2>());
  }

  // ****************************************************************************
  void MolDraw2DCairo::drawLine( const std::pair<float,float> &cds1 ,
                               const std::pair<float,float> &cds2 ) {

    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );

    unsigned int width=2;
    std::string dashString="";

    cairo_set_line_width(d_cr,width);

    const DashPattern &dashes=dash();
    if(dashes.size()){
      double dd[dashes.size()];
      for(unsigned int di=0;di<dashes.size();++di) dd[di]=dashes[di]*1.5;
      cairo_set_dash(d_cr,dd,dashes.size(),0);
    } else {
      cairo_set_dash(d_cr,0,0,0);
    }

    cairo_move_to(d_cr,c1.first,c1.second);
    cairo_line_to(d_cr,c2.first,c2.second);
    cairo_stroke(d_cr);
  }

  // ****************************************************************************
  // draw the char, with the bottom left hand corner at cds
  void MolDraw2DCairo::drawChar( char c , const std::pair<float,float> &cds ) {
    char txt[2];
    txt[0]=c;
    txt[1]=0;

    cairo_text_extents_t extents;
    cairo_text_extents(d_cr,txt,&extents);
    double twidth=extents.width,theight=extents.height;

    unsigned int fontSz=scale()*fontSize();
    std::pair<float,float> c1 = cds ;// getDrawCoords( cds );

    cairo_move_to(d_cr,c1.first,c1.second+theight);
    cairo_show_text(d_cr,txt);
    cairo_stroke(d_cr);
  }

  // ****************************************************************************
  void MolDraw2DCairo::drawTriangle( const std::pair<float , float> &cds1 ,
                                   const std::pair<float , float> &cds2 ,
                                   const std::pair<float, float> &cds3 ) {
    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );
    std::pair<float,float> c3 = getDrawCoords( cds3 );

    unsigned int width=1;
    cairo_set_line_width(d_cr,width);
    cairo_set_dash(d_cr,0,0,0);
    cairo_move_to(d_cr,c1.first,c1.second);
    cairo_line_to(d_cr,c2.first,c2.second);
    cairo_line_to(d_cr,c3.first,c3.second);
    cairo_close_path(d_cr);
    cairo_fill_preserve(d_cr);
    cairo_stroke(d_cr);
  }

  // ****************************************************************************
  void MolDraw2DCairo::clearDrawing() {
    
    cairo_set_source_rgb (d_cr, 1.0, 1.0, 1.0);
    cairo_rectangle(d_cr,0,0,width(),height());
    cairo_fill(d_cr);
  }


  
  // ****************************************************************************
  void MolDraw2DCairo::setFontSize( float new_size ) {
    MolDraw2D::setFontSize( new_size );
    float font_size_in_points = fontSize() * scale();
    cairo_set_font_size (d_cr, font_size_in_points);
  }

  // ****************************************************************************
  // using the current scale, work out the size of the label in molecule coordinates
  void MolDraw2DCairo::getStringSize( const std::string &label , float &label_width ,
                                    float &label_height ) const {
    label_width = 0.0;
    label_height = 0.0;

    int draw_mode = 0; // 0 for normal, 1 for superscript, 2 for subscript

    bool had_a_super = false;

    char txt[2];
    txt[1]=0;
    for( int i = 0 , is = label.length() ; i < is ; ++i ) {

      // setStringDrawMode moves i along to the end of any <sub> or <sup>
      // markup
      if( '<' == label[i] && setStringDrawMode( label , draw_mode , i ) ) {
        continue;
      }

      txt[0]=label[i];
      cairo_text_extents_t extents;
      cairo_text_extents(d_cr,txt,&extents);
      double twidth=extents.x_advance,theight=extents.height;

      label_height = theight/scale();
      float char_width = twidth/scale();

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
    label_height *= 1.2; // empirical
  }
} // EO namespace RDKit
