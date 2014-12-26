//
// file MolDraw2DSVG.cc
// Greg Landrum
// 26 Dec 2014
//
// derived from Dave Cosgrove's MolDraw2DQT
//

#include "MolDraw2DSVG.H"

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
    // QColor this_col( int( 255.0 * col.get<0>() ) , int( 255.0 * col.get<1>() ) ,
    //                  int( 255.0 * col.get<2>() ) );

    // QPen pen( this_col );
    // pen.setJoinStyle( Qt::RoundJoin );
    // pen.setColor( this_col );
    // qp_.setPen( pen );

    // QBrush brush( this_col );
    // brush.setStyle( Qt::SolidPattern );
    // qp_.setBrush( brush );

  }

  // ****************************************************************************
  void MolDraw2DSVG::drawLine( const std::pair<float,float> &cds1 ,
                               const std::pair<float,float> &cds2 ) {

    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );
    std::string col=DrawColourToSVG(getColour());
    unsigned int width=2;
    std::string dashString="";
    d_os<<"<svg:path ";
    d_os<< "d='M " << c1.first << "," << c1.second << " " << c2.first << "," << c2.second << "' ";
    d_os<<"style='fill:none;fill-rule:evenodd;stroke:"<<col<<";stroke-width:"<<width<<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"<<dashString<<"'";
    d_os<<" />\n";
  }

  // ****************************************************************************
  // draw the char, with the bottom left hand corner at cds
  void MolDraw2DSVG::drawChar( char c , const std::pair<float,float> &cds ) {
    unsigned int fontSz=scale()*fontSize();
    std::string col = DrawColourToSVG(getColour());

    d_os<<"<svg:text";
    d_os<<" x='" << cds.first;
    d_os<< "' y='" << cds.second + fontSz <<"'"; // doesn't seem like this should be necessary, but vertical text alignment seems impossible
    d_os<<" style='font-size:"<<fontSz<<"px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;"<<"fill:"<<col<<"'";
    d_os<<" >";
    d_os<<c;
    d_os<<"</svg:text>";

    // QRectF br = qp_.boundingRect( 0 , 0 , 100 , 100 , Qt::AlignLeft | Qt::AlignBottom ,
    //                               QString( c ) );
    // qp_.drawText( QRectF( cds.first , cds.second , br.width() , br.height() ) ,
    //               Qt::AlignLeft | Qt::AlignBottom , QString( c ) , &br );

  }

  // ****************************************************************************
  void MolDraw2DSVG::drawTriangle( const std::pair<float , float> &cds1 ,
                                   const std::pair<float , float> &cds2 ,
                                   const std::pair<float, float> &cds3 ) {
    std::pair<float,float> c1 = getDrawCoords( cds1 );
    std::pair<float,float> c2 = getDrawCoords( cds2 );
    std::pair<float,float> c3 = getDrawCoords( cds3 );

    std::string col=DrawColourToSVG(getColour());
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
    d_os << " style='opacity:1.0;fill:#eeeeee;stroke:none'";
    d_os << " width='" << width() << "' height='" << height() << "'";
    d_os << " x='0' y='0'";
    d_os << "> </svg:rect>\n";
  }


  
  // ****************************************************************************
  void MolDraw2DSVG::setFontSize( float new_size ) {
    std::cerr<<"sfs: "<<new_size<<" "<<scale()<<std::endl;
    MolDraw2D::setFontSize( new_size );
    float font_size_in_points = fontSize() * scale();

    // QFont font( qp_.font() );
    // font.setPointSizeF( font_size_in_points );
    // qp_.setFont( font );

    // while( 1 ) {
    //   float old_font_size_in_points = font_size_in_points;
    //   float font_size_in_points = fontSize() * scale();
    //   if( fabs( font_size_in_points - old_font_size_in_points ) < 0.1 ) {
    //     break;
    //   }
    //   QFont font( qp_.font() );
    //   font.setPointSizeF( font_size_in_points );
    //   qp_.setFont( font );
    //   calculateScale();
    // }

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

      label_height = fontSize()*scale() / scale();
      float char_width = fontSize()*scale() / scale();
      char_width *= 0.75; // extremely empirical
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

} // EO namespace RDKit
