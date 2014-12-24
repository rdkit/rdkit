//
// file MolDraw2DQt.cc
// David Cosgrove
// AstraZeneca
// 19th June 2014
//

#include "MolDraw2DQt.H"

#include <QPainter>
#include <QString>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDraw2DQt::MolDraw2DQt( int width , int height , QPainter &qp ) :
MolDraw2D( width , height ) , qp_( qp ) {

}

// ****************************************************************************
void MolDraw2DQt::setColour( const DrawColour &col ) {

  MolDraw2D::setColour( col );
  QColor this_col( int( 255.0 * col.get<0>() ) , int( 255.0 * col.get<1>() ) ,
                   int( 255.0 * col.get<2>() ) );

  QPen pen( this_col );
  pen.setJoinStyle( Qt::RoundJoin );
  pen.setColor( this_col );
  qp_.setPen( pen );

  QBrush brush( this_col );
  brush.setStyle( Qt::SolidPattern );
  qp_.setBrush( brush );

}

// ****************************************************************************
void MolDraw2DQt::drawLine( const pair<float,float> &cds1 ,
                            const pair<float,float> &cds2 ) {

  pair<float,float> c1 = getDrawCoords( cds1 );
  pair<float,float> c2 = getDrawCoords( cds2 );

  qp_.drawLine( QPointF( c1.first , c1.second ) , QPointF( c2.first , c2.second ) );

}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DQt::drawChar( char c , const std::pair<float,float> &cds ) {

  QRectF br = qp_.boundingRect( 0 , 0 , 100 , 100 , Qt::AlignLeft | Qt::AlignBottom ,
                                QString( c ) );
  qp_.drawText( QRectF( cds.first , cds.second , br.width() , br.height() ) ,
                Qt::AlignLeft | Qt::AlignBottom , QString( c ) , &br );

}

// ****************************************************************************
void MolDraw2DQt::drawTriangle( const pair<float , float> &cds1 ,
                                const pair<float , float> &cds2 ,
                                const pair<float, float> &cds3 ) {

#ifdef NOTYET
  QBrush brush( "Black" );
  brush.setStyle( Qt::SolidPattern );
  DrawColour cc = getColour();
  brush.setColor( QColor( 255.0 * cc.get<0>() , 255.0 * cc.get<1>() , 255.0 * cc.get<2>() ) );
#endif

  qp_.save();
//   qp_.setBrush( brush );

  pair<float,float> c1 = getDrawCoords( cds1 );
  pair<float,float> c2 = getDrawCoords( cds2 );
  pair<float,float> c3 = getDrawCoords( cds3 );

  QPointF points[3] = { QPointF( c1.first , c1.second ) ,
                        QPointF( c2.first , c2.second ) ,
                        QPointF( c3.first , c3.second ) };

  qp_.drawConvexPolygon( points , 3 );

  qp_.restore();

}

// ****************************************************************************
void MolDraw2DQt::clearDrawing() {

  qp_.setBackground( QBrush( "white" ) );
  qp_.fillRect( 0 , 0 , width() , height() , QColor( "white") );

}

// ****************************************************************************
void MolDraw2DQt::setFontSize( float new_size ) {

  MolDraw2D::setFontSize( new_size );
  float font_size_in_points = fontSize() * scale();
#ifdef NOTYET
  cout << "initial font size in points : " << qp_.font().pointSizeF() << endl;
  cout << "font_size_in_points : " << font_size_in_points << endl;
#endif
  QFont font( qp_.font() );
  font.setPointSizeF( font_size_in_points );
  qp_.setFont( font );

  while( 1 ) {
    float old_font_size_in_points = font_size_in_points;
    float font_size_in_points = fontSize() * scale();
    if( fabs( font_size_in_points - old_font_size_in_points ) < 0.1 ) {
      break;
    }
    QFont font( qp_.font() );
    font.setPointSizeF( font_size_in_points );
    qp_.setFont( font );
    calculateScale();
  }

}

// ****************************************************************************
// using the current scale, work out the size of the label in molecule coordinates
void MolDraw2DQt::getStringSize( const string &label , float &label_width ,
                                 float &label_height ) const {

  label_width = 0.0;
  label_height = 0.0;

  int draw_mode = 0; // 0 for normal, 1 for superscript, 2 for subscript
  QString next_char( " " );
  bool had_a_super = false;

  for( int i = 0 , is = label.length() ; i < is ; ++i ) {

    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if( '<' == label[i] && setStringDrawMode( label , draw_mode , i ) ) {
      continue;
    }

    next_char[0] = label[i];
    QRectF br = qp_.boundingRect( 0 , 0 , 100 , 100 ,
                                  Qt::AlignBottom | Qt::AlignLeft , next_char );
    label_height = br.height() / scale();
    float char_width = br.width() / scale();
    if( 2 == draw_mode ) {
      char_width *= 0.5;
    } else if( 1 == draw_mode ) {
      char_width *= 0.5;
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
