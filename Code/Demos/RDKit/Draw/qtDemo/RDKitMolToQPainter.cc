
//
// file RDKitMolToQPainter.cc
//
// Contributed by:
// David Cosgrove
// AstraZeneca
// 28th June 2012
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Draw an RDKit Drawing into a QPainter object.
// Hacked together by analogy with the SVG drawer in
// $RDBASE/Code/Demos/RDKit/Draw/demo.cpp
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDrawing/MolDrawing.h> // this is what does the 2D drawing!

#include <QFont>
#include <QPainter>
#include <QString>

#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include <boost/scoped_ptr.hpp>

using namespace boost;
using namespace std;
using namespace RDKit;

typedef scoped_ptr<RWMol> pRWMol;

static const int LINE_WIDTH_MULT = 4;
static const float ORIGINAL_FONT_SIZE = 35.0; // this size seems to work well with the relative scale of the Drawing

// ****************************************************************************
QColor &colour_for_atomic_num( int atnum ) {

  static QColor colours[10] = { QColor( "Black" ) , QColor( "Blue" ) ,
                                QColor( "Red" ) , QColor( "LightGreen" ) ,
                                QColor( "Orange" ) , QColor( "Yellow" ) ,
                                QColor( "Green" ) , QColor( "Brown" ) ,
                                QColor( "Purple" ) , QColor( "Cyan" ) };

  switch( atnum ) {
  case 6 : return colours[0];
  case 7 : return colours[1];
  case 8 : return colours[2];
  case 9 : return colours[3];
  case 15 : return colours[4];
  case 16 : return colours[5];
  case 17 : return colours[6];
  case 37 : return colours[7];
  case 53 : return colours[8];
  default : return colours[9];
  }

}

#ifdef NOTYET
// ****************************************************************************
string colour_for_atomic_num( int atnum ) {

  static map<int,string> colours;
  if( colours.empty() ) {
    colours.insert( make_pair( 6 , "Black" ) ); //carbon
    colours.insert( make_pair( 7 , "Blue" ) ); // nitrogen
    colours.insert( make_pair( 8 , "Red" ) ); // oxygen
    colours.insert( make_pair( 9 , "LightGreen" ) ); // fluorine
    colours.insert( make_pair( 15 , "Orange" ) ); // phosphorous
    colours.insert( make_pair( 16 , "Yellow" ) ); // sulfur
    colours.insert( make_pair( 17 , "Green" ) ); // chlorine
    colours.insert( make_pair( 37 , "Brown" ) ); // bromine
    colours.insert( make_pair( 53 , "Purple" ) ); // iodine
  }

  map<int,string>::iterator q = colours.find( atnum );
  if( q == colours.end() ) {
    return "Cyan";
  } else {
    return q->second;
  }

}
#endif

// ****************************************************************************
void get_scale( int width , int height , int qp_width , int qp_height ,
                qreal &gscale ) {

  qreal xscale = qreal( qp_width ) / qreal( width );
  qreal yscale = qreal( qp_height ) / qreal( height );

  gscale = std::min( xscale , yscale );

#ifdef NOTYET
  cout << "Scale for picture " << width << " by " << height << endl
       << " qpainter " << qp_width << " by " << qp_height << endl
       << "x and y scales : " << xscale << " , " << yscale << endl;
  cout << "Final scale : " << gscale << endl;
#endif

}

// ****************************************************************************
void set_scale( int width , int height , int min_x , int min_y ,
                int qp_width , int qp_height ,
                QPainter &qp , qreal &gscale ) {

  qp.scale( gscale , gscale );
  qreal qp_width_in_scale = qreal( qp_width ) / gscale;
  qreal qp_height_in_scale = qreal( qp_height ) / gscale;

  // move it to a central point in the window
  qp.translate( -min_x + 0.5 * ( qp_width_in_scale - width ) , -min_y + 0.5 * ( qp_height_in_scale - height ) );

}

// ****************************************************************************
void set_font( QPainter &qp ) {

  static QFont *font = 0;
  if( !font ) {
    font = new QFont( qp.font() );
  }

  font->setPointSizeF( ORIGINAL_FONT_SIZE );
  qp.setFont( *font );

}

// ****************************************************************************
void draw_line( vector<int>::iterator &pos , bool bow , QPainter &qp ) {

  int width = *pos++;
  pos++; // skip 'dashed' flag
  int atom1 = *pos++;
  int atom2 = *pos++;

  int xp1 = *pos++;
  int yp1 = *pos++;
  int xp2 = *pos++;
  int yp2 = *pos++;

  static QPen pen( "Black" );
  pen.setJoinStyle( Qt::RoundJoin );
  pen.setWidth( width * LINE_WIDTH_MULT );
  qp.setPen( pen );

  if( atom1 == atom2 ) {
    if( !bow ) {
      pen.setColor( colour_for_atomic_num( atom1 ) );
      qp.setPen( pen );
    }
    qp.drawLine( xp1 , yp1 , xp2 , yp2 );
  } else  {
    int mx = xp1 + ( xp2 - xp1 ) / 2;
    int my = yp1 + ( yp2 - yp1 ) / 2;
    if( !bow ) {
      pen.setColor( colour_for_atomic_num( atom1 ) );
      qp.setPen( pen );
    }
    qp.drawLine( xp1 , yp1 , mx , my );
    if( !bow ) {
      pen.setColor( colour_for_atomic_num( atom2 ) );
      qp.setPen( pen );
    }
    qp.drawLine( mx , my , xp2 , yp2 );
  }
}

// ****************************************************************************
void draw_atom( vector<int>::iterator &pos , bool bow , QPainter &qp ) {

  int atnum = *pos++;
  int xp = *pos++;
  int yp = *pos++;
  int slen = *pos++;
  QString label;
  for(int i = 0 ; i < slen ; ++i ){
    label += static_cast<char>( *pos++ );
  }
  pos++; // move past orient flag, which is not being used

  char letter[2] = { '\0' , '\0' };

  qp.save();

  QRectF whole_br = qp.boundingRect( 0 , 0 , 100 , 100 , Qt::AlignHCenter , label );
  qp.translate( xp - whole_br.width() / 2.0 , yp - whole_br.height() / 2.0 );

  static QPen pen( "Black" );
  qp.setPen( pen );
  pen.setJoinStyle( Qt::RoundJoin );
  if( !bow ) {
    pen.setColor( colour_for_atomic_num( atnum ) );
    qp.setPen( pen );
  }

  // blank out behind the letters, so bonds don't go through where the label is. I think
  // circles work better than rectangles for this, but it means going through it all twice,
  // as the circles overlap so would obscure the letter that came afterwards
  // This assumes the bonds are drawn before the atoms, which is how it is set up in MolDrawing.h.
  static QBrush letter_brush( qp.background() );
  letter_brush.setStyle( Qt::SolidPattern );

  qp.save();
  qp.setPen( qp.background().color() ); // so the ellipse doesn't have a border
  qp.setBrush( letter_brush );

  // do things 1 letter at a time. We needed to stuff the characters into label first, so as to
  // get the bounding rectangle for the whole label.
  for( int i = 0 ; i < slen ; ++i ) {

    letter[0] = label.toLocal8Bit().data()[i];
    QString at_lab( letter );
    QRectF br = qp.boundingRect( 0 , 0 , 1000 , 1000 , Qt::AlignHCenter , at_lab );

    // subscript for numbers, superscript for charge. Stretch the ellipse out so the bond
    // is obscured over the whole of its extent.
    qreal letter_top = 0.0 , letter_height = br.height();
    if( isdigit( letter[0] ) ) {
      qp.translate( 0.0 , br.height() / 2.0 );
      letter_top = -br.height() / 2.0;
      letter_height += br.height() / 2.0;
    } else if( '+' == letter[0] || '-' == letter[0] ) {
      qp.translate( 0.0 , -br.height() / 2.0 );
      letter_height += br.height() / 2.0;
    }
    qp.drawEllipse( QRectF( -2.0 , -2.0 + letter_top , br.width() + 2.0 , letter_height + 2.0 ) );
    qp.translate( br.width() , 0.0 );

    // remove superscript/subscript
    if( isdigit( letter[0] ) ) {
      qp.translate( 0.0 , -br.height() / 2.0 );
    } else if( '+' == letter[0] || '-' == letter[0] ) {
      qp.translate( 0.0 , br.height() / 2.0 );
    }
  }
  qp.setPen( pen );
  qp.restore();

  // now the letters
  for( int i = 0 ; i < slen ; ++i ) {

    letter[0] = label.toLocal8Bit().data()[i];
    QString at_lab( letter );
    QRectF br = qp.boundingRect( 0 , 0 , 1000 , 1000 , Qt::AlignHCenter , at_lab );

    // subscript for numbers, superscript for charge
    if( isdigit( letter[0] ) ) {
      qp.translate( 0.0 , br.height() / 2.0 );
    } else if( '+' == letter[0] || '-' == letter[0] ) {
      qp.translate( 0.0 , -br.height() / 2.0 );
    }
    qp.drawText( QRectF( 0.0 , 0.0 , br.width() , br.height() ) , Qt::AlignHCenter , at_lab , &br );
    qp.translate( br.width() , 0.0 );

    // remove superscript/subscript
    if( isdigit( letter[0] ) ) {
      qp.translate( 0.0 , -br.height() / 2.0 );
    } else if( '+' == letter[0] || '-' == letter[0] ) {
      qp.translate( 0.0 , br.height() / 2.0 );
    }
  }

  qp.restore();

}

// ****************************************************************************
void line_extremes( vector<int>::const_iterator &pos ,
                    int &min_x , int &max_x ,
                    int &min_y , int &max_y ) {

  pos += 4; // past the line width, type and atom numbers
  int xp1 = *pos++;
  int yp1 = *pos++;
  int xp2 = *pos++;
  int yp2 = *pos++;

  min_x = xp1 < min_x ? xp1 : min_x;
  max_x = xp1 > max_x ? xp1 : max_x;
  min_y = yp1 < min_y ? yp1 : min_y;
  max_y = yp1 > max_y ? yp1 : max_y;

  min_x = xp2 < min_x ? xp2 : min_x;
  max_x = xp2 > max_x ? xp2 : max_x;
  min_y = yp2 < min_y ? yp2 : min_y;
  max_y = yp2 > max_y ? yp2 : max_y;

}

// ****************************************************************************
void atom_extremes( vector<int>::const_iterator &pos ,
                    QPainter &qp ,
                    int &min_x , int &max_x ,
                    int &min_y , int &max_y ) {

  pos++; // skip the atomic number, not needed for this
  int xp = *pos++;
  int yp = *pos++;
  int slen = *pos++;
  QString label;
  for(int i = 0 ; i < slen ; ++i ){
    label += static_cast<char>( *pos++ );
  }
  pos++; // move past orient flag, which is not being used, and indeed is not understood!

  QRectF whole_br = qp.boundingRect( 0 , 0 , 1000 , 1000 , Qt::AlignHCenter , label );

  int labx1 = xp - int( whole_br.width() / 2.0 );
  int labx2 = xp + int( whole_br.width() );
  // subscripts and superscripts may occur in the atom label. Just allow for that willy-nilly.
  int laby1 = yp - int( 0.5 * whole_br.height() );
  int laby2 = yp + int( 0.5 * whole_br.height() );

  min_x = labx1 < min_x ? labx1 : min_x;
  max_x = labx2 > max_x ? labx2 : max_x;
  min_y = laby1 < min_y ? laby1 : min_y;
  max_y = laby2 > max_y ? laby2 : max_y;

}

// ****************************************************************************
void setup_scale( vector<int> &drawing , QPainter &qp ,
                  int qp_width , int qp_height ) {

  // calculate the width of the drawing. The values for width and height as
  // supplied in drawing by MolDrawing.h don't allow for the width of atom labels
  // so if there are non-carbon atoms at the extremes of the picture, these
  // will be drawn off the page.
  vector<int>::const_iterator pos = drawing.begin() + 7;
  int min_x = numeric_limits<int>::max();
  int max_x = numeric_limits<int>::min();
  int min_y = numeric_limits<int>::max();
  int max_y = numeric_limits<int>::min();

  while( pos != drawing.end() ) {
    int token = *pos++;
    switch( token ) {
    case Drawing::LINE :
      line_extremes( pos , min_x , max_x , min_y , max_y );
      break;
    case Drawing::ATOM :
      atom_extremes( pos , qp , min_x , max_x , min_y , max_y );
      break;
    default :
      break;
    }
  }

  // allow for width of lines
  min_x -= LINE_WIDTH_MULT;
  max_x += LINE_WIDTH_MULT;
  min_y -= LINE_WIDTH_MULT;
  max_y += LINE_WIDTH_MULT;

  int width = max_x - min_x;
  int height = max_y - min_y;

  qreal gscale;
  get_scale( width , height , qp_width , qp_height , gscale );
  set_scale( width , height , min_x , min_y , qp_width , qp_height , qp , gscale );

}

// ****************************************************************************
void MolToQPainter( vector<int> &drawing , int qp_width , int qp_height ,
                    bool bow , QPainter &qp ) {

  vector<int>::iterator pos = drawing.begin() + 2;
  if( *pos != Drawing::BOUNDS ) {
    // no drawing made by RDKit
    QString msg( "Bad Drawing" );
    qp.drawText( 0 , qp_height / 2 , msg );
    return;
  }
  pos += 3;

  // save current tranformation matrix()
  qp.save();

  set_font( qp );
  setup_scale( drawing , qp, qp_width , qp_height );

  pos += 2; // for the width and height in the drawing
  while( pos != drawing.end() ) {
    int token = *pos++;
    switch( token ) {
    case Drawing::LINE :
      draw_line( pos , bow , qp );
      break;
    case Drawing::ATOM :
      draw_atom( pos , bow , qp );
      break;
    default :
      cerr << "Duff token, expect a crash." << token << endl;
      break;
    }
  }

  // put the transformation matrix back to how it was on entry
  qp.restore();

}

// ****************************************************************************
void label_molecule( const QString &label , int width , int height ,
                     int label_height , QPainter &qp ) {

  int label_y = height - label_height;
  qp.drawText( QRect( 0 , label_y , width , label_height ) , Qt::AlignHCenter , label );

}

// ****************************************************************************
// width and height are the size of the rectangle in qp we are drawing into.
// bow is black on white- when false,  colours are used in drawing
void smiles_to_qpainter( const string &smiles , const QString &label ,
                         int width , int height , bool bow ,
                         QPainter &qp ) {

  // render into an RDKit drawing
  pRWMol mol( SmilesToMol( smiles ) );
  MolOps::Kekulize( *mol );
  RDDepict::compute2DCoords( *mol );
  vector<int> drawing = Drawing::DrawMol( *mol );

  qp.setRenderHint( QPainter::Antialiasing , true );
  qp.setRenderHint( QPainter::TextAntialiasing , true );

  // get the height that will be needed for the label, and allow the picture to accommodate it
  QRect br;
  if( !label.isEmpty() ) {
    br = qp.boundingRect( QRect( 0 , 0 , width , height ) , Qt::AlignHCenter , label );
  }

  int label_height = br.height() + 5;
  int pict_height = height - label_height;
  MolToQPainter( drawing , width , pict_height , bow , qp );

  if( !label.isEmpty() ) {
    label_molecule( label , width , height , label_height , qp );
  }

}
