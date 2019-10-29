//
// file qtdraw_demo.cc
// Contributed by:
// David Cosgrove
// AstraZeneca
// 24th July 2012
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Uses RDKit draw funtionality to draw a 2D molecule into a Qt widget

#include <QApplication>
#include <QPainter>
#include <QWidget>

// *****************************************************************************
// in file RDKitMolToQPainter.cc
void smiles_to_qpainter( const std::string &smiles , const QString &label ,
                         int width , int height , bool bow ,
                         QPainter &qp );

// *****************************************************************************

class RDKitDraw : public QWidget {

public :

  RDKitDraw( QWidget *parent = 0 ) : QWidget( parent ) {}

  void set_smiles( const std::string &smi ) {
    smiles_ = smi;
    update();
  }

  void set_name( const std::string &name ) {
    name_ = name;
    update();
  }

protected :

  void paintEvent( QPaintEvent *event ) {

    int w = width();
    int h = height();

    QPainter qp;
    qp.begin( this );
    qp.setBackground( QColor( "White" ) );
    qp.fillRect( 0 , 0 , w , h , QColor( "White" ) );

    if( smiles_.empty() ) {
      qp.end();
      return;
    }

    smiles_to_qpainter( smiles_ , QString( name_.c_str() ) , w , h , true , qp );

    qp.end();

  }

private :

  std::string smiles_;
  std::string name_;

};

// *****************************************************************************
int main( int argc , char **argv ) {

  QApplication a( argc , argv );

  RDKitDraw *rdkd = new RDKitDraw;

  rdkd->setGeometry( 50 , 50 , 500 , 500 );
  rdkd->show();

  if( argc > 1 ) {
    rdkd->set_smiles( argv[1] );
  }
  if( argc > 2 ) {
    rdkd->set_name( argv[2] );
  }
  if( argc == 1 ) {
    rdkd->set_smiles( "[Mg]c1c(C#N)cc(C(=O)NCc2sc([NH3+])c([NH3+])c2)cc1" );
    rdkd->set_name( "Molecule 2" );
  }

  return a.exec();

}
