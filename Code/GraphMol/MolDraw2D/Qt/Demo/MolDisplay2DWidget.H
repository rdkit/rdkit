//
// file MolDisplay2DWidget.H
// David Cosgrove
// AstraZeneca
// 6th June 2014
//
// Qt Widget that draws am RDKit molecule using MolDraw2D.

#ifndef MOLDISPLAY2DWIDGET_H
#define MOLDISPLAY2DWIDGET_H

#include <QWidget>

#include <boost/shared_ptr.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

// ****************************************************************************

class QColor;
class QPaintEvent;
class QPainter;

// ****************************************************************************

namespace RDKit {

class MolDraw2D;

class MolDisplay2DWidget : public QWidget {

  Q_OBJECT

public :

  MolDisplay2DWidget( QWidget *parent = 0 , Qt::WindowFlags f = 0);

  void set_display_mol( ROMOL_SPTR new_mol );
  ROMOL_SPTR display_mol() { return disp_mol_; }
  boost::shared_ptr<MolDraw2D> mol_drawer() { return mol_drawer_; }
  const boost::shared_ptr<MolDraw2D> mol_drawer() const { return mol_drawer_; }

  void set_atom_picking( bool new_val );
  std::vector<int> selected_atoms() const { return picked_atoms_; }
  void set_selected_atoms( const std::vector<int> &sa );

  int pick_circle_rad() const;

  QSize minimumSize() const;
  QSize sizeHint() const;

protected :

  void paintEvent( QPaintEvent *event );
  void mousePressEvent( QMouseEvent *event );

  void select_atom( QMouseEvent *event );
  void draw_molecule( QPainter &qp );
  void add_molecule_title( QPainter &qp , const std::string &mol_name ,
                           int label_box_height );
  void identify_selected_atoms( QPainter &qp );

private :

  RWMOL_SPTR disp_mol_;
  boost::shared_ptr<MolDraw2D> mol_drawer_;
  mutable int pick_circle_rad_;
  bool atom_picking_;
  std::vector<int> picked_atoms_;

  int find_nearest_atom( int x_screen_pos , int y_screen_pos ) const;

};

} // EO namespace RDKit

#endif // MOLDISPLAY2DWIDGET_H
