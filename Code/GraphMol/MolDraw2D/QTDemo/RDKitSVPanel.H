//
// file RDKitSVPanel.H
// David Cosgrove
// AstraZeneca
// 20th June 2014
//
// This class does a 2D display of 1 molecule selected from a vector of them
// using a QSlider.

#ifndef RDKITSVPANEL_H
#define RDKITSVPANEL_H

#include <GraphMol/MolDraw2D/MolDraw2D.h>

#include <QWidget>

// ****************************************************************************

class QLabel;
class QSlider;

namespace RDKit {
class MolDisplay2DWidget;
}

// ****************************************************************************

namespace RDKitSV {

class RDKitSVPanel : public QWidget {

  Q_OBJECT

public :

  RDKitSVPanel( bool left_slider = true , QWidget *parent = 0 , Qt::WindowFlags f = 0 );

  void set_molecules( const std::vector<RDKit::ROMOL_SPTR> &new_mols ,
                      const std::vector<std::vector<int> > &highlight_atoms = std::vector<std::vector<int> >() );
  std::vector<RDKit::ROMOL_SPTR> get_molecules() const { return mols_; }

  void set_label( const QString &new_label );

private :

  RDKit::MolDisplay2DWidget *mol_draw_;
  QSlider *mol_slider_;
  QLabel *label_;

  std::vector<RDKit::ROMOL_SPTR> mols_;
  // highlight_atoms_ is either empty or the same size as mols_, and gives
  // for each molecule a set of atoms to be highlighted.
  std::vector<std::vector<int> > highlight_atoms_;

  void build_widget( bool left_slider );

public slots :

  void slot_slider_changed();

};

} // EO namespace RDKitSV

#endif // RDKITSVPANEL_H
