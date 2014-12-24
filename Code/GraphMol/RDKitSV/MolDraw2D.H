//
// file MolDraw2D.H
// David Cosgrove
// AstraZeneca
// 27th May 2014
//
// This class makes a 2D drawing of an RDKit molecule.
// It draws heavily on $RDBASE/GraphMol/MolDrawing/MolDrawing.h.
// One purpose of this is to make it easier to overlay annotations on top of
// the molecule drawing, which is difficult to do from the output of
// MolDrawing.h
// The class design philosophy echoes a standard one:
// a virtual base class defines the interface and does all
// the heavy lifting and concrete derived classes implement
// library-specific drawing code such as drawing lines, writing strings
// etc.

#ifndef RDKITMOLDRAW2D_H
#define RDKITMOLDRAW2D_H

#include <vector>

#include <GraphMol/RDKitBase.h>

#include <boost/tuple/tuple.hpp>

// ****************************************************************************

namespace RDKit {

class MolDraw2D {

public :

  typedef enum { C = 0 , N , E , S , W } OrientType;
  typedef boost::tuple<float,float,float> DrawColour;

  MolDraw2D( int width , int height );
  virtual ~MolDraw2D() {}

  virtual void drawMolecule( const ROMol &mol ,
                             const std::vector<int> &highlight_atoms = std::vector<int>() ,
                             const std::map<int,DrawColour> &highlight_map = std::map<int,DrawColour>() );

  // transform a set of coords in the molecule's coordinate system
  // to drawing system coordinates and vice versa. Note that the coordinates have
  // the origin in the top left corner, which is how Qt and Cairo have it, no
  // doubt a holdover from X Windows. This means that a higher y value will be
  // nearer the bottom of the screen. This doesn't really matter except when
  // doing text superscripts and subscripts.
  virtual std::pair<float,float> getDrawCoords( const std::pair<float,float> &mol_cds ) const;
  virtual std::pair<float,float> getDrawCoords( int at_num ) const;
  virtual std::pair<float,float> getAtomCoords( const std::pair<int,int> &screen_cds ) const;
  virtual std::pair<float,float> getAtomCoords( int at_num ) const;


  virtual int width() const { return width_; }
  virtual int height() const { return height_; }
  virtual float scale() const { return scale_; }
  virtual float fontSize() const { return font_size_; }
  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit.
  virtual void setFontSize( float new_size );
  virtual void calculateScale();
  virtual void setColour( const DrawColour &col ) { curr_colour_ = col; }
  virtual DrawColour getColour() const { return curr_colour_; }

  // establishes whether to put string draw mode into super- or sub-script
  // mode based on contents of instring from i onwards. Increments i appropriately
  // and returns true or false depending on whether it did something or not
  bool setStringDrawMode( const std::string &instring , int &draw_mode ,
                             int &i ) const;

private :

  int width_ , height_;
  float scale_;
  float x_min_ , y_min_ , x_range_ , y_range_;
  float x_trans_ , y_trans_;
  // font_size_ in molecule coordinate units. Default 0.5 (a bit bigger
  // than the default width of a double bond)
  float font_size_;
  DrawColour curr_colour_;

  std::vector<std::pair<float,float> > at_cds_; // from mol
  std::vector<int> atomic_nums_;
  std::vector<std::pair<std::string,OrientType> > atom_syms_;

  virtual void drawLine( const std::pair<float,float> &cds1 ,
                         const std::pair<float,float> &cds2 ) = 0;
  virtual void drawLine( const std::pair<float,float> &cds1 ,
                         const std::pair<float,float> &cds2 ,
                         const DrawColour &col1 ,
                         const DrawColour &col2 );
  // draw the char, with the bottom left hand corner at cds
  virtual void drawChar( char c , const std::pair<float,float> &cds ) = 0;

  // drawString centres the string on cds.
  virtual void drawString( const std::string &str ,
                           const std::pair<float,float> &cds );
  // draw a filled triangle
  virtual void drawTriangle( const std::pair<float,float> &cds1 ,
                             const std::pair<float,float> &cds2 ,
                             const std::pair<float,float> &cds3 ) = 0;

  virtual void clearDrawing() = 0;

  // using the current scale, work out the size of the label in molecule coordinates.
  // Bear in mind when implementing this, that, for example, NH2 will appear as
  // NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to accounted
  // for in the width and height.
  virtual void getStringSize( const std::string &label , float &label_width ,
                              float &label_height ) const = 0;

  // return a DrawColour based on the contents of highlight_atoms or
  // highlight_map, falling back to atomic number by default
  DrawColour getColour( int atom_idx ,
                         const std::vector<int> &highlight_atoms ,
                         const std::map<int,DrawColour> &highlight_map );
  DrawColour getColourByAtomicNum( int atomic_num );

  void extractAtomCoords( const ROMol &mol );
  void extractAtomSymbols( const ROMol &mol );

  void drawBond( const ROMol &mol , const BOND_SPTR &bond ,
                  int at1_idx , int at2_idx ,
                  const std::vector<int> &highlight_atoms ,
                  const std::map<int,DrawColour> &highlight_map );
  void drawWedgedBond( const std::pair<float,float> &cds1 ,
                         const std::pair<float,float> &cds2 ,
                         bool draw_dashed , const DrawColour &col1 ,
                         const DrawColour &col2 );
  void drawAtomLabel( int atom_num ,
                        const std::vector<int> &highlight_atoms ,
                        const std::map<int,DrawColour> &highlight_map );

  // calculate normalised perpendicular to vector between two coords
  std::pair<float,float> calcPerpendicular( const std::pair<float,float> &cds1 ,
                                             const std::pair<float,float> &cds2 );
  // cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing into
  // the ring.
  std::pair<float,float> bondInsideRing( const ROMol &mol , const BOND_SPTR &bond ,
                                           const std::pair<float,float> &cds1 ,
                                           const std::pair<float,float> &cds2 );
  // cds1 and cds2 are 2 atoms in a chain double bond.  Returns the perpendicular
  // pointing into the inside of the bond
  std::pair<float,float> bondInsideDoubleBond( const ROMol &mol , const BOND_SPTR &bond );
  // calculate normalised perpendicular to vector between two coords, such that
  // it's inside the angle made between (1 and 2) and (2 and 3).
  std::pair<float,float> calcInnerPerpendicular( const std::pair<float,float> &cds1 ,
                                                   const std::pair<float,float> &cds2 ,
                                                   const std::pair<float,float> &cds3 );

  // take the coords for atnum, with neighbour nbr_cds, and move cds out to accommodate
  // the label associated with it.
  void adjustBondEndForLabel( int atnum , const std::pair<float,float> &nbr_cds ,
                                  std::pair<float,float> &cds ) const;

  // adds LaTeX-like annotation for super- and sub-script.
  std::pair<std::string,OrientType> getAtomSymbolAndOrientation( const Atom &atom ,
                                                                     const std::pair<float,float> &nbr_sum );
};

}

#endif // RDKITMOLDRAW2D_H
