//
// file MolDraw2D.cc
// David Cosgrove
// AstraZeneca
// 27th May 2014
//

#include "MolDraw2D.H"

#include <cstdlib>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDraw2D::MolDraw2D( int width, int height ) :
width_( width ) , height_( height ) , scale_( 1.0 ) , x_trans_( 0.0 ) ,
  y_trans_( 0.0 ) , font_size_( 0.5 ) {

}

// ****************************************************************************
void MolDraw2D::DrawMolecule( const ROMol &mol ,
                              const vector<int> &highlight_atoms ,
                              const map<int,DrawColour> &highlight_map ) {

  clearDrawing();
  extract_atom_coords( mol );
  extract_atom_symbols( mol );
  calculateScale();
  setFontSize( font_size_ );

  ROMol::VERTEX_ITER this_at , end_at;
  boost::tie( this_at , end_at ) = mol.getVertices();
  while( this_at != end_at ) {
    int this_idx = mol[*this_at]->getIdx();
    ROMol::OEDGE_ITER nbr , end_nbr;
    boost::tie( nbr , end_nbr ) = mol.getAtomBonds( mol[*this_at].get() );
    while( nbr != end_nbr ) {
      const BOND_SPTR bond = mol[*nbr];
      ++nbr;
      int nbr_idx = bond->getOtherAtomIdx( this_idx );
      if( nbr_idx < static_cast<int>( at_cds_.size() ) && nbr_idx > this_idx ) {
        draw_bond( mol , bond , this_idx , nbr_idx , highlight_atoms ,
                   highlight_map );
      }
    }
    ++this_at;
  }

  for( int i = 0 , is = atom_syms_.size() ; i < is ; ++i ) {
    if( !atom_syms_[i].first.empty() ) {
      draw_atom_label( i , highlight_atoms , highlight_map );
    }
  }

}

// ****************************************************************************
// transform a set of coords in the molecule's coordinate system
// to drawing system coordinates
pair<float,float> MolDraw2D::getDrawCoords( const std::pair<float, float> &mol_cds ) const {

  float x = scale_ * ( mol_cds.first - x_min_ + x_trans_ );
  float y = scale_ * ( mol_cds.second - y_min_ + y_trans_ );

  return make_pair( x , y );

}

// ****************************************************************************
pair<float,float> MolDraw2D::getDrawCoords( int at_num ) const {

  return getDrawCoords( at_cds_[at_num] );

}

// ****************************************************************************
pair<float,float> MolDraw2D::getAtomCoords( const pair<int, int> &screen_cds) const {

  int x = int( float( screen_cds.first ) / scale_ + x_min_ - x_trans_ );
  int y = int( float( screen_cds.second ) / scale_ + y_min_ - y_trans_ );

  return make_pair( x , y );

}

// ****************************************************************************
pair<float,float> MolDraw2D::getAtomCoords( int at_num ) const {

  return at_cds_[at_num];

}

// ****************************************************************************
void MolDraw2D::setFontSize( float new_size ) {

  font_size_ = new_size;

}

// ****************************************************************************
void MolDraw2D::calculateScale() {

  x_min_ = y_min_ = numeric_limits<float>::max();
  float x_max( -numeric_limits<float>::max() ) , y_max( -numeric_limits<float>::max() );
  for( int i = 0 , is = at_cds_.size() ; i < is ; ++i ) {
    const pair<float,float> &pt = at_cds_[i];
    x_min_ = std::min( pt.first , x_min_ );
    y_min_ = std::min( pt.second , y_min_ );
    x_max = std::max( pt.first , x_max );
    y_max = std::max( pt.second , y_max );
  }

  x_range_ = x_max - x_min_;
  y_range_ = y_max - y_min_;
  scale_ = std::min( float( width_ ) / x_range_ , float( height_ ) / y_range_ );

  // we may need to adjust the scale if there are atom symbols that go off
  // the edges, and we probably need to do it iteratively because get_string_size
  // uses the current value of scale_.
  while( 1 ) {
    for( int i = 0 , is = atom_syms_.size() ; i < is ; ++i ) {
      if( !atom_syms_[i].first.empty() ) {
        float atsym_width , atsym_height;
        getStringSize( atom_syms_[i].first , atsym_width , atsym_height );
        float this_x = at_cds_[i].first;
        float this_y = at_cds_[i].second + atsym_height;
        if( W == atom_syms_[i].second ) {
          this_x -= atsym_width;
        } else if( E == atom_syms_[i].second ) {
          this_x += atsym_width;
        }
        x_max = std::max( x_max , this_x );
        x_min_ = std::min( x_min_ , this_x );
        y_max = std::max( y_max , this_y );
      }
    }
    float old_scale = scale_;
    x_range_ = x_max - x_min_;
    y_range_ = y_max - y_min_;
    scale_ = std::min( float( width_ ) / x_range_ , float( height_ ) / y_range_ );
    if( fabs( scale_ - old_scale ) < 0.1 ) {
      break;
    }
  }

  // put a 5% buffer round the drawing and calculate a final scale
  x_min_ -= 0.05 * x_range_;
  x_range_ *= 1.1;
  y_min_ -= 0.05 * y_range_;
  y_range_ *= 1.1;

  scale_ = std::min( float( width_ ) / x_range_ , float( height_ ) / y_range_ );

  float x_mid = x_min_ + 0.5 * x_range_;
  float y_mid = y_min_ + 0.5 * y_range_;
  pair<float,float> mid = getDrawCoords( make_pair( x_mid , y_mid ) );
  x_trans_ = ( width_ / 2 - mid.first ) / scale_;
  y_trans_ = ( height_ / 2 - mid.second ) / scale_;

}

// ****************************************************************************
// establishes whether to put string draw mode into super- or sub-script
// mode based on contents of instring from i onwards. Increments i appropriately
// and returns true or false depending on whether it did something or not.
bool MolDraw2D::set_string_draw_mode( const string &instring , int &draw_mode ,
                                      int &i ) const {

  string bit1 = instring.substr( i , 5 );
  string bit2 = instring.substr( i , 6 );

  // could be markup for super- or sub-script
  if( string( "<sub>" ) == bit1 ) {
    draw_mode = 2;
    i += 4;
    return true;
  } else if( string( "<sup>" ) == bit1 ) {
    draw_mode = 1;
    i += 4;
    return true;
  } else if( string( "</sub>") == bit2 ) {
    draw_mode = 0;
    i += 5;
    return true;
  } else if( string( "</sup>") == bit2 ) {
    draw_mode = 0;
    i += 5;
    return true;
  }

  return false;

}

// ****************************************************************************
void MolDraw2D::drawLine( const pair<float, float> &cds1 ,
                          const pair<float, float> &cds2 ,
                          const DrawColour &col1 , const DrawColour &col2) {

  if( col1 == col2 ) {
    setColour( col1 );
    drawLine( cds1 , cds2 );
  } else {
    pair<float,float> mid( 0.5 * ( cds1.first + cds2.first ) ,
                           0.5 * ( cds1.second + cds2.second ) );
    setColour( col1 );
    drawLine( cds1 , mid );
    setColour( col2 );
    drawLine( mid , cds2 );
  }

}

// ****************************************************************************
// draws the string centred on cds
void MolDraw2D::drawString( const string &str ,
                            const pair<float,float> &cds ) {

  float string_width , string_height;
  getStringSize( str , string_width , string_height );

  float draw_x = cds.first - string_width / 2.0;
  float draw_y = cds.second - string_height / 2.0;

  float full_font_size = fontSize();
  int draw_mode = 0; // 0 for normal, 1 for superscript, 2 for subscript
  string next_char( " " );

  for( int i = 0 , is = str.length() ; i < is ; ++i ) {

    // set_string_draw_mode moves i along to the end of any <sub> or <sup>
    // markup
    if( '<' == str[i] && set_string_draw_mode( str , draw_mode , i ) ) {
      continue;
    }

    char next_c = str[i];
    next_char[0] = next_c;
    float char_width , char_height;
    getStringSize( next_char , char_width , char_height );

    // these font sizes and positions work best for Qt, IMO. They may want
    // tweaking for a more general solution.
    if( 2 == draw_mode ) {
      // y goes from top to bottom, so add for a subscript!
      setFontSize( 0.5 * full_font_size );
      char_width *= 0.5;
      drawChar( next_c , getDrawCoords( make_pair( draw_x , draw_y + 0.5 * char_height ) ) );
      setFontSize( full_font_size );
    } else if( 1 == draw_mode ) {
      setFontSize( 0.5 * full_font_size );
      char_width *= 0.5;
      drawChar( next_c , getDrawCoords( make_pair( draw_x , draw_y - 0.25 * char_height ) ) );
      setFontSize( full_font_size );
    } else {
      drawChar( next_c , getDrawCoords( make_pair( draw_x , draw_y ) ) );
    }
    draw_x += char_width;
  }

}

// ****************************************************************************
MolDraw2D::DrawColour MolDraw2D::get_colour( int atom_idx ,
                                             const std::vector<int> &highlight_atoms ,
                                             const std::map<int,DrawColour> &highlight_map ) {

  DrawColour retval = get_colour_by_atomic_num( atomic_nums_[atom_idx] );

  // set contents of highlight_atoms to red
  if( highlight_atoms.end() != find( highlight_atoms.begin() ,
                                     highlight_atoms.end() , atom_idx )  ) {
    retval = DrawColour( 1.0 , 0.0 , 0.0 );
  }
  // over-ride with explicit colour from highlight_map if there is one
  map<int,DrawColour>::const_iterator p = highlight_map.find( atom_idx );
  if( p != highlight_map.end() ) {
    retval = p->second;
  }

  return retval;

}

// ****************************************************************************
MolDraw2D::DrawColour MolDraw2D::get_colour_by_atomic_num( int atomic_num ) {

  // RGB values taken from Qt's QColor. The seem to work pretty well on my
  // machine. Using them as fractions of 255, as that's the way Cairo does it.
  float this_col[3] = { 0 , 1.0 , 1.0 }; // default to cyan

  switch( atomic_num ) {
  case 6 : this_col[0] = 0.0; this_col[1] = 0.0; this_col[2] = 0.0; break;
  case 7 : this_col[0] = 0.0; this_col[1] = 0.0; this_col[2] = 1.0; break;
  case 8 : this_col[0] = 1.0; this_col[1] = 0.0; this_col[2] = 0.0; break;
  case 9 : this_col[0] = 0.565; this_col[1] = 0.933; this_col[2] = 0.565; break;
  case 15 : this_col[0] = 1.0; this_col[1] = 0.647; this_col[2] = 0.0; break;
  case 16 : this_col[0] = 1.0; this_col[1] = 1.0; this_col[2] = 0.0; break;
  case 17 : this_col[0] = 0.0; this_col[1] = 0.502; this_col[2] = 0.0; break;
  case 37 : this_col[0] = 0.647; this_col[1] = 0.165; this_col[2] = 0.165; break;
  case 53 : this_col[0] = 0.502; this_col[1] = 0.0; this_col[2] = 0.502; break;
  default : break;
  }

  return DrawColour( this_col[0] , this_col[1] , this_col[2] );

}

// ****************************************************************************
void MolDraw2D::extract_atom_coords( const ROMol &mol ) {

  at_cds_.clear();
  atomic_nums_.clear();
  const RDGeom::POINT3D_VECT &locs = mol.getConformer( -1 ).getPositions();
  ROMol::VERTEX_ITER this_at , end_at;
  tie( this_at , end_at ) = mol.getVertices();
  while( this_at != end_at ) {
    int this_idx = mol[*this_at]->getIdx();
    at_cds_.push_back( make_pair( locs[this_idx].x , locs[this_idx].y ) );
    ++this_at;
  }

}

// ****************************************************************************
void MolDraw2D::extract_atom_symbols( const ROMol &mol ) {

  ROMol::VERTEX_ITER atom , end_atom;
  tie( atom , end_atom ) = mol.getVertices();
  while( atom != end_atom ) {
    ROMol::OEDGE_ITER nbr , end_nbrs;
    const Atom *at1 = mol[*atom].get();
    tie( nbr , end_nbrs ) = mol.getAtomBonds( at1 );
    pair<float,float> &at1_cds = at_cds_[at1->getIdx()];
    pair<float,float> nbr_sum( 0.0 , 0.0 );
    while( nbr != end_nbrs ) {
      const BOND_SPTR bond = mol[*nbr];
      ++nbr;
      pair<float,float> &at2_cds = at_cds_[bond->getOtherAtomIdx( at1->getIdx() )];
      nbr_sum.first += at2_cds.first - at1_cds.first;
      nbr_sum.second += at2_cds.second - at1_cds.second;
    }
    atom_syms_.push_back( get_atom_symbol_and_orientation( *at1 , nbr_sum ) );
    atomic_nums_.push_back( at1->getAtomicNum() );
    ++atom;
  }
}

// ****************************************************************************
void MolDraw2D::draw_bond( const ROMol &mol , const BOND_SPTR &bond ,
                           int at1_idx , int at2_idx ,
                           const vector<int> &highlight_atoms ,
                           const map<int,DrawColour> &highlight_map ) {

  const Atom *at1 = mol.getAtomWithIdx( at1_idx );
  const Atom *at2 = mol.getAtomWithIdx( at2_idx );
  const float double_bond_offset = 0.1;

  pair<float,float> at1_cds = at_cds_[at1_idx];
  pair<float,float> at2_cds = at_cds_[at2_idx];

  adjust_bond_end_for_label( at1_idx , at2_cds , at1_cds );
  adjust_bond_end_for_label( at2_idx , at1_cds , at2_cds );

  DrawColour col1 = get_colour( at1_idx , highlight_atoms , highlight_map );
  DrawColour col2 = get_colour( at2_idx , highlight_atoms , highlight_map );

  // it it's a double bond and one end is 1-connected, do two lines parallel
  // to the atom-atom line.
  if( ( bond->getBondType() == Bond::DOUBLE ) &&
      ( 1 == at1->getDegree() || 1 == at2->getDegree() ) ) {
    pair<float,float> perp = calc_perpendicular( at1_cds , at2_cds );
    drawLine( make_pair( at1_cds.first + double_bond_offset * perp.first ,
                         at1_cds.second + double_bond_offset * perp.second ) ,
              make_pair( at2_cds.first + double_bond_offset * perp.first ,
                         at2_cds.second + double_bond_offset * perp.second ) , col1 , col2 );
    drawLine( make_pair( at1_cds.first - double_bond_offset * perp.first ,
                         at1_cds.second - double_bond_offset * perp.second ) ,
              make_pair( at2_cds.first - double_bond_offset * perp.first ,
                         at2_cds.second - double_bond_offset * perp.second ) , col1 , col2 );
    if( bond->getBondType() == Bond::TRIPLE ) {
      drawLine( at1_cds , at2_cds , col1 , col2 );
    }
  } else if( Bond::SINGLE == bond->getBondType() &&
             ( Bond::BEGINWEDGE == bond->getBondDir() || Bond::BEGINDASH == bond->getBondDir() ) ) {
    if( bond->getBeginAtom()->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
        bond->getBeginAtom()->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW ) {
      swap( at1_cds , at2_cds );
      swap( col1 , col2 );
    }
    if( Bond::BEGINWEDGE == bond->getBondDir() ) {
      draw_wedged_bond( at1_cds , at2_cds , false , col1 , col2 );
    } else {
      draw_wedged_bond( at1_cds , at2_cds , true , col1 , col2 );
    }
  } else {
    // in all other cases, we will definitely want to draw a line between the
    // two atoms
    drawLine( at1_cds , at2_cds , col1 , col2 );
    if( Bond::TRIPLE == bond->getBondType() ) {
      // 2 further lines, a bit shorter and offset on the perpendicular
      pair<float,float> perp = calc_perpendicular( at1_cds , at2_cds );
      float dbo = 2.0 * double_bond_offset;
      float end1_trunc = 1 == at1->getDegree() ? 0.0 : 0.1;
      float end2_trunc = 1 == at2->getDegree() ? 0.0 : 0.1;
      float bv[2] = { at1_cds.first - at2_cds.first , at1_cds.second - at2_cds.second };
      float px1 = at1_cds.first - end1_trunc * bv[0] + dbo * perp.first;
      float py1 = at1_cds.second - end1_trunc * bv[1] + dbo * perp.second;
      float px2 = at2_cds.first + end2_trunc * bv[0] + dbo * perp.first;
      float py2 = at2_cds.second + end2_trunc * bv[1] + dbo * perp.second;
      drawLine( make_pair( px1 , py1 ) , make_pair( px2 , py2 ) , col1 , col2 );
      px1 = at1_cds.first - end1_trunc * bv[0] - dbo * perp.first;
      py1 = at1_cds.second - end1_trunc * bv[1] - dbo * perp.second;
      px2 = at2_cds.first + end2_trunc * bv[0] - dbo * perp.first;
      py2 = at2_cds.second + end2_trunc * bv[1] - dbo * perp.second;
      drawLine( make_pair( px1 , py1 ) , make_pair( px2 , py2 ) , col1 , col2 );
    }
    // all we have left now are double bonds in a ring or not in a ring
    // and multiply connected
    if( Bond::DOUBLE == bond->getBondType() ) {
      pair<float,float> perp;
      if( mol.getRingInfo()->numBondRings( bond->getIdx() ) ) {
        // in a ring, we need to draw the bond inside the ring
        perp = bond_inside_ring( mol , bond , at1_cds , at2_cds );
      } else {
        perp = bond_inside_double_bond( mol , bond );
      }
      float dbo = 2.0 * double_bond_offset;
      float bv[2] = { at1_cds.first - at2_cds.first , at1_cds.second - at2_cds.second };
      float px1 = at1_cds.first - 0.1 * bv[0] + dbo * perp.first;
      float py1 = at1_cds.second - 0.1 * bv[1] + dbo * perp.second;
      float px2 = at2_cds.first + 0.1 * bv[0] + dbo * perp.first;
      float py2 = at2_cds.second + 0.1 * bv[1] + dbo * perp.second;
      drawLine( make_pair( px1 , py1 ) , make_pair( px2 , py2 ) , col1 , col2 );
    }
  }

}

// ****************************************************************************
void MolDraw2D::draw_wedged_bond( const pair<float, float> &cds1 ,
                                  const pair<float, float> &cds2 ,
                                  bool draw_dashed , const DrawColour &col1 ,
                                  const DrawColour &col2) {

  pair<float,float> perp = calc_perpendicular( cds1 , cds2 );
  pair<float,float> disp( 0.1 * perp.first , 0.1 * perp.second );
  pair<float,float> end1 , end2;
  end1.first = cds2.first + disp.first;
  end1.second = cds2.second + disp.second;
  end2.first = cds2.first - disp.first;
  end2.second = cds2.second - disp.second;

  setColour( col1 );
  if( draw_dashed ) {
    pair<float,float> e1( end1.first - cds1.first , end1.second - cds1.second );
    pair<float,float> e2( end2.first - cds1.first , end2.second - cds1.second );
    for( int i = 1 ; i < 11 ; ++i ) {
      if( 5 == i ) {
        setColour( col2 );
      }
      pair<float,float> e11( cds1.first + float( i ) * 0.1 * e1.first ,
                             cds1.second + float( i ) * 0.1 * e1.second );
      pair<float,float> e22( cds1.first + float( i ) * 0.1 * e2.first ,
                             cds1.second + float( i ) * 0.1 * e2.second );
      drawLine( e11 , e22 );
    }
  } else {
    if( col1 == col2 ) {
      drawTriangle( cds1 , end1 , end2 );
    } else {
      pair<float,float> e1( end1.first - cds1.first , end1.second - cds1.second );
      pair<float,float> e2( end2.first - cds1.first , end2.second - cds1.second );
      pair<float,float> mid1( cds1.first + 0.5 * e1.first ,
                              cds1.second + 0.5 * e1.second );
      pair<float,float> mid2( cds1.first + 0.5 * e2.first ,
                              cds1.second + 0.5 * e2.second );
      drawTriangle( cds1 , mid1 , mid2 );
      setColour( col2 );
      drawTriangle( mid1 , end2 , end1 );
      drawTriangle( mid1 , mid2 , end2 );
    }
  }

}

// ****************************************************************************
void MolDraw2D::draw_atom_label( int atom_num ,
                                 const std::vector<int> &highlight_atoms ,
                                 const std::map<int,DrawColour> &highlight_map ) {

  setColour( get_colour( atom_num , highlight_atoms , highlight_map ) );
  drawString( atom_syms_[atom_num].first , at_cds_[atom_num] );

}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords
pair<float,float> MolDraw2D::calc_perpendicular( const pair<float,float> &cds1 ,
                                                 const pair<float,float> &cds2 ) {

  float bv[2] = { cds1.first - cds2.first , cds1.second - cds2.second };
  float perp[2] = { -bv[1] , bv[0] };
  float perp_len = sqrt( perp[0] * perp[0] + perp[1] * perp[1] );
  perp[0] /= perp_len; perp[1] /= perp_len;

  return make_pair( perp[0] , perp[1] );

}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing into
// the ring
pair<float,float> MolDraw2D::bond_inside_ring( const ROMol &mol , const BOND_SPTR &bond ,
                                               const pair<float,float> &cds1 ,
                                               const pair<float,float> &cds2 ) {

  Atom *bgn_atom = bond->getBeginAtom();
  ROMol::OEDGE_ITER nbr2 , end_nbrs2;
  tie( nbr2 , end_nbrs2 ) = mol.getAtomBonds( bgn_atom );
  while( nbr2 != end_nbrs2 ) {
    const BOND_SPTR bond2 = mol[*nbr2];
    ++nbr2;
    if( bond2->getIdx() == bond->getIdx() ||
        !mol.getRingInfo()->numBondRings( bond2->getIdx() ) ) {
      continue;
    }
    bool same_ring = false;
    BOOST_FOREACH( const INT_VECT &ring , mol.getRingInfo()->bondRings() ) {
      if( find( ring.begin() , ring.end() , bond->getIdx() ) != ring.end() &&
          find( ring.begin() , ring.end() , bond2->getIdx() ) != ring.end() ) {
        same_ring = true;
        break;
      }
    }
    if( same_ring ) {
      // bond and bond2 are in the same ring, so use their vectors to define
      // the sign of the perpendicular.
      int atom3 = bond2->getOtherAtomIdx( bond->getBeginAtomIdx() );
      return calc_inner_perpendicular( cds1 , cds2 , at_cds_[atom3] );
    }
  }

  return calc_perpendicular( cds1 , cds2 );

}

// ****************************************************************************
// cds1 and cds2 are 2 atoms in a chain double bond.  Returns the perpendicular
// pointing into the inside of the bond
pair<float,float> MolDraw2D::bond_inside_double_bond( const ROMol &mol , const BOND_SPTR &bond ) {

  // a chain double bond, were it looks nicer IMO if the 2nd line is inside
  // the angle of outgoing bond. Unless it's an allene, where nothing
  // looks great.
  const Atom *at1 = bond->getBeginAtom();
  const Atom *at2 = bond->getEndAtom();
  const Atom *bond_atom , *end_atom;
  if( at1->getDegree() > 1 ) {
    bond_atom = at1;
    end_atom = at2;
  } else {
    bond_atom = at2;
    end_atom = at1;
  }
  int at3 = -1; // to stop the compiler whinging.
  ROMol::OEDGE_ITER nbr2 , end_nbrs2;
  tie( nbr2 , end_nbrs2 ) = mol.getAtomBonds( bond_atom );
  while( nbr2 != end_nbrs2 ) {
    const BOND_SPTR bond2 = mol[*nbr2];
    ++nbr2;
    if( bond != bond2 ) {
      at3 = bond2->getOtherAtomIdx( bond_atom->getIdx() );
      break;
    }
  }

  return calc_inner_perpendicular( at_cds_[end_atom->getIdx()] ,
      at_cds_[bond_atom->getIdx()] , at_cds_[at3] );

}

// ****************************************************************************
// calculate normalised perpendicular to vector between two coords, such that
// it's inside the angle made between (1 and 2) and (2 and 3).
pair<float,float> MolDraw2D::calc_inner_perpendicular( const pair<float,float> &cds1 ,
                                                       const pair<float,float> &cds2 ,
                                                       const pair<float,float> &cds3 ) {

  pair<float,float> perp = calc_perpendicular( cds1 , cds2 );
  float v1[2] = { cds1.first - cds2.first , cds1.second - cds2.second };
  float v2[2] = { cds2.first - cds3.first , cds2.second - cds3.second };
  float obv[2] = { v1[0] - v2[0] , v1[1] - v2[1] };

  // if dot product of centre_dir and perp < 0.0, they're pointing in opposite
  // directions, so reverse perp
  if( obv[0] * perp.first + obv[1] * perp.second < 0.0 ) {
    perp.first *= -1.0;
    perp.second *= - 1.0;
  }

  return perp;

}

// ****************************************************************************
// take the coords for atnum, with neighbour nbr_cds, and move cds out to accommodate
// the label associated with it.
void MolDraw2D::adjust_bond_end_for_label( int atnum , const std::pair<float,float> &nbr_cds ,
                                           std::pair<float,float> &cds ) const {

  if( atom_syms_[atnum].first.empty() ) {
    return;
  }

  float label_width , label_height;
  getStringSize( atom_syms_[atnum].first , label_width , label_height );

  float lw2 = label_width / 2.0;
  float lh2 = label_height / 2.0;

  float x_offset = 0.0 , y_offset = 0.0;
  if( fabs( nbr_cds.second - cds.second ) < 1.0e-5 ) {
    // if the bond is horizontal
    x_offset = lw2;
  } else {
    x_offset = fabs( lh2 * ( nbr_cds.first - cds.first ) / ( nbr_cds.second - cds.second ) );
    if( x_offset >= lw2 ) {
      x_offset = lw2;
    }
  }
  if( nbr_cds.first < cds.first ) {
    x_offset *= -1.0;
  }

  if( fabs( nbr_cds.first - cds.first ) < 1.0e-5 ) {
    // if the bond is vertical
    y_offset = lh2;
  } else {
    y_offset = fabs( lw2 * ( cds.second - nbr_cds.second ) / ( nbr_cds.first - cds.first ) );
    if( y_offset >= lh2 ) {
      y_offset = lh2;
    }
  }
  if( nbr_cds.second < cds.second ) {
    y_offset *= -1.0;
  }

  cds.first += x_offset;
  cds.second += y_offset;

}

// ****************************************************************************
// adds XML-like annotation for super- and sub-script, in the same manner
// as MolDrawing.py. My first thought was for a LaTeX-like system, obviously...
pair<string,MolDraw2D::OrientType> MolDraw2D::get_atom_symbol_and_orientation( const Atom &atom ,
                                                                               const pair<float, float> &nbr_sum ) {

  string symbol( "" );
  OrientType orient = C;

  if( 6 != atom.getAtomicNum() ) {
    symbol = atom.getSymbol();
  }

  if( 0 != atom.getIsotope() ) {
    symbol = lexical_cast<string>( atom.getIsotope() ) + symbol;
  }

  if( atom.hasProp( "molAtomMapNumber" ) ) {
    string map_num;
    atom.getProp( "molAtomMapNumber" , map_num );
    symbol += string( ":" ) + map_num;
  }

  int num_h = 6 == atom.getAtomicNum() ? 0 : atom.getTotalNumHs();
  if( num_h > 0 ) {
    string h( "H" );
    if( num_h > 1 ) {
      // put the number as a subscript
      h += string( "<sub>" ) + lexical_cast<string>( num_h ) + string( "</sub>" );
    }
    symbol += h;
  }

  if( 0 != atom.getFormalCharge() ) {
    int chg = atom.getFormalCharge();
    string sgn = chg > 0 ? string( "+" ) : string( "-" );
    chg = abs( chg );
    if( chg > 1 ) {
      sgn += lexical_cast<string>( chg );
    }
    // put the charge as a superscript
    symbol += string( "<sup>" ) + sgn + string( "</sup>" );
  }

  if( 1 == atom.getDegree() ) {
    float islope = 0.0;
    if( fabs( nbr_sum.second ) > 1.0 ) {
      islope = nbr_sum.first / fabs( nbr_sum.second );
    } else {
      islope = nbr_sum.first;
    }
    if( fabs( islope ) > 0.85 ) {
      if( islope > 0.0 ) {
        orient = W;
      } else {
        orient = E;
      }
    } else {
      if( nbr_sum.second > 0.0 ) {
        orient = N;
      } else {
        orient = S;
      }
    }
  }

  return make_pair( symbol , orient );

}

} // EO namespace RDKit
