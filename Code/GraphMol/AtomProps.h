//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
/*! \file AtomProps.h

  \brief No user-serviceable parts inside.

  Contains internals for using Atoms with the BGL

*/  
#ifndef _RD_ATOM_PROPS_H
#define _RD_ATOM_PROPS_H

// boost stuff
#include <boost/property_map.hpp>
#include <boost/smart_ptr.hpp>

namespace RDKit{
  class Atom;
  typedef boost::shared_ptr<Atom>    ATOM_SPTR;
  
  struct vertex_atom_t {
    enum { num=1001 };
    typedef boost::vertex_property_tag kind;
  };

  typedef boost::property<vertex_atom_t,Atom *> AtomProperty;

};


#endif
