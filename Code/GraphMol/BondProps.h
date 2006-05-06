//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
/*! \file BondProps.h

  \brief No user-serviceable parts inside.

  Contains internals for using Bonds with the BGL

*/  
#ifndef _RD_BOND_PROPS_H
#define _RD_BOND_PROPS_H

// boost stuff
#include <boost/property_map.hpp>
#include <boost/smart_ptr.hpp>

#include "Bond.h"

namespace RDKit {
  typedef boost::shared_ptr<Bond>  BOND_SPTR;
  struct edge_bond_t {
    enum { num = 1002 };
    typedef boost::edge_property_tag kind;
  };
  struct edge_wght_t {
    enum { num = 1003 };
    typedef boost::edge_property_tag kind;
  };

  typedef boost::property<edge_wght_t,double> BondWeight;
  typedef boost::property<edge_bond_t,Bond *,BondWeight> BondProperty;
}

#endif
