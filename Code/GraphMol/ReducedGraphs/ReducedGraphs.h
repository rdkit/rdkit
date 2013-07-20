//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_REDUCEDGRAPHS_H_
#define _RD_REDUCEDGRAPHS_H_

#include <vector>
#include <boost/cstdint.hpp>
#include <boost/dynamic_bitset.hpp>

namespace RDKit{
  class ROMol;

  //! \brief Generates a reduced graph representation of a molecule
  /*!

    \param mol:          the molecule to be fingerprinted

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
    
  */
  ROMol *createMolExtendedReducedGraph(const ROMol &mol,
                                       std::vector<boost::dynamic_bitset<> > *atomTypes=0
                                                   );
}

#endif
