//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef RDDEPICTOR_H
#define RDDEPICTOR_H

#include <RDGeneral/types.h>
#include <list>
#include <GraphMol/ROMol.h>


namespace RDKit {
  class ROMol;
}

namespace RDDepict {
  class EmbeddedFrag;

  //! \brief This the function the user gets to compute coodinates
  /*! 
    \param mol          the molecule were are interested in
    \param coordMap     a map of int to Point2D, between atom IDs and their locations. 
    This is the container the user need to fill if he/she wants to 
    specify coordinates for a portion of the molecule, defualts to 0
    \param canonOrient  canonicalize the orientation so that the the long axes 
    align with the x-axis etc.
    \param clearConfs   clear all existing conformations on the molecule them adding the
    2D coordinates or simple add to the list
    
                          
    \return ID of the conformation added to the molecule cotaining the 2D coordinates
  */
  unsigned int compute2DCoords(RDKit::ROMol &mol, const RDGeom::INT_POINT2D_MAP *coordMap=0,
                               bool canonOrient=false, bool clearConfs=true);

  namespace DepictorLocal {
    void embedFusedSystems(const RDKit::ROMol &mol, const RDKit::VECT_INT_VECT &arings,
                           std::list<EmbeddedFrag> &efrags);
    void embedCisTransSystems(const RDKit::ROMol &mol, 
                              std::list<EmbeddedFrag> &efrags);
    RDKit::INT_LIST getNonEmbeddedAtoms(const RDKit::ROMol &mol, const std::list<EmbeddedFrag> &efrags);
  };
};

#endif
