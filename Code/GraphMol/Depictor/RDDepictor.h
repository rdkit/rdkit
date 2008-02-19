//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef _RDDEPICTOR_H_
#define _RDDEPICTOR_H_

#include <RDGeneral/types.h>
#include <list>
#include <GraphMol/ROMol.h>
#include "EmbeddedFrag.h"

namespace RDKit {
  class ROMol;
}

namespace RDDepict {

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
                               bool canonOrient=false, bool clearConfs=true,
                               unsigned int nFlipsPerSample=0, unsigned int nSamples=0,
                               int sampleSeed=0, bool permuteDeg4Nodes=false);

  //! \brief Compute the 2D coordinates such the interatom distances mimic those in ditance matrix
  /*!
    This function generates 2D coordinates such that the inter atom distance mimic those specified via
    dmat. This is done by randomly sampling(flipping) the rotatable bonds in the molecule and evaluating 
    a cost function which contains two components. The first component is the sum of inverse of the squared 
    inter-atom distances, this helps in spreading the atoms far from each other. The second component is
    the sum of squares of the difference in distance between those in dmat and the generated structure. 
    The user can adjust the relative importance of the two components via a adjustable paramter (see below)

    ARGUMENTS:
    \param mol - molecule involved in the frgament
    \param dmat - the distance matrix we want to mimic, this is symmteric N by N matrix when N is the number
                  of atoms in mol. All ngative entries in dmat are ignored.
    \param canonOrient - canonicalze the orientation after the 2D embedding is done
    \param clearConfs - clear any previously existing conformations on mol before adding a conformation 
    \param weightDistMat - A value between 0.0 and 1.0, this determines the importance of mimicing the 
                           the inter atoms distances in dmat. (1.0 - weightDistMat) is the weight associated
                           to spreading out the structure (density) in the cost function
    \param nFlipsPerSample - the number of rotatable bonds that are flips at random for each sample
    \param nSample - the number of samples
    \param sampleSeed - seed for the random sampling process
  */
  unsigned int compute2DCoordsMimicDistMat(RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat=0,
                                           bool canonOrient=true, bool clearConfs=true, double weightDistMat=0.5,
                                           unsigned int nFlipsPerSample=3, unsigned int nSamples=100,
                                           int sampleSeed=25, bool permuteDeg4Nodes=true);
  namespace DepictorLocal {
    void embedFusedSystems(const RDKit::ROMol &mol, const RDKit::VECT_INT_VECT &arings,
                           std::list<EmbeddedFrag> &efrags);
    void embedCisTransSystems(const RDKit::ROMol &mol, 
                              std::list<EmbeddedFrag> &efrags);
    RDKit::INT_LIST getNonEmbeddedAtoms(const RDKit::ROMol &mol, const std::list<EmbeddedFrag> &efrags);
  };
};

#endif
