//
//  Copyright (C) 2011-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file USRDescriptor.h

  \brief Contains the USR descriptor. Use MolDescriptors.h in client code.

*/
#ifndef __RD_USR_H__
#define __RD_USR_H__

#include <Geometry/point.h>
#include <Numerics/Vector.h>

namespace RDKit{
  class ROMol;
  class Conformer;
  namespace Descriptors {
    /*!
      Calculates the ultra-fast shape recognition (USR)

      Reference: P. J. Ballester, W. G. Richards, JCC (2007), 28, 1711 - 1723.

      Derived from RDKit Python implementation of Jan Domanski
      who derived his code from Adrian Schreyer's code:
      http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd25370e7ef93301d0e28d49a0851c2/usrcat/geometry.py?at=default

      \param mol          the molecule of interest
      \param descriptor   storage for the computed USR descriptor
      \param confId       the conformer Id

    */
    void USR(const ROMol &mol, std::vector<double> &descriptor, int confId = -1);

    /*!
      Calculates the four distance distributions for the USR descriptor

      \param coords      the atom coordinates
      \param dist    storage for the four distance distributions

    */
    void calcUSRDistributions(const RDGeom::Point3DConstPtrVect &coords,
                              std::vector<std::vector<double> > &dist);

    /*!
      Calculates the four distance distributions for the USR descriptor

      \param coords      the atom coordinates
      \param points      vector with the points
      \param dist   storage for the distance distributions

    */
    void calcUSRDistributionsFromPoints(const RDGeom::Point3DConstPtrVect &coords,
                                        const std::vector<RDGeom::Point3D> &points,
                                        std::vector<std::vector<double> > &dist);

    /*!
      Calculates the USR descriptor from the four distance distributions

      \param dist        vector with the four distance distributions
      \param descriptor  storage for the computed USR descriptor

    */
    void calcUSRFromDistributions(const std::vector<std::vector<double> > &dist,
                                  std::vector<double> &descriptor);

    /*!
     Calculates the score between two USR descriptors

     \param d1    descriptor 1
     \param d2    descriptor 2

     \return the score
     */
    double calcUSRScore(const std::vector<double> &d1, const std::vector<double> &d2);

  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
