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
      Calculates the USR descriptor for a single conformer

      \param coords      the 3D coordinates of the atoms of a conformer
      \param descriptor  storage for the computed USR descriptor

    */
    void calcUSRForPoints(const RDGeom::Point3DConstPtrVect &coords, std::vector<double> &descriptor);

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
