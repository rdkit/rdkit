// 
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_MOLTRANSFORMS_H_
#define _RD_MOLTRANSFORMS_H_

#include <Geometry/point.h>
#include <Numerics/SymmMatrix.h>

namespace RDKit{
  class ROMol;
  class Atom;
  class Conformer;
}

namespace RDGeom {
  class Transform3D;
}

namespace MolTransforms{
  void transformMolsAtoms(RDKit::ROMol *mol,RDGeom::Transform3D &tform);
  void transformAtom(RDKit::Atom *atom,RDGeom::Transform3D &tform);

  //! Compute the centroid of a conformer
  /*! 
    This is simple the average of the heavy atom locations in the conformer,
    not attention is paid to hydrogens or the differences in atomic radii
    
    \param conf     Conformer of interest
    \param ignoreHs If true, ignore hydrogen atoms
  */
  RDGeom::Point3D computeCentroid(const RDKit::Conformer &conf, bool ignoreHs=true);

  //! Compute the covariance matrix for a conformer
  /*!
    \param conf       Conformer of interest
    \param center     Center to be used for covariance matrix calculation
    \param normalize  If true, normalize the covariance matrix by the number of atoms
    \param ignoreHs   If true, ignore hydrogen atoms
  */
  RDNumeric::DoubleSymmMatrix *computeCovarianceMatrix(const RDKit::Conformer &conf,
                                                       const RDGeom::Point3D &center,
                                                       bool normalize=false, 
                                                       bool ignoreHs=true);
   

  //! Compute the transformation require to orient the conformation
  //! along the principal axes about the center; i.e. center is made to coincide with the
  //! origin, the largest princiapl axis with the x-axis, the next largest with the y-axis
  //! and the smallest with the z-axis
  /*!
    If center is not specified the the centroid of the conformer will be used
    \param conf                Conformer of interest
    \param center              Center to be used for canonicalization, defaults to the centroid of the 
                               conformation
    \param normalizeCovar      Normalize the covariance matrix with the number of atoms
    \param ignoreHs            Optinally ignore hydrogens
  */
  RDGeom::Transform3D *computeCanonicalTransform(const RDKit::Conformer &conf,
                                                 const RDGeom::Point3D *center=0,
                                                 bool normalizeCovar=false,
                                                 bool ignoreHs=true);

  //! Transform the conformation using the specified transformation
  void transformConformer(RDKit::Conformer &conf, const RDGeom::Transform3D &trans);

  //! Canonicalize the orientation of a conformer so that its principal axes
  //! around the specified center point coincide with the x, y, z axes
  /*!
    \param conf           The conformer of interest
    \param center         Optional center point about which the principal axes are computed 
                          if not specified the centroid of the conformer will be used
    \param normalizeCovar Optionally normalize the covariance matrix by the number of atoms
    \param ignoreHs       If true, ignore hydrogen atoms
    
  */
  void canonicalizeConformer(RDKit::Conformer &conf, const RDGeom::Point3D *center=0,
                             bool normalizeCovar=false, bool ignoreHs=true);

  //! Canonicalize all the conformations in a molecule
  /*!
    \param mol            the molecule of interest
    \param normalizeCovar Optionally normalize the covariance matrix by the number of atoms
    \param ignoreHs       If true, ignore hydrogens
  */
  void canonicalizeMol(RDKit::ROMol &mol, bool normalizeCovar=false, bool ignoreHs=true);

}
#endif
