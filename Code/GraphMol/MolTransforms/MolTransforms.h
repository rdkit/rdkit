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

#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>
#endif

namespace RDKit {
class ROMol;
class Atom;
class Conformer;
}

namespace RDGeom {
class Transform3D;
}

namespace MolTransforms {
void transformMolsAtoms(RDKit::ROMol *mol, RDGeom::Transform3D &tform);
void transformAtom(RDKit::Atom *atom, RDGeom::Transform3D &tform);

//! Compute the centroid of a conformer
/*!
  This is simple the average of the heavy atom locations in the conformer,
  not attention is paid to hydrogens or the differences in atomic radii

  \param conf     Conformer of interest
  \param ignoreHs If true, ignore hydrogen atoms
*/
RDGeom::Point3D computeCentroid(const RDKit::Conformer &conf,
                                bool ignoreHs = true);

#ifdef RDK_HAS_EIGEN3
//! Compute principal axes and moments for a conformer
/*!
  \param conf       Conformer of interest
  \param axes       used to return the principal axes
  \param  moments    used to return the principal moments
  \param ignoreHs   If true, ignore hydrogen atoms
  \param force      If true, the calculation will be carried out even if a cached value is present
  \param weights    If present used to weight the atomic coordinates

  \returns whether or not the calculation was successful
*/
bool computePrincipalAxesAndMoments(
    const RDKit::Conformer &conf,
    Eigen::Matrix3d &axes,
    Eigen::Vector3d &moments,
    bool ignoreHs = true,
    bool force = false,
    const std::vector<double> *weights = NULL);
#endif



//! Compute the transformation require to orient the conformation
//! along the principal axes about the center; i.e. center is made to coincide
// with the
//! origin, the largest princiapl axis with the x-axis, the next largest with
// the y-axis
//! and the smallest with the z-axis
/*!
  If center is not specified the the centroid of the conformer will be used
  \param conf                Conformer of interest
  \param center              Center to be used for canonicalization, defaults to
  the centroid of the
                             conformation
  \param normalizeCovar      Normalize the covariance matrix with the number of
  atoms
  \param ignoreHs            Optinally ignore hydrogens
*/
RDGeom::Transform3D *computeCanonicalTransform(
    const RDKit::Conformer &conf, const RDGeom::Point3D *center = 0,
    bool normalizeCovar = false, bool ignoreHs = true);

//! Transform the conformation using the specified transformation
void transformConformer(RDKit::Conformer &conf,
                        const RDGeom::Transform3D &trans);

//! Canonicalize the orientation of a conformer so that its principal axes
//! around the specified center point coincide with the x, y, z axes
/*!
  \param conf           The conformer of interest
  \param center         Optional center point about which the principal axes are
  computed
                        if not specified the centroid of the conformer will be
  used
  \param normalizeCovar Optionally normalize the covariance matrix by the number
  of atoms
  \param ignoreHs       If true, ignore hydrogen atoms

*/
void canonicalizeConformer(RDKit::Conformer &conf,
                           const RDGeom::Point3D *center = 0,
                           bool normalizeCovar = false, bool ignoreHs = true);

//! Canonicalize all the conformations in a molecule
/*!
  \param mol            the molecule of interest
  \param normalizeCovar Optionally normalize the covariance matrix by the number
  of atoms
  \param ignoreHs       If true, ignore hydrogens
*/
void canonicalizeMol(RDKit::ROMol &mol, bool normalizeCovar = false,
                     bool ignoreHs = true);

//! Get the bond length between the specified atoms i, j
double getBondLength(const RDKit::Conformer &conf, unsigned int iAtomId,
                     unsigned int jAtomId);

//! Set the bond length between the specified atoms i, j
//! (all atoms bonded to atom j are moved)
void setBondLength(RDKit::Conformer &conf, unsigned int iAtomId,
                   unsigned int jAtomId, double value);

//! Get the angle in radians among the specified atoms i, j, k
double getAngleRad(const RDKit::Conformer &conf, unsigned int iAtomId,
                   unsigned int jAtomId, unsigned int kAtomId);

//! Get the angle in degrees among the specified atoms i, j, k
inline double getAngleDeg(const RDKit::Conformer &conf, unsigned int iAtomId,
                          unsigned int jAtomId, unsigned int kAtomId) {
  return (180. / M_PI * getAngleRad(conf, iAtomId, jAtomId, kAtomId));
}

//! Set the angle in radians among the specified atoms i, j, k
//! (all atoms bonded to atom k are moved)
void setAngleRad(RDKit::Conformer &conf, unsigned int iAtomId,
                 unsigned int jAtomId, unsigned int kAtomId, double value);

//! Set the angle in degrees among the specified atoms i, j, k
//! (all atoms bonded to atom k are moved)
inline void setAngleDeg(RDKit::Conformer &conf, unsigned int iAtomId,
                        unsigned int jAtomId, unsigned int kAtomId,
                        double value) {
  setAngleRad(conf, iAtomId, jAtomId, kAtomId, value / 180. * M_PI);
}

//! Get the dihedral angle in radians among the specified atoms i, j, k, l
double getDihedralRad(const RDKit::Conformer &conf, unsigned int iAtomId,
                      unsigned int jAtomId, unsigned int kAtomId,
                      unsigned int lAtomId);

//! Get the dihedral angle in degrees among the specified atoms i, j, k, l
inline double getDihedralDeg(const RDKit::Conformer &conf, unsigned int iAtomId,
                             unsigned int jAtomId, unsigned int kAtomId,
                             unsigned int lAtomId) {
  return (180. / M_PI *
          getDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId));
}

//! Set the dihedral angle in radians among the specified atoms i, j, k, l
//! (all atoms bonded to atom l are moved)
void setDihedralRad(RDKit::Conformer &conf, unsigned int iAtomId,
                    unsigned int jAtomId, unsigned int kAtomId,
                    unsigned int lAtomId, double value);

//! Set the dihedral angle in degrees among the specified atoms i, j, k, l
//! (all atoms bonded to atom l are moved)
inline void setDihedralDeg(RDKit::Conformer &conf, unsigned int iAtomId,
                           unsigned int jAtomId, unsigned int kAtomId,
                           unsigned int lAtomId, double value) {
  setDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId, value / 180. * M_PI);
}
}
#endif
