//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_FORCEFIELD_H__
#define __RD_FORCEFIELD_H__

#include <vector>
#include <boost/smart_ptr.hpp>
#include <Geometry/point.h>
#include <GraphMol/Trajectory/Snapshot.h>

namespace RDKit {
namespace ForceFieldsHelper {
void RDKIT_FORCEFIELD_EXPORT normalizeAngleDeg(double &angleDeg);
void RDKIT_FORCEFIELD_EXPORT computeDihedral(
    const RDGeom::PointPtrVect &pos, unsigned int idx1, unsigned int idx2,
    unsigned int idx3, unsigned int idx4, double *dihedral = nullptr,
    double *cosPhi = nullptr, RDGeom::Point3D r[4] = nullptr,
    RDGeom::Point3D t[2] = nullptr, double d[2] = nullptr);
void RDKIT_FORCEFIELD_EXPORT computeDihedral(
    const double *pos, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, double *dihedral = nullptr, double *cosPhi = nullptr,
    RDGeom::Point3D r[4] = nullptr, RDGeom::Point3D t[2] = nullptr,
    double d[2] = nullptr);
void RDKIT_FORCEFIELD_EXPORT
computeDihedral(const RDGeom::Point3D *p1, const RDGeom::Point3D *p2,
                const RDGeom::Point3D *p3, const RDGeom::Point3D *p4,
                double *dihedral = nullptr, double *cosPhi = nullptr,
                RDGeom::Point3D r[4] = nullptr, RDGeom::Point3D t[2] = nullptr,
                double d[2] = nullptr);
}  // namespace ForceFieldsHelper
}  // namespace RDKit

namespace ForceFields {
class ForceFieldContrib;
typedef std::vector<int> INT_VECT;
typedef boost::shared_ptr<const ForceFieldContrib> ContribPtr;
typedef std::vector<ContribPtr> ContribPtrVect;

//-------------------------------------------------------
//! A class to store forcefields and handle minimization
/*!
   A force field is used like this (schematically):

   \verbatim
     ForceField ff;

     // add contributions:
     for contrib in contribs:
       ff.contribs().push_back(contrib);

     // set up the points:
     for positionPtr in positions:
       ff.positions().push_back(point);

     // initialize:
     ff.initialize()

     // and minimize:
     needsMore = ff.minimize();

   \endverbatim

   <b>Notes:</b>
     - The ForceField owns its contributions, which are stored using smart
       pointers.
     - Distance calculations are currently lazy; the full distance matrix is
       never generated.  In systems where the distance matrix is not sparse,
       this is almost certainly inefficient.

*/
class RDKIT_FORCEFIELD_EXPORT ForceField {
 public:
  //! construct with a dimension
  ForceField(unsigned int dimension = 3) : d_dimension(dimension) {}

  ~ForceField();

  //! copy ctor, copies contribs.
  ForceField(const ForceField &other);

  //! does initialization
  void initialize();

  //! calculates and returns the energy (in kcal/mol) based on existing
  /// positions in the forcefield
  /*!

  \return the current energy

    <b>Note:</b>
      This function is less efficient than calcEnergy with postions passed in as
  double *
      the positions need to be converted to double * here
  */
  double calcEnergy(std::vector<double> *contribs = nullptr) const;

  // these next two aren't const because they may update our
  // distance matrix

  //! calculates and returns the energy of the position passed in
  /*!
    \param pos an array of doubles.  Should be \c 3*this->numPoints() long.

    \return the current energy

    <b>Side effects:</b>
      - Calling this resets the current distance matrix
      - The individual contributions may further update the distance matrix
  */
  double calcEnergy(double *pos);

  //! calculates the gradient of the energy at the current position
  /*!

    \param forces an array of doubles.  Should be \c 3*this->numPoints() long.

    <b>Note:</b>
      This function is less efficient than calcGrad with positions passed in
      the positions need to be converted to double * here
   */
  void calcGrad(double *forces) const;

  //! calculates the gradient of the energy at the provided position
  /*!

    \param pos      an array of doubles.  Should be \c 3*this->numPoints() long.
    \param forces   an array of doubles.  Should be \c 3*this->numPoints() long.

    <b>Side effects:</b>
      - The individual contributions may modify the distance matrix
   */
  void calcGrad(double *pos, double *forces);

  //! minimizes the energy of the system by following gradients
  /*!
    \param maxIts    the maximum number of iterations to try
    \param forceTol  the convergence criterion for forces
    \param energyTol the convergence criterion for energies

    \return an integer value indicating whether or not the convergence
            criteria were achieved:
      - 0: indicates success
      - 1: the minimization did not converge in \c maxIts iterations.
  */
  int minimize(unsigned int snapshotFreq, RDKit::SnapshotVect *snapshotVect,
               unsigned int maxIts = 200, double forceTol = 1e-4,
               double energyTol = 1e-6);

  //! minimizes the energy of the system by following gradients
  /*!
    \param maxIts            the maximum number of iterations to try
    \param forceTol          the convergence criterion for forces
    \param energyTol         the convergence criterion for energies
    \param snapshotFreq      a snapshot of the minimization trajectory
                             will be stored after as many steps as indicated
                             through this parameter; defaults to 0 (no
                             trajectory stored)
    \param snapshotVect      a pointer to a std::vector<Snapshot> where
                             coordinates and energies will be stored

    \return an integer value indicating whether or not the convergence
            criteria were achieved:
      - 0: indicates success
      - 1: the minimization did not converge in \c maxIts iterations.
  */
  int minimize(unsigned int maxIts = 200, double forceTol = 1e-4,
               double energyTol = 1e-6);

  // ---------------------------
  // setters and getters

  //! returns a reference to our points (a PointPtrVect)
  RDGeom::PointPtrVect &positions() { return d_positions; }
  const RDGeom::PointPtrVect &positions() const { return d_positions; }

  //! returns a reference to our contribs (a ContribPtrVect)
  ContribPtrVect &contribs() { return d_contribs; }
  const ContribPtrVect &contribs() const { return d_contribs; }

  //! returns the distance between two points
  /*!
    \param i point index
    \param j point index
    \param pos (optional) If this argument is provided, it will be used
        to provide the positions of points. \c pos should be
        \c 3*this->numPoints() long.

    \return the distance

    <b>Side effects:</b>
      - if the distance between i and j has not previously been calculated,
        our internal distance matrix will be updated.
  */
  double distance(unsigned int i, unsigned int j, double *pos = nullptr);

  //! returns the distance between two points
  /*!
    \param i point index
    \param j point index
    \param pos (optional) If this argument is provided, it will be used
        to provide the positions of points. \c pos should be
        \c 3*this->numPoints() long.

    \return the distance

    <b>Note:</b>
      The internal distance matrix is not updated in this case
  */
  double distance(unsigned int i, unsigned int j, double *pos = nullptr) const;

  //! returns the squared distance between two points
  /*!
    \param i point index
    \param j point index
    \param pos (optional) If this argument is provided, it will be used
        to provide the positions of points. \c pos should be
        \c 3*this->numPoints() long.

    \return the squared distance

    <b>Note:</b>
      The internal distance matrix is not updated
  */
  double distance2(unsigned int i, unsigned int j, double *pos = nullptr) const;

  //! returns the dimension of the forcefield
  unsigned int dimension() const { return d_dimension; }

  //! returns the number of points the ForceField is handling
  unsigned int numPoints() const { return d_numPoints; }

  INT_VECT &fixedPoints() { return d_fixedPoints; }
  const INT_VECT &fixedPoints() const { return d_fixedPoints; }

 protected:
  unsigned int d_dimension;
  bool df_init{false};               //!< whether or not we've been initialized
  unsigned int d_numPoints{0};       //!< the number of active points
  double *dp_distMat{nullptr};       //!< our internal distance matrix
  RDGeom::PointPtrVect d_positions;  //!< pointers to the points we're using
  ContribPtrVect d_contribs;         //!< contributions to the energy
  INT_VECT d_fixedPoints;
  unsigned int d_matSize = 0;
  //! scatter our positions into an array
  /*!
      \param pos     should be \c 3*this->numPoints() long;
  */
  void scatter(double *pos) const;

  //! update our positions from an array
  /*!
      \param pos     should be \c 3*this->numPoints() long;
  */
  void gather(double *pos);

  //! initializes our internal distance matrix
  void initDistanceMatrix();
};
}  // namespace ForceFields
#endif
