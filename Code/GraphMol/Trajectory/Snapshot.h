//
// Copyright (C) 2003-2016 Sereina Riniker, Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef __RD_SNAPSHOT_H__
#define __RD_SNAPSHOT_H__
#include <Geometry/point.h>
#include <boost/shared_array.hpp>

namespace RDKit {
class Snapshot;
class Trajectory;
typedef std::vector<Snapshot> SnapshotVect;
}  // namespace RDKit

namespace RDKit {

class RDKIT_TRAJECTORY_EXPORT Snapshot {
  friend class Trajectory;

 public:
  /*! \brief Constructor
      \param pos is a pointer to an array of (numPoints * dimension) doubles;
      numPoints and dimension must match the Trajectory which is going to
      contain this Snapshot
      \param energy is the energy associated with this set of coordinates
   */
  Snapshot(boost::shared_array<double> pos, double energy = 0.0)
      : d_trajectory(nullptr), d_energy(energy), d_pos(pos) {}
  /*! \return a const pointer to the parent Trajectory
   */
  const Trajectory *trajectory() const { return d_trajectory; }
  /*! \param pointNum is the atom number whose coordinates will be retrieved
      \return the coordinates at pointNum as a Point2D object;
      requires the Trajectory dimension to be == 2
   */
  RDGeom::Point2D getPoint2D(unsigned int pointNum) const;
  /*! \param pointNum is the atom number whose coordinates will be retrieved
      \return the coordinates at pointNum as a Point3D object;
      requires the Trajectory dimension to be >= 2
   */
  RDGeom::Point3D getPoint3D(unsigned int pointNum) const;
  /*! \return the energy for this Snapshot
   */
  double getEnergy() const { return d_energy; };
  /*! \brief Sets the energy for this Snapshot
      \param energy the energy value assigned to this Snapshot
   */
  void setEnergy(double energy) { d_energy = energy; }

 private:
  // Pointer to the parent Trajectory object
  const Trajectory *d_trajectory;
  // Energy for this set of coordinates
  double d_energy;
  // shared array to Snapshot coordinates
  boost::shared_array<double> d_pos;
};

}  // namespace RDKit
#endif
