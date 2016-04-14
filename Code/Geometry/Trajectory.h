//
// Copyright (C) 2003-2016 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef __RD_TRAJECTORY_H__
#define __RD_TRAJECTORY_H__
#include <vector>
#include "point.h"

namespace RDGeom {

class Trajectory;

class Snapshot {
  public:
    /*! \brief Constructor
        \param pos is a pointer to an array of (numPoints * dimension) doubles;
        numPoints and dimension must match the Trajectory which is going to
        contain this Snapshot
        If the FREE_POS_ON_DESTROY flag is set on the Trajectory object which
        contains this Snapshot, this pointer will be freed upon:
        1) destruction of the Trajectory object which contains the Snapshot
        2) call of the removeSnapshot() method on the Trajectory object
           which contains the Snapshot
        \param energy is the energy associated with this set of coordinates
        \param data is a pointer to user data
     */
    Snapshot(double *pos = NULL, double energy = 0.0, void *data = NULL);
    /*! \brief Returns a pointer to the parent Trajectory
     */
    const Trajectory *trajectory() const {
      return d_trajectory;
    }
    /*! \brief Sets the pointer to the parent Trajectory
     */
    void setTrajectory(const Trajectory *traj);
    /*! \brief Gets the coordinates at pointNum as a Point2D object;
        requires the Trajectory dimension to be == 2
     */
    Point2D getPoint2D(unsigned int pointNum) const;
    /*! \brief Gets the coordinates at pointNum as a Point3D object;
        requires the Trajectory dimension to be >= 2
     */
    Point3D getPoint3D(unsigned int pointNum) const;
    /*! \brief Gets the energy for this snapshot
     */
    double getEnergy() const {
      return d_energy;
    };
    /*! \brief Sets the energy for this snapshot
     */
    void setEnergy(double energy) {
      d_energy = energy;
    }
    /*! \brief Returns the pointer to user data
     */
    void *getData() const {
      return d_data;
    }
    /*! \brief Sets the pointer to user data
     */
    void setData(void *data) {
      d_data = data;
    }
    /*! \brief Frees the pointer to the array of doubles where the
        coordinates for this snapshot are stored
     */
    void freePos();
  private:
    // Pointer to the parent Trajectory object
    const Trajectory *d_trajectory;
    // Energy for this set of coordinates
    double d_energy;
    // Pointer to snapshot coordinates. If the FREE_POS_ON_DESTROY flag is set
    // on the Trajectory object which contains this Snapshot, this pointer
    // will be freed upon:
    // 1) destruction of the Trajectory object which contains the Snapshot
    // 2) call of the removeSnapshot() method on the Trajectory object
    //    which contains the Snapshot
    double *d_pos;
    // Pointer to user data. Please note that this pointer
    // will not be freed upon object destruction
    void *d_data;
};

class Trajectory {
  public:
    enum Flags {
      FREE_POS_ON_DESTROY = (1 << 0)
    };
    Trajectory(unsigned int dimension, unsigned int numPoints, unsigned int flags = 0);
    ~Trajectory();
    unsigned int dimension() const {
      return d_dimension;
    }
    unsigned int numPoints() const {
      return d_numPoints;
    }
    unsigned int size() const {
      return d_snapshotVect.size();
    }
    unsigned int addSnapshot(Snapshot s);
    const Snapshot &getSnapshot(unsigned int snapshotNum) const;
    unsigned int insertSnapshot(unsigned int snapshotNum, Snapshot s);
    unsigned int removeSnapshot(unsigned int snapshotNum);
    bool getFreePosOnDestroy() const {
      return (d_flags & FREE_POS_ON_DESTROY);
    }
    void setFreePosOnDestroy(bool free) {
      if (free)
        d_flags |= FREE_POS_ON_DESTROY;
      else
        d_flags &= ~FREE_POS_ON_DESTROY;
    }
  private:
    const unsigned int d_dimension;
    const unsigned int d_numPoints;
    unsigned int d_flags;
    std::vector<Snapshot> d_snapshotVect;
};

}
#endif
