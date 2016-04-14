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
        \param traj is the parent Trajectory
     */
    void setTrajectory(const Trajectory *traj);
    /*! \brief Gets the coordinates at pointNum as a Point2D object;
        requires the Trajectory dimension to be == 2
        \param pointNum is the atom number whose coordinates will be retrieved
     */
    Point2D getPoint2D(unsigned int pointNum) const;
    /*! \brief Gets the coordinates at pointNum as a Point3D object;
        requires the Trajectory dimension to be >= 2
        \param pointNum is the atom number whose coordinates will be retrieved
     */
    Point3D getPoint3D(unsigned int pointNum) const;
    /*! \brief Gets the energy for this snapshot
     */
    double getEnergy() const {
      return d_energy;
    };
    /*! \brief Sets the energy for this snapshot
        \param energy the energy value assigned to this snapshot
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
        \param data is a (void *), needs to be cast to the appropriate data type
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
      /*! if this flag is set, upon destruction of the Trajectory object
          or call of the removeSnapshot() method, the d_pos array
          of the relevant Snapshot object(s) will be freed
       */
      FREE_POS_ON_DESTROY = (1 << 0)
    };
    /*! \brief Constructor
        \param dimension represents the dimensionality of this Trajectory's coordinate tuples;
        this is normally 2 (2D coordinates) or 3 (3D coordinates)
        \param numPoints is the number of coordinate tuples associated to each Snapshot
        \param flags see the Flags enum
     */
    Trajectory(unsigned int dimension, unsigned int numPoints, unsigned int flags = 0);
    /*! \brief Destructor. If the FREE_POS_ON_DESTROY flag is set , upon destruction all d_pos
        arrays in this Trajectory's Snapshot objects are freed
     */
    ~Trajectory();
    /*! \brief Returns the dimensionality of this Trajectory's coordinate tuples
     */
    unsigned int dimension() const {
      return d_dimension;
    }
    /*! \brief Returns the number of coordinate tuples associated to each Snapshot
     */
    unsigned int numPoints() const {
      return d_numPoints;
    }
    /*! \brief Returns the number of Snapshots asociated to this Trajectory
     */
    unsigned int size() const {
      return d_snapshotVect.size();
    }
    /*! \brief Appends a Snapshot to this Trajectory
        \param s is the Snapshot to be added
     */
    unsigned int addSnapshot(Snapshot s);
    /*! \brief Retrieves a const reference to a Snapshot in the Trajectory
        \param snapshotNum is the zero-based index of the retrieved Snapshot
     */
    const Snapshot &getSnapshot(unsigned int snapshotNum) const;
    /*! \brief Inserts a Snapshot into this Trajectory
        \param snapshotNum is the zero-based index of the Trajectory's Snapshot
               before which the Snapshot s will be inserted
        \param s is the Snapshot to be inserted
     */
    unsigned int insertSnapshot(unsigned int snapshotNum, Snapshot s);
    /*! \brief Removes a Snapshot from this Trajectory
        \param snapshotNum is the zero-based index of Snapshot to be removed
     */
    unsigned int removeSnapshot(unsigned int snapshotNum);
    /*! \brief Retrieves the status of the FREE_POS_ON_DESTROY flag
        true: set, false: not set
     */
    bool getFreePosOnDestroy() const {
      return ((d_flags & FREE_POS_ON_DESTROY) ? true : false);
    }
    /*! \brief Sets the status of the FREE_POS_ON_DESTROY flag
        true: set, false: not set
     */
    void setFreePosOnDestroy(bool free) {
      if (free)
        d_flags |= FREE_POS_ON_DESTROY;
      else
        d_flags &= ~FREE_POS_ON_DESTROY;
    }
  private:
    // dimensionality of this Trajectory's coordinates;
    // this is normally 2 (2D coordinates) or 3 (3D coordinates)
    const unsigned int d_dimension;
    // number of coordinate tuples associated to each Snapshot
    const unsigned int d_numPoints;
    // flags
    unsigned int d_flags;
    // vector holding the Snapshots for this Trajectory
    std::vector<Snapshot> d_snapshotVect;
};

}
#endif
