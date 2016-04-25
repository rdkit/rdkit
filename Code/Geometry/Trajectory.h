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
#include <map>
#include "point.h"
#include "boost/shared_array.hpp"

namespace RDGeom {

class Trajectory;

class Snapshot {
  friend class Trajectory;
  public:
    /*! \brief Constructor
        \param pos is a pointer to an array of (numPoints * dimension) doubles;
        numPoints and dimension must match the Trajectory which is going to
        contain this Snapshot
        \param energy is the energy associated with this set of coordinates
     */
    Snapshot(boost::shared_array<double> pos = boost::shared_array<double>(), double energy = 0.0);
    /*! \return a const pointer to the parent Trajectory
     */
    const Trajectory *trajectory() const {
      return d_trajectory;
    }
    /*! \param pointNum is the atom number whose coordinates will be retrieved
        \return the coordinates at pointNum as a Point2D object;
        requires the Trajectory dimension to be == 2
     */
    Point2D getPoint2D(unsigned int pointNum) const;
    /*! \param pointNum is the atom number whose coordinates will be retrieved
        \return the coordinates at pointNum as a Point3D object;
        requires the Trajectory dimension to be >= 2
     */
    Point3D getPoint3D(unsigned int pointNum) const;
    /*! \return the energy for this Snapshot
     */
    double getEnergy() const {
      return d_energy;
    };
    /*! \brief Sets the energy for this Snapshot
        \param energy the energy value assigned to this Snapshot
     */
    void setEnergy(double energy) {
      d_energy = energy;
    }
    /*! \brief Frees the pointer to the array of doubles where the
        coordinates for this Snapshot are stored
     */
  private:
    // Pointer to the parent Trajectory object
    const Trajectory *d_trajectory;
    // Energy for this set of coordinates
    double d_energy;
    // shared array to Snapshot coordinates
    boost::shared_array<double> d_pos;
};

class Trajectory {
  public:
    typedef std::vector<Snapshot *> SnapshotPtrVect;
    /*! \brief Constructor
        \param dimension represents the dimensionality of this Trajectory's coordinate tuples;
        this is normally 2 (2D coordinates) or 3 (3D coordinates)
        \param numPoints is the number of coordinate tuples associated to each Snapshot
     */
    Trajectory(unsigned int dimension, unsigned int numPoints);
    /*! \brief Copy constructor
     */
    Trajectory(const Trajectory &other);
    /*! \brief Destructor
        as the Trajectory is the owner of the Snapshot, all Snapshots are destroyed
        upon destruction of the Trajectory
     */
    ~Trajectory();
    /*! \return the dimensionality of this Trajectory's coordinate tuples
     */
    unsigned int dimension() const {
      return d_dimension;
    }
    /*! \return the number of coordinate tuples associated to each Snapshot
     */
    unsigned int numPoints() const {
      return d_numPoints;
    }
    /*! \return the number of Snapshots associated to this Trajectory
     */
    size_t size() const {
      return d_snapshotVect.size();
    }
    /*! \brief Appends a Snapshot to this Trajectory
        \param s is the Snapshot to be added; the Trajectory
               takes ownership of the snapshot coordinates
        \return the zero-based index position of the added Snapshot
     */
    unsigned int addSnapshot(Snapshot *s);
    /*! \param snapshotNum is the zero-based index of the retrieved Snapshot
        \return a const reference to the relevant Snapshot in the Trajectory
     */
    Snapshot *getSnapshot(unsigned int snapshotNum) const;
    /*! \brief Inserts a Snapshot into this Trajectory
        \param snapshotNum is the zero-based index of the Trajectory's Snapshot
               before which the Snapshot s will be inserted
        \param s is the Snapshot to be inserted; the Trajectory takes ownership of the Snapshot
        \return the zero-based index position of the inserted Snapshot
     */
    unsigned int insertSnapshot(unsigned int snapshotNum, Snapshot *s);
    /*! \brief Removes a Snapshot from this Trajectory
        \param snapshotNum is the zero-based index of Snapshot to be removed
        \return the zero-based index position of the Snapshot after the
                removed one; if the last snapsot was removed, it returns the
                size of the trajectory
     */
    unsigned int removeSnapshot(unsigned int snapshotNum);
    /*! \brief Reads coordinates from an AMBER trajectory file
               into the Trajectory object
        \return the number of Snapshot objects read in
     */
    unsigned int readAmber(const std::string &fName);
    /*! \brief Reads coordinates from a GROMOS trajectory file
               into the Trajectory object
        \return the number of Snapshot objects read in
     */
    unsigned int readGromos(const std::string &fName);
  private:
    // dimensionality of this Trajectory's coordinates;
    // this is normally 2 (2D coordinates) or 3 (3D coordinates)
    const unsigned int d_dimension;
    // number of coordinate tuples associated to each Snapshot
    const unsigned int d_numPoints;
    // vector holding the Snapshots for this Trajectory
    SnapshotPtrVect d_snapshotVect;
};

}
#endif
