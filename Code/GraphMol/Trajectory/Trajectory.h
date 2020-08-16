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
#ifndef RD_TRAJECTORY_H
#define RD_TRAJECTORY_H
#include <vector>
#include "Snapshot.h"

namespace RDKit {

class ROMol;

class RDKIT_TRAJECTORY_EXPORT Trajectory {
 public:
  /*! \brief Constructor
      \param dimension represents the dimensionality of this Trajectory's
     coordinate tuples; this is normally 2 (2D coordinates) or 3 (3D
     coordinates) \param numPoints is the number of coordinate tuples associated
     to each Snapshot \param snapshotVect (optional, defaults to NULL) is a
     pointer to a SnapshotVect used to initialize the Trajectory; if not NULL,
     the Trajectory takes ownership of the SnapshotVect
   */
  Trajectory(unsigned int dimension, unsigned int numPoints,
             SnapshotVect *snapshotVect = nullptr);
  /*! \brief Copy constructor
   */
  Trajectory(const Trajectory &other);
  /*! \return the dimensionality of this Trajectory's coordinate tuples
   */
  unsigned int dimension() const { return d_dimension; }
  /*! \return the number of coordinate tuples associated to each Snapshot
   */
  unsigned int numPoints() const { return d_numPoints; }
  /*! \return the number of Snapshots associated to this Trajectory
   */
  size_t size() const { return d_snapshotVect->size(); }
  /*! \brief Appends a Snapshot to this Trajectory
      \param s is the Snapshot to be added; the Trajectory
             takes ownership of the snapshot coordinates
      \return the zero-based index position of the added Snapshot
   */
  unsigned int addSnapshot(const Snapshot &s);
  /*! \param snapshotNum is the zero-based index of the retrieved Snapshot
      \return a const reference to the relevant Snapshot in the Trajectory
   */
  const Snapshot &getSnapshot(unsigned int snapshotNum) const;
  /*! \brief Inserts a Snapshot into this Trajectory
      \param snapshotNum is the zero-based index of the Trajectory's Snapshot
             before which the Snapshot s will be inserted
      \param s is the Snapshot to be inserted; the Trajectory
             takes ownership of the snapshot coordinates
      \return the zero-based index position of the inserted Snapshot
   */
  unsigned int insertSnapshot(unsigned int snapshotNum, Snapshot s);
  /*! \brief Removes a Snapshot from this Trajectory
      \param snapshotNum is the zero-based index of Snapshot to be removed
      \return the zero-based index position of the Snapshot after the
              removed one; if the last Snapshot was removed, it returns the
              size of the trajectory
   */
  unsigned int removeSnapshot(unsigned int snapshotNum);
  //! Clear all Snapshots from a Trajectory
  void clear() { d_snapshotVect->clear(); };
  //! Add conformations from the Trajectory to a molecule
  /*!
    \param mol - ROMol to which Conformers with coordinates from the Trajectory
    will be added; the Trajectory must have numPoints() == mol.getNumAtoms()
    \param from - the first Snapshot that will be added as a Conformer; defaults
    to -1 (first available) \param to - the last Snapshot that will be added as
    a Conformer; defaults to -1 (all) \return the number of conformations added
  */
  unsigned int addConformersToMol(ROMol &mol, int from = -1, int to = -1);

 private:
  // dimensionality of this Trajectory's coordinates;
  // this is normally 2 (2D coordinates) or 3 (3D coordinates)
  const unsigned int d_dimension;
  // number of coordinate tuples associated to each Snapshot
  const unsigned int d_numPoints;
  // smart_ptr to vector holding the Snapshots for this Trajectory
  boost::shared_ptr<SnapshotVect> d_snapshotVect;
};
/*! \brief Reads coordinates from an AMBER trajectory file
           into the traj Trajectory object
    \return the number of Snapshot objects read in
 */
RDKIT_TRAJECTORY_EXPORT unsigned int readAmberTrajectory(
    const std::string &fName, Trajectory &traj);
/*! \brief Reads coordinates from a GROMOS trajectory file
           into the traj Trajectory object
    \return the number of Snapshot objects read in
 */
RDKIT_TRAJECTORY_EXPORT unsigned int readGromosTrajectory(
    const std::string &fName, Trajectory &traj);

}  // namespace RDKit
#endif
