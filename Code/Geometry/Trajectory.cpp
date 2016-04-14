//
// Copyright (C) 2003-2016 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include "Trajectory.h"

namespace RDGeom {

Snapshot::Snapshot(double *pos, double energy, void *data) :
  d_trajectory(NULL),
  d_energy(energy),
  d_pos(pos),
  d_data(data)
{
}

Point2D Snapshot::getPoint2D(unsigned int pointNum) const {
  PRECONDITION(d_pos, "pos must not be NULL");
  PRECONDITION(trajectory()->dimension() == 2, "d_dimension must be == 2");
  unsigned int i = pointNum * trajectory()->dimension();
  return Point2D(d_pos[i], d_pos[i + 1]);
}

Point3D Snapshot::getPoint3D(unsigned int pointNum) const {
  PRECONDITION(trajectory()->dimension() >= 2, "d_dimension must be >= 2");
  unsigned int i = pointNum * trajectory()->dimension();
  return (d_pos ? Point3D(d_pos[i], d_pos[i + 1],
          ((trajectory()->dimension() == 3) ? d_pos[i + 2] : 0.0)) : Point3D());
}

void Snapshot::setTrajectory(const Trajectory *traj) {
  PRECONDITION(traj, "traj must not be NULL");
  d_trajectory = traj;
}

void Snapshot::freePos() {
  delete d_pos;
  d_pos = NULL;
}

Trajectory::Trajectory(unsigned int dimension, unsigned int numPoints, unsigned int flags) :
  d_dimension(dimension),
  d_numPoints(numPoints),
  d_flags(flags)
{
}

Trajectory::~Trajectory() {
  for (std::vector<Snapshot>::iterator it = d_snapshotVect.begin();
    (d_flags & FREE_POS_ON_DESTROY) && (it != d_snapshotVect.end()); ++it)
    it->freePos();
}

unsigned int Trajectory::addSnapshot(Snapshot s) {
  unsigned int size = d_snapshotVect.size();
  s.setTrajectory(this);
  d_snapshotVect.push_back(s);
  return size;
}

const Snapshot &Trajectory::getSnapshot(unsigned int snapshotNum) const {
  PRECONDITION(snapshotNum < d_snapshotVect.size(), "snapshotNum out of bounds");
  return d_snapshotVect[snapshotNum];
}

unsigned int Trajectory::insertSnapshot(unsigned int snapshotNum, Snapshot s) {
  PRECONDITION(snapshotNum < d_snapshotVect.size(), "snapshotNum out of bounds");
  s.setTrajectory(this);
  return (d_snapshotVect.insert(d_snapshotVect.begin() + snapshotNum, s)
          - d_snapshotVect.begin());
}

unsigned int Trajectory::removeSnapshot(unsigned int snapshotNum) {
  PRECONDITION(snapshotNum < d_snapshotVect.size(), "snapshotNum out of bounds");
  if (d_flags & FREE_POS_ON_DESTROY)
    d_snapshotVect[snapshotNum].freePos();
  return (d_snapshotVect.erase(d_snapshotVect.begin() + snapshotNum) - d_snapshotVect.begin());
}

}
