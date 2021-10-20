//
// Copyright (C) 2003-2016 Sereina Riniker, Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <fstream>
#include <sstream>
#include <set>
#include "Trajectory.h"

namespace RDKit {

RDGeom::Point2D Snapshot::getPoint2D(unsigned int pointNum) const {
  PRECONDITION(d_pos, "d_pos must not be NULL");
  PRECONDITION(d_trajectory, "d_trajectory must not be NULL");
  PRECONDITION(d_trajectory->dimension() == 2, "d_dimension must be == 2");
  PRECONDITION(d_trajectory->numPoints(), "d_numPoints must be > 0");
  URANGE_CHECK(pointNum, d_trajectory->numPoints());
  unsigned int i = pointNum * d_trajectory->dimension();
  return RDGeom::Point2D(d_pos[i], d_pos[i + 1]);
}

RDGeom::Point3D Snapshot::getPoint3D(unsigned int pointNum) const {
  PRECONDITION(d_pos, "d_pos must not be NULL");
  PRECONDITION(d_trajectory, "d_trajectory must not be NULL");
  PRECONDITION(d_trajectory->dimension() >= 2, "d_dimension must be >= 2");
  PRECONDITION(d_trajectory->numPoints(), "d_numPoints must be > 0");
  URANGE_CHECK(pointNum, d_trajectory->numPoints());
  unsigned int i = pointNum * d_trajectory->dimension();
  return (
      RDGeom::Point3D(d_pos[i], d_pos[i + 1],
                      (d_trajectory->dimension() == 3) ? d_pos[i + 2] : 0.0));
}

Trajectory::Trajectory(unsigned int dimension, unsigned int numPoints,
                       SnapshotVect *snapshotVect)
    : d_dimension(dimension), d_numPoints(numPoints) {
  if (!snapshotVect) {
    snapshotVect = new SnapshotVect;
  }
  d_snapshotVect.reset(snapshotVect);
  for (auto &vectIt : *d_snapshotVect) {
    vectIt.d_trajectory = this;
  }
}

Trajectory::Trajectory(const Trajectory &other)
    : d_dimension(other.d_dimension),
      d_numPoints(other.d_numPoints),
      d_snapshotVect(new SnapshotVect) {
  for (SnapshotVect::const_iterator vectIt = other.d_snapshotVect->begin();
       vectIt != other.d_snapshotVect->end(); ++vectIt) {
    addSnapshot(*vectIt);
  }
}

unsigned int Trajectory::addSnapshot(const Snapshot &s) {
  return insertSnapshot(d_snapshotVect->size(), s);
}

const Snapshot &Trajectory::getSnapshot(unsigned int snapshotNum) const {
  URANGE_CHECK(snapshotNum, d_snapshotVect->size());
  return (*d_snapshotVect)[snapshotNum];
}

unsigned int Trajectory::insertSnapshot(unsigned int snapshotNum, Snapshot s) {
  URANGE_CHECK(snapshotNum, d_snapshotVect->size() + 1);
  s.d_trajectory = this;
  return (d_snapshotVect->insert(d_snapshotVect->begin() + snapshotNum, s) -
          d_snapshotVect->begin());
}

unsigned int Trajectory::removeSnapshot(unsigned int snapshotNum) {
  URANGE_CHECK(snapshotNum, d_snapshotVect->size());
  return (d_snapshotVect->erase(d_snapshotVect->begin() + snapshotNum) -
          d_snapshotVect->begin());
}

unsigned int Trajectory::addConformersToMol(ROMol &mol, int from, int to) {
  PRECONDITION(d_numPoints == mol.getNumAtoms(),
               "Number of atom mismatch between ROMol and Trajectory");
  PRECONDITION(from < static_cast<int>(size()), "from must be < size()");
  PRECONDITION(to < static_cast<int>(size()), "to must be < size()");
  if (from < 0) {
    from = 0;
  }
  if (to < 0) {
    to = size() - 1;
  }
  PRECONDITION(!size() || (from <= to), "from must be <= to");
  int n;
  unsigned int nConf;
  for (n = from, nConf = 0; size() && (n <= to); ++n, ++nConf) {
    auto *conf = new Conformer(mol.getNumAtoms());
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      conf->setAtomPos(i, getSnapshot(n).getPoint3D(i));
    }
    mol.addConformer(conf, true);
  }
  return nConf;
}

unsigned int readAmberTrajectory(const std::string &fName, Trajectory &traj) {
  PRECONDITION(traj.dimension() == 3,
               "The trajectory must have dimension == 3");
  std::ifstream inStream(fName.c_str());
  if (!inStream || inStream.bad()) {
    std::stringstream ss;
    ss << "Bad input file: " << fName;
    throw RDKit::BadFileException(ss.str());
  }
  std::string tempStr;
  // title
  std::getline(inStream, tempStr);
  // read coordinates
  unsigned int nCoords = traj.numPoints() * 3;
  unsigned int nSnapshots = 0;
  while (inStream.good() && !inStream.eof()) {
    boost::shared_array<double> c(new double[nCoords]());
    unsigned int i = 0;
    while (i < nCoords) {
      if (!(inStream >> c[i])) {
        if (!inStream.eof()) {
          std::stringstream ss;
          ss << "Error while reading file: " << fName;
          throw ValueErrorException(ss.str());
        } else if (i && (i < (nCoords - 1))) {
          std::stringstream ss;
          ss << "Premature end of file: " << fName;
          throw ValueErrorException(ss.str());
        }
        break;
      }
      ++i;
    }
    if (!inStream.eof()) {
      traj.addSnapshot(Snapshot(c));
      ++nSnapshots;
    }
  }
  return nSnapshots;
}

unsigned int readGromosTrajectory(const std::string &fName, Trajectory &traj) {
  PRECONDITION(traj.dimension() == 3,
               "The trajectory must have dimension == 3");
  std::ifstream inStream(fName.c_str());
  if (!inStream || inStream.bad()) {
    std::stringstream ss;
    ss << "Bad input file: " << fName;
    throw RDKit::BadFileException(ss.str());
  }
  std::string tempStr;
  unsigned int nCoords = traj.numPoints() * 3;
  unsigned int nSnapshots = 0;
  const static char *ignoredKeywordArray[] = {
      "TITLE", "TIMESTEP", "VELOCITYRED", "VELOCITY", "GENBOX", "BOX"};
  std::set<std::string> ignoredKeywordSet;
  for (auto &kw : ignoredKeywordArray) {
    ignoredKeywordSet.insert(std::string(kw));
  }
  while (inStream.good() && !inStream.eof()) {
    std::getline(inStream, tempStr);
    if (inStream.bad() || inStream.eof()) {
      continue;
    }
    if (ignoredKeywordSet.find(tempStr) != ignoredKeywordSet.end()) {
      // ignored block
      while (inStream.good() && !inStream.eof() && (tempStr != "END")) {
        std::getline(inStream, tempStr);
      }
    } else if ((tempStr == "POSITIONRED") || (tempStr == "POSITION")) {
      // these are the positions
      boost::shared_array<double> c(new double[nCoords]());
      unsigned int j = 0;
      for (unsigned int i = 0; i < traj.numPoints();) {
        std::getline(inStream, tempStr);
        if (inStream.bad() || inStream.eof() || (tempStr == "END")) {
          throw ValueErrorException("Wrong number of coordinates");
        }
        // ignore comments
        if (tempStr.find("#") != std::string::npos) {
          continue;
        }
        std::stringstream ls(tempStr);
        double x, y, z;
        if (!(ls >> x >> y >> z)) {
          throw ValueErrorException("Error while reading file");
        }
        // store the coordinates (convert to Angstrom!)
        c[j++] = x * 10.0;
        c[j++] = y * 10.0;
        c[j++] = z * 10.0;
        ++i;
      }
      std::getline(inStream, tempStr);  // the END line
      if (inStream.bad() || inStream.eof() || (tempStr != "END")) {
        throw ValueErrorException("Wrong number of coordinates");
      }
      traj.addSnapshot(Snapshot(c));
      ++nSnapshots;
    } else {
      std::string supportedBlocks("POSITIONRED, POSITION");
      for (const auto &it : ignoredKeywordSet) {
        supportedBlocks += ", " + it;
      }
      throw ValueErrorException("Unsupported block: " + tempStr +
                                ". Supported blocks are " + supportedBlocks);
    }
  }  // read file
  if (inStream.bad()) {
    std::stringstream ss;
    ss << "Bad input file: " << fName;
    throw RDKit::BadFileException(ss.str());
  }
  return nSnapshots;
}
}  // namespace RDKit
