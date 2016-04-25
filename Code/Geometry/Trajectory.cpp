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
#include <RDGeneral/BadFileException.h>
#include <fstream>
#include <sstream>
#include <set>
#include "Trajectory.h"

namespace RDGeom {

Snapshot::Snapshot(boost::shared_array<double> pos, double energy) :
  d_trajectory(NULL),
  d_energy(energy),
  d_pos(pos) {}

Point2D Snapshot::getPoint2D(unsigned int pointNum) const {
  PRECONDITION(d_pos, "pos must not be NULL");
  PRECONDITION(trajectory()->dimension() == 2, "d_dimension must be == 2");
  unsigned int i = pointNum * trajectory()->dimension();
  return Point2D(d_pos[i], d_pos[i + 1]);
}

Point3D Snapshot::getPoint3D(unsigned int pointNum) const {
  PRECONDITION(d_pos, "pos must not be NULL");
  PRECONDITION(trajectory()->dimension() >= 2, "d_dimension must be >= 2");
  unsigned int i = pointNum * trajectory()->dimension();
  return (Point3D(d_pos[i], d_pos[i + 1],
          (trajectory()->dimension() == 3) ? d_pos[i + 2] : 0.0));
}

Trajectory::Trajectory(unsigned int dimension, unsigned int numPoints) :
  d_dimension(dimension),
  d_numPoints(numPoints) {}

Trajectory::Trajectory(const Trajectory &other) :
  d_dimension(other.d_dimension),
  d_numPoints(other.d_numPoints) {
  for (SnapshotPtrVect::const_iterator vectIt = other.d_snapshotVect.begin();
    vectIt != other.d_snapshotVect.end(); ++vectIt) {
    Snapshot *s = new Snapshot(*(*vectIt));
    addSnapshot(s);
  }
}

Trajectory::~Trajectory()
{
  for (SnapshotPtrVect::const_iterator vectIt = d_snapshotVect.begin();
    vectIt != d_snapshotVect.end(); ++ vectIt)
    delete *vectIt;
}

unsigned int Trajectory::addSnapshot(Snapshot *s) {
  return insertSnapshot(d_snapshotVect.size(), s);
}

Snapshot *Trajectory::getSnapshot(unsigned int snapshotNum) const {
  PRECONDITION(snapshotNum < d_snapshotVect.size(), "snapshotNum out of bounds");
  return d_snapshotVect[snapshotNum];
}

unsigned int Trajectory::insertSnapshot(unsigned int snapshotNum, Snapshot *s) {
  PRECONDITION(snapshotNum <= d_snapshotVect.size(), "snapshotNum out of bounds");
  s->d_trajectory = this;
  return (d_snapshotVect.insert(d_snapshotVect.begin() + snapshotNum, s)
          - d_snapshotVect.begin());
}

unsigned int Trajectory::removeSnapshot(unsigned int snapshotNum) {
  PRECONDITION(snapshotNum < d_snapshotVect.size(), "snapshotNum out of bounds");
  delete d_snapshotVect[snapshotNum];
  return (d_snapshotVect.erase(d_snapshotVect.begin() + snapshotNum) - d_snapshotVect.begin());
}

unsigned int Trajectory::readAmber(const std::string &fName) {
  PRECONDITION(d_dimension == 3, "The trajectory must have dimension == 3");
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
  unsigned int nCoords = d_numPoints * 3;
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
        }
        else if (i && (i < (nCoords - 1))) {
          std::stringstream ss;
          ss << "Premature end of file: " << fName;
          throw ValueErrorException(ss.str());
        }
        break;
      }
      ++i;
    }
    if (!inStream.eof()) {
      addSnapshot(new Snapshot(c));
      ++nSnapshots;
    }
  }
  return nSnapshots;
}

unsigned int Trajectory::readGromos(const std::string &fName) {
  PRECONDITION(d_dimension == 3, "The trajectory must have dimension == 3");
  std::ifstream inStream(fName.c_str());
  if (!inStream || inStream.bad()) {
    std::stringstream ss;
    ss << "Bad input file: " << fName;
    throw RDKit::BadFileException(ss.str());
  }
  std::string tempStr;
  unsigned int nCoords = d_numPoints * 3;
  unsigned int nSnapshots = 0;
  const static char *ignoredKeywordArray[] = {
    "TITLE",
    "TIMESTEP",
    "VELOCITYRED",
    "VELOCITY",
    "GENBOX",
    "BOX"
  };
  std::set<std::string> ignoredKeywordSet;
  for (unsigned int i = 0; i < (sizeof(ignoredKeywordArray) / sizeof(char *)); ++i)
    ignoredKeywordSet.insert(std::string(ignoredKeywordArray[i]));
  while (inStream.good() && !inStream.eof()) {
    std::getline(inStream, tempStr);
    if (inStream.bad() || inStream.eof())
      continue;
    if (ignoredKeywordSet.find(tempStr) != ignoredKeywordSet.end()) {
      // ignored block
      while (inStream.good() && !inStream.eof() && (tempStr != "END"))
        std::getline(inStream, tempStr);
    }
    else if ((tempStr == "POSITIONRED") || (tempStr == "POSITION")) {
      // these are the positions
      boost::shared_array<double> c(new double[nCoords]());
      unsigned int j = 0;
      for (unsigned int i = 0; i < d_numPoints;) {
        std::getline(inStream, tempStr);
        if (inStream.bad() || inStream.eof() || (tempStr == "END"))
          throw ValueErrorException("Wrong number of coordinates");
        // ignore comments
        if (tempStr.find("#") != std::string::npos)
          continue;
        std::stringstream ls(tempStr);
        double x, y, z;
        if (!(ls >> x >> y >> z))
          throw ValueErrorException("Error while reading file");
        // store the coordinates (convert to Angstrom!)
        c[j++] = x * 10.0;
        c[j++] = y * 10.0;
        c[j++] = z * 10.0;
        ++i;
      }
      std::getline(inStream, tempStr); // the END line
      if (inStream.bad() || inStream.eof() || (tempStr != "END"))
        throw ValueErrorException("Wrong number of coordinates");
      addSnapshot(new Snapshot(c));
      ++nSnapshots;
    }
    else {
      std::string supportedBlocks("POSITIONRED, POSITION");
      for (std::set<std::string>::const_iterator it = ignoredKeywordSet.begin();
        it != ignoredKeywordSet.end(); ++it)
        supportedBlocks += ", " + *it;
      throw ValueErrorException("Unsupported block: "
        + tempStr + ". Supported blocks are " + supportedBlocks);
    }
  } // read file
  if (inStream.bad()) {
    std::stringstream ss;
		ss << "Bad input file: " << fName;
		throw RDKit::BadFileException(ss.str());
  }
  return nSnapshots;
}

}
