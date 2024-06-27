// $Id$
//
//   Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "UniformGrid3D.h"
#include <DataStructs/DiscreteValueVect.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/Exceptions.h>
#include "point.h"
#include <fstream>
#include <cstdint>

constexpr double OFFSET_TOL = 1.e-8;
constexpr double SPACING_TOL = 1.e-8;

using namespace RDKit;

namespace RDGeom {
unsigned int ci_GRIDPICKLE_VERSION = 0x1;

UniformGrid3D::UniformGrid3D(const UniformGrid3D &other) : Grid3D(other) {
  PRECONDITION(other.dp_storage, "cannot copy an uninitialized grid");
  auto *data = new RDKit::DiscreteValueVect(*other.dp_storage);
  initGrid(other.d_numX * other.d_spacing, other.d_numY * other.d_spacing,
           other.d_numZ * other.d_spacing, other.d_spacing,
           other.dp_storage->getValueType(), other.d_offSet, data);
}

UniformGrid3D &UniformGrid3D::operator=(const UniformGrid3D &other) {
  if (&other == this) {
    return *this;
  }
  PRECONDITION(other.dp_storage, "cannot copy an uninitialized grid");
  delete dp_storage;
  auto *data = new RDKit::DiscreteValueVect(*other.dp_storage);
  initGrid(other.d_numX * other.d_spacing, other.d_numY * other.d_spacing,
           other.d_numZ * other.d_spacing, other.d_spacing,
           other.dp_storage->getValueType(), other.d_offSet, data);
  return *this;
}

void UniformGrid3D::initGrid(
    double dimX, double dimY, double dimZ, double spacing,
    RDKit::DiscreteValueVect::DiscreteValueType valType,
    const RDGeom::Point3D &offSet, RDKit::DiscreteValueVect *data) {
  PRECONDITION(dimX > 0.0, "Invalid x-dimension for grid");
  PRECONDITION(dimY > 0.0, "Invalid y-dimension for grid");
  PRECONDITION(dimZ > 0.0, "Invalid z-dimension for grid");
  PRECONDITION(spacing > 0.0, "Invalid spacing for grid");
  d_numX = static_cast<unsigned int>(floor(dimX / spacing + 0.5));
  d_numY = static_cast<unsigned int>(floor(dimY / spacing + 0.5));
  d_numZ = static_cast<unsigned int>(floor(dimZ / spacing + 0.5));
  PRECONDITION((!data) || data->getValueType() == valType,
               "grid data type mismatch");
  PRECONDITION((!data) || data->getLength() == d_numX * d_numY * d_numZ,
               "grid data size mismatch");

  d_spacing = spacing;
  d_offSet = offSet;
  if (!data) {
    dp_storage =
        new RDKit::DiscreteValueVect(valType, d_numX * d_numY * d_numZ);
  } else {
    dp_storage = data;
  }
}

UniformGrid3D::UniformGrid3D(const std::string &pkl) {
  dp_storage = nullptr;
  this->initFromText(pkl.c_str(), pkl.size());
}
UniformGrid3D::UniformGrid3D(const char *pkl, const unsigned int len) {
  dp_storage = nullptr;
  this->initFromText(pkl, len);
}

UniformGrid3D::~UniformGrid3D() {
  delete dp_storage;
  dp_storage = nullptr;
}

int UniformGrid3D::getGridIndex(unsigned int xi, unsigned int yi,
                                unsigned int zi) const {
  if (xi >= d_numX || yi >= d_numY || zi >= d_numZ) {
    return -1;
  }

  return (zi * d_numX * d_numY + yi * d_numX + xi);
}

void UniformGrid3D::getGridIndices(unsigned int idx, unsigned int &xi,
                                   unsigned int &yi, unsigned int &zi) const {
  if (idx >= d_numX * d_numY * d_numZ) {
    throw IndexErrorException(idx);
  }
  xi = idx % d_numX;
  yi = (idx % (d_numX * d_numY)) / d_numX;
  zi = idx / (d_numX * d_numY);
}

int UniformGrid3D::getGridPointIndex(const Point3D &point) const {
  Point3D tPt(point);
  tPt -= d_offSet;  // d_origin;
  tPt /= d_spacing;
  double move = 0.5;
  int xi = static_cast<int>(floor(tPt.x + move));
  int yi = static_cast<int>(floor(tPt.y + move));
  int zi = static_cast<int>(floor(tPt.z + move));

  if ((xi < 0) || (xi >= static_cast<int>(d_numX))) {
    return -1;
  }
  if ((yi < 0) || (yi >= static_cast<int>(d_numY))) {
    return -1;
  }
  if ((zi < 0) || (zi >= static_cast<int>(d_numZ))) {
    return -1;
  }

  return (zi * d_numX * d_numY + yi * d_numX + xi);
}

int UniformGrid3D::getVal(const Point3D &point) const {
  int id = getGridPointIndex(point);
  if (id < 0) {
    return -1;
  }
  return dp_storage->getVal(static_cast<unsigned int>(id));
}

unsigned int UniformGrid3D::getVal(unsigned int pointId) const {
  return dp_storage->getVal(pointId);
}

void UniformGrid3D::setVal(const Point3D &point, unsigned int val) {
  int id = getGridPointIndex(point);
  if (id < 0) {
    return;
  }
  dp_storage->setVal(static_cast<unsigned int>(id), val);
}

Point3D UniformGrid3D::getGridPointLoc(unsigned int pointId) const {
  if (pointId >= d_numX * d_numY * d_numZ) {
    throw IndexErrorException(pointId);
  }
  Point3D res;
  res.x = (pointId % d_numX) * d_spacing;
  // the rounding here is intentional, we want the coordinates of a grid point
  res.y =
      static_cast<double>((pointId % (d_numX * d_numY)) / d_numX) * d_spacing;
  res.z = static_cast<double>(pointId / (d_numX * d_numY)) * d_spacing;
  res += d_offSet;  // d_origin;
  return res;
}

void UniformGrid3D::setVal(unsigned int pointId, unsigned int val) {
  dp_storage->setVal(pointId, val);
}

bool UniformGrid3D::compareParams(const UniformGrid3D &other) const {
  if (d_numX != other.getNumX() || d_numY != other.getNumY() ||
      d_numZ != other.getNumZ()) {
    return false;
  }

  if (fabs(d_spacing - other.getSpacing()) > SPACING_TOL) {
    return false;
  }
  Point3D dOffset = d_offSet;
  dOffset -= other.getOffset();
  return dOffset.lengthSq() <= OFFSET_TOL;
}

void UniformGrid3D::setSphereOccupancy(const Point3D &center, double radius,
                                       double stepSize, int maxNumLayers,
                                       bool ignoreOutOfBound) {
  int ptIndex = this->getGridPointIndex(center);
  if (ptIndex == -1) {
    if (ignoreOutOfBound) {
      return;
    } else {
      throw GridException("Center outside the grid boundary");
    }
  }
  Point3D gPt(center);  // point on the grid
  gPt -= d_offSet;

  gPt /= d_spacing;

  // unsigned int z = ptIndex/(d_numX*d_numY);
  // unsigned int y = (ptIndex%(d_numX*d_numY))/d_numX;
  // unsigned int x = ptIndex%d_numX;

  unsigned int bPerVal = dp_storage->getNumBitsPerVal();
  unsigned int maxVal = (1 << bPerVal) - 1;
  unsigned int nLayers = maxVal;
  unsigned int valStep = 1;
  if ((maxNumLayers > 0) && (maxNumLayers <= static_cast<int>(nLayers))) {
    nLayers = maxNumLayers;
    valStep = (maxVal + 1) / nLayers;
  }
  double bgRad = radius / d_spacing;        // base radius in grid coords
  double gStepSize = stepSize / d_spacing;  // step size in grid coords
  double gRadius =
      bgRad + nLayers * gStepSize;  // largest radius in grid coords
  double gRad2 = gRadius * gRadius;
  double bgRad2 = bgRad * bgRad;
  double dx, dy, dz, d, d2, dy2z2, dz2;
  int xmax = static_cast<int>(floor(gPt.x + gRadius));
  int xmin = static_cast<int>(ceil(gPt.x - gRadius));
  int ymax = static_cast<int>(floor(gPt.y + gRadius));
  int ymin = static_cast<int>(ceil(gPt.y - gRadius));
  int zmax = static_cast<int>(floor(gPt.z + gRadius));
  int zmin = static_cast<int>(ceil(gPt.z - gRadius));

  unsigned int oval, val, valChange;
  int ptId1, ptId2;
  for (int k = zmin; k <= zmax; ++k) {
    if ((k >= 0) &&
        (k < static_cast<int>(d_numZ))) {  // we are inside the grid in the z-direction
      dz = static_cast<double>(k) - gPt.z;
      dz2 = dz * dz;
      ptId1 = k * d_numX * d_numY;
      for (int j = ymin; j <= ymax; ++j) {
        if ((j >= 0) &&
            (j < static_cast<int>(d_numY))) {  // inside the grid in the y-direction
          dy = static_cast<double>(j) - gPt.y;
          dy2z2 = dy * dy + dz2;
          if (dy2z2 < gRad2) {  // we are within the radius at least from the
                                // y,z coordinates
            ptId2 = ptId1 + j * d_numX;
            for (int i = xmin; i <= xmax; ++i) {
              if ((i >= 0) && (i < static_cast<int>(d_numX))) {
                oval = dp_storage->getVal(static_cast<unsigned int>(ptId2 + i));
                if (oval < maxVal) {  // if we are already at maxVal we will not
                                      // change that
                  dx = static_cast<double>(i) - gPt.x;
                  d2 = dx * dx + dy2z2;
                  if (d2 < gRad2) {  // we are inside the sphere
                    if (d2 < bgRad2) {
                      val = maxVal;
                    } else {
                      d = sqrt(d2);
                      valChange = (static_cast<unsigned int>(
                                      (d - bgRad) / gStepSize + 1)) *
                                  (valStep);
                      if (valChange < maxVal) {
                        val = maxVal - valChange;
                      } else {
                        val = 0;
                      }
                    }
                    if (val > oval) {
                      dp_storage->setVal(ptId2 + i, val);
                    }
                  }  // we are inside the sphere
                }    // grid point does not already have maxVal
              }      // inside the grid in x-direction
            }        // loop over points in x-direction
          }          // inside the sphere based on only z and y coords
        }            // we are inside the grid in the y-direction
      }              // loop over points in y-direction
    }                // inside grid in z-direction
  }                  // loop over points in z-direction
}

UniformGrid3D &UniformGrid3D::operator|=(const UniformGrid3D &other) {
  PRECONDITION(dp_storage, "uninitialized grid");
  PRECONDITION(other.dp_storage, "uninitialized grid");
  PRECONDITION(compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  RDKit::DiscreteValueVect *newData =
      new RDKit::DiscreteValueVect((*dp_storage) | (*other.dp_storage));
  delete dp_storage;
  dp_storage = newData;
  return *this;
}

UniformGrid3D &UniformGrid3D::operator&=(const UniformGrid3D &other) {
  PRECONDITION(dp_storage, "uninitialized grid");
  PRECONDITION(other.dp_storage, "uninitialized grid");
  PRECONDITION(compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  RDKit::DiscreteValueVect *newData =
      new RDKit::DiscreteValueVect((*dp_storage) & (*other.dp_storage));
  delete dp_storage;
  dp_storage = newData;
  return *this;
}

UniformGrid3D &UniformGrid3D::operator+=(const UniformGrid3D &other) {
  PRECONDITION(dp_storage, "uninitialized grid");
  PRECONDITION(other.dp_storage, "uninitialized grid");
  PRECONDITION(compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  *dp_storage += *other.dp_storage;
  return *this;
}

UniformGrid3D &UniformGrid3D::operator-=(const UniformGrid3D &other) {
  PRECONDITION(dp_storage, "uninitialized grid");
  PRECONDITION(other.dp_storage, "uninitialized grid");
  PRECONDITION(compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  *dp_storage -= *other.dp_storage;
  return *this;
}

std::string UniformGrid3D::toString() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  std::int32_t tVers = ci_GRIDPICKLE_VERSION * -1;
  streamWrite(ss, tVers);
  std::uint32_t tInt;
  tInt = d_numX;
  streamWrite(ss, tInt);
  tInt = d_numY;
  streamWrite(ss, tInt);
  tInt = d_numZ;
  streamWrite(ss, tInt);
  streamWrite(ss, d_spacing);
  streamWrite(ss, d_offSet.x);
  streamWrite(ss, d_offSet.y);
  streamWrite(ss, d_offSet.z);

  std::string storePkl = dp_storage->toString();
  std::uint32_t pklSz = storePkl.size();
  streamWrite(ss, pklSz);
  ss.write(storePkl.c_str(), pklSz * sizeof(char));

  std::string res(ss.str());
  return (res);
}
void UniformGrid3D::initFromText(const char *pkl, const unsigned int length) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pkl, length);
  std::int32_t tVers;
  streamRead(ss, tVers);
  tVers *= -1;
  if (tVers == 0x1) {
  } else {
    throw ValueErrorException("bad version in UniformGrid3D pickle");
  }
  std::uint32_t tInt;
  streamRead(ss, tInt);
  d_numX = tInt;
  streamRead(ss, tInt);
  d_numY = tInt;
  streamRead(ss, tInt);
  d_numZ = tInt;
  streamRead(ss, d_spacing);
  double oX, oY, oZ;
  streamRead(ss, oX);
  streamRead(ss, oY);
  streamRead(ss, oZ);
  d_offSet = Point3D(oX, oY, oZ);

  std::uint32_t pklSz;
  streamRead(ss, pklSz);
  auto *buff = new char[pklSz];
  ss.read(buff, pklSz * sizeof(char));
  delete dp_storage;
  dp_storage = new RDKit::DiscreteValueVect(buff, pklSz);
  delete[] buff;
}

void writeGridToStream(const UniformGrid3D &grid, std::ostream &outStrm) {
  int dimX = static_cast<int>(grid.getNumX());  //+2;
  int dimY = static_cast<int>(grid.getNumY());  //+2;
  int dimZ = static_cast<int>(grid.getNumZ());  //+2;
  double spacing = grid.getSpacing();
  double lenX = dimX * spacing;
  double lenY = dimY * spacing;
  double lenZ = dimZ * spacing;
  Point3D offSet = grid.getOffset();
  offSet /= spacing;
  outStrm << "Grid file representing a Shape \n\n";
  outStrm << lenX << " " << lenY << " " << lenZ << " 90.0 90.0 90.0"
          << std::endl;
  outStrm << dimX - 1 << " " << dimY - 1 << " " << dimZ - 1 << std::endl;

  int outX1 = static_cast<int>(floor(offSet.x + 0.5));
  int outX2 = static_cast<int>(floor(offSet.x + 0.5)) + static_cast<int>(dimX - 1);
  // REVIEW: ok - here is a fix to try and make the grid closer to the molecule
  // when displayed
  // (at least in PyMol). The difference between the pair of values (outX1,
  // outX2) is (dimX-1),
  // that is not the case with pairs the (outY1, outY2) a (outZ1, outZ2). In
  // these cases, the difference is
  // dimY and dimZ respectively. Not sure why this should be the case, but
  // almost always we get a
  // better display in PyMol.
  int outY1 = static_cast<int>(floor(offSet.y + 0.5));
  int outY2 = static_cast<int>(floor(offSet.y + 0.5)) + static_cast<int>(dimY);
  int outZ1 = static_cast<int>(floor(offSet.z + 0.5));
  int outZ2 = static_cast<int>(floor(offSet.z + 0.5)) + static_cast<int>(dimZ);

  outStrm << "1"
          << " " << outX1 << " " << outX2 << " " << outY1 << " " << outY2 << " "
          << outZ1 << " " << outZ2 << "\n";
  unsigned int i, nPts = grid.getSize();
  for (i = 0; i < nPts; i++) {
    outStrm << static_cast<double>(grid.getVal(i)) << std::endl;
  }
}

void writeGridToFile(const UniformGrid3D &grid, const std::string &filename) {
  // std::ofstream ofStrm(filename.c_str());
  auto *ofStrm = new std::ofstream(filename.c_str());
  auto *oStrm = static_cast<std::ostream *>(ofStrm);
  writeGridToStream(grid, *oStrm);
  delete ofStrm;
}
}  // namespace RDGeom
