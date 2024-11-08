//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "UniformRealValueGrid3D.h"
#include <DataStructs/RealValueVect.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/Exceptions.h>
#include <Geometry/point.h>
#include <fstream>
#include <iomanip>
#include <boost/cstdint.hpp>

#define OFFSET_TOL 1.e-4
#define SPACING_TOL 1.e-4
using namespace RDKit;

namespace RDGeom {
boost::int32_t ci_RealValueGrid3DPICKLE_VERSION = 0x1;
UniformRealValueGrid3D::UniformRealValueGrid3D(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(other.dp_storage,
               "cannot copy an uninitialized RealValueGrid3D");
  RDKit::RealValueVect *data = new RDKit::RealValueVect(*other.dp_storage);
  UniformRealValueGrid3D::initGrid(
      other.d_numX * other.d_spacing, other.d_numY * other.d_spacing,
      other.d_numZ * other.d_spacing, other.d_spacing, other.d_offSet, data);
}

void UniformRealValueGrid3D::initGrid(double dimX, double dimY, double dimZ,
                                      double spacing,
                                      const RDGeom::Point3D &offSet,
                                      RDKit::RealValueVect *data) {
  PRECONDITION(dimX > 0.0, "Invalid x-dimension for grid");
  PRECONDITION(dimY > 0.0, "Invalid y-dimension for grid");
  PRECONDITION(dimZ > 0.0, "Invalid z-dimension for grid");
  PRECONDITION(spacing > 0.0, "Invalid spacing for grid");
  d_numX = static_cast<unsigned int>(floor(dimX / spacing + 0.5));
  d_numY = static_cast<unsigned int>(floor(dimY / spacing + 0.5));
  d_numZ = static_cast<unsigned int>(floor(dimZ / spacing + 0.5));
  // PRECONDITION((!data)||data->getValueType()==double,"grid data type
  // mismatch");
  PRECONDITION((!data) || data->size() == d_numX * d_numY * d_numZ,
               "grid data size mismatch");

  d_spacing = spacing;
  d_offSet = offSet;
  if (!data) {
    dp_storage = new RDKit::RealValueVect(d_numX * d_numY * d_numZ);
  } else {
    dp_storage = data;
  }
}

UniformRealValueGrid3D::UniformRealValueGrid3D(const std::string &pkl) {
  dp_storage = 0;
  this->initFromText(pkl.c_str(), pkl.size());
}

UniformRealValueGrid3D::UniformRealValueGrid3D(const char *pkl,
                                               unsigned int len) {
  dp_storage = 0;
  this->initFromText(pkl, len);
}

UniformRealValueGrid3D::~UniformRealValueGrid3D() {
  delete dp_storage;
  dp_storage = NULL;
}

int UniformRealValueGrid3D::getGridIndex(unsigned int xi, unsigned int yi,
                                         unsigned int zi) const {
  if (xi >= d_numX) {
    return -1;
  }
  if (yi >= d_numY) {
    return -1;
  }
  if (zi >= d_numZ) {
    return -1;
  }
  return (zi * d_numX * d_numY + yi * d_numX + xi);
}

void UniformRealValueGrid3D::getGridIndices(unsigned int idx, unsigned int &xi,
                                            unsigned int &yi,
                                            unsigned int &zi) const {
  if (idx >= d_numX * d_numY * d_numZ) {
    throw IndexErrorException(idx);
  }
  xi = idx % d_numX;
  yi = (idx % (d_numX * d_numY)) / d_numX;
  zi = idx / (d_numX * d_numY);
}

int UniformRealValueGrid3D::getGridPointIndex(
    const RDGeom::Point3D &point) const {
  RDGeom::Point3D tPt(point);
  tPt -= d_offSet;  // d_origin;
  tPt /= d_spacing;
  int xi, yi, zi;
  double move = 0.5;
  xi = static_cast<int>(floor(tPt.x + move));
  yi = static_cast<int>(floor(tPt.y + move));
  zi = static_cast<int>(floor(tPt.z + move));

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

double UniformRealValueGrid3D::getVal(const RDGeom::Point3D &point) const {
  int id = getGridPointIndex(point);
  if (id < 0) {
    return -1;
  }
  return dp_storage->getVal(static_cast<unsigned int>(id));
}

double UniformRealValueGrid3D::getVal(unsigned int pointId) const {
  return dp_storage->getVal(pointId);
}

void UniformRealValueGrid3D::setVal(const RDGeom::Point3D &point, double val) {
  int id = getGridPointIndex(point);
  if (id < 0) {
    return;
  }
  dp_storage->setVal(static_cast<unsigned int>(id), val);
}

RDGeom::Point3D UniformRealValueGrid3D::getGridPointLoc(
    unsigned int pointId) const {
  if (pointId >= d_numX * d_numY * d_numZ) {
    throw IndexErrorException(pointId);
  }
  RDGeom::Point3D res;
  res.x = (pointId % d_numX) * d_spacing;
  res.y = ((pointId % (d_numX * d_numY)) / d_numX) * d_spacing;
  res.z = (pointId / (d_numX * d_numY)) * d_spacing;
  res += d_offSet;  // d_origin;
  return res;
}

void UniformRealValueGrid3D::setVal(unsigned int pointId, double val) {
  dp_storage->setVal(pointId, val);
}

bool UniformRealValueGrid3D::compareParams(
    const UniformRealValueGrid3D &other) const {
  if (d_numX != other.getNumX()) return false;
  if (d_numY != other.getNumY()) return false;
  if (d_numZ != other.getNumZ()) return false;
  if (fabs(d_spacing - other.getSpacing()) > SPACING_TOL) return false;
  RDGeom::Point3D dOffset = d_offSet;
  dOffset -= other.getOffset();
  if (dOffset.lengthSq() > OFFSET_TOL) {
    return false;
  }

  return true;
}

bool UniformRealValueGrid3D::compareVectors(
    const UniformRealValueGrid3D &other) const {
  return dp_storage->compareVectors(*(other.dp_storage));
}

bool UniformRealValueGrid3D::compareGrids(
    const UniformRealValueGrid3D &other) const {
  if (!(compareParams(other))) {
    return false;
  }
  return dp_storage->compareVectors(*(other.dp_storage));
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator=(
    const UniformRealValueGrid3D &other) {
  d_numX = other.d_numX;
  d_numY = other.d_numY;
  d_numZ = other.d_numZ;
  d_spacing = other.d_spacing;
  d_offSet = other.d_offSet;
  *dp_storage = *other.dp_storage;

  return *this;
}
UniformRealValueGrid3D &UniformRealValueGrid3D::operator|=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(dp_storage, "unintialized grid");
  PRECONDITION(other.dp_storage, "unintialized grid");
  PRECONDITION(this->compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  RDKit::RealValueVect *newData =
      new RDKit::RealValueVect((*dp_storage) | (*other.dp_storage));
  delete dp_storage;
  dp_storage = newData;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator&=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(dp_storage, "unintialized grid");
  PRECONDITION(other.dp_storage, "unintialized grid");
  PRECONDITION(this->compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  RDKit::RealValueVect *newData =
      new RDKit::RealValueVect((*dp_storage) & (*other.dp_storage));
  delete dp_storage;
  dp_storage = newData;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator+=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(dp_storage, "unintialized grid");
  PRECONDITION(other.dp_storage, "unintialized grid");
  PRECONDITION(this->compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  *dp_storage += *other.dp_storage;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator-=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(dp_storage, "unintialized grid");
  PRECONDITION(other.dp_storage, "unintialized grid");
  PRECONDITION(this->compareParams(other), "incompatible grids");

  // EFF: we're probably doing too much copying here:
  *dp_storage -= *other.dp_storage;
  return *this;
}

std::string UniformRealValueGrid3D::toString() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  boost::int32_t tVers = ci_RealValueGrid3DPICKLE_VERSION * -1;
  streamWrite(ss, tVers);
  boost::uint32_t tInt;
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
  boost::uint32_t pklSz = storePkl.size();
  streamWrite(ss, pklSz);
  ss.write(storePkl.c_str(), pklSz * sizeof(char));

  std::string res(ss.str());
  return (res);
}

void UniformRealValueGrid3D::initFromText(const char *pkl,
                                          const unsigned int length) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pkl, length);
  boost::int32_t tVers;
  streamRead(ss, tVers);
  tVers *= -1;
  if (tVers == 0x1) {
  } else {
    throw ValueErrorException("bad version in UniformRealValueGrid3D pickle");
  }
  boost::uint32_t tInt;
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
  d_offSet = RDGeom::Point3D(oX, oY, oZ);

  boost::uint32_t pklSz;
  streamRead(ss, pklSz);
  char *buff = new char[pklSz];
  ss.read(buff, pklSz * sizeof(char));
  delete dp_storage;
  dp_storage = new RDKit::RealValueVect(buff, pklSz);
  delete[] buff;
}

UniformRealValueGrid3D operator|(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D res(grd1);
  res |= grd2;
  return res;
}

UniformRealValueGrid3D operator&(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D res(grd1);
  res &= grd2;
  return res;
}

UniformRealValueGrid3D operator+(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D res(grd1);
  res += grd2;
  return res;
}

UniformRealValueGrid3D operator-(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D res(grd1);
  res -= grd2;
  return res;
}

}  // namespace RDGeom
