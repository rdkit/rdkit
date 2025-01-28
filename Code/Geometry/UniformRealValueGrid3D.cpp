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

using namespace RDKit;

namespace {
constexpr double OFFSET_TOL = 1.e-4;
constexpr double SPACING_TOL = 1.e-4;
constexpr const char *INCOMPATIBLE_GRIDS = "incompatible grids";
}  // end anonymous namespace

namespace RDGeom {
std::int32_t ci_RealValueGrid3DPICKLE_VERSION = 0x1;
UniformRealValueGrid3D::UniformRealValueGrid3D(
    const UniformRealValueGrid3D &other) {
  UniformRealValueGrid3D::initGrid(
      other.d_numX * other.d_spacing, other.d_numY * other.d_spacing,
      other.d_numZ * other.d_spacing, other.d_spacing, other.d_offSet,
      &other.d_storage);
}

void UniformRealValueGrid3D::initGrid(double dimX, double dimY, double dimZ,
                                      double spacing,
                                      const RDGeom::Point3D &offSet,
                                      const RDKit::RealValueVect *data) {
  PRECONDITION(dimX > 0.0, "Invalid x-dimension for grid");
  PRECONDITION(dimY > 0.0, "Invalid y-dimension for grid");
  PRECONDITION(dimZ > 0.0, "Invalid z-dimension for grid");
  PRECONDITION(spacing > 0.0, "Invalid spacing for grid");
  d_numX = static_cast<unsigned int>(floor(dimX / spacing + 0.5));
  d_numY = static_cast<unsigned int>(floor(dimY / spacing + 0.5));
  d_numZ = static_cast<unsigned int>(floor(dimZ / spacing + 0.5));
  // PRECONDITION((!data)||data->getValueType()==double,"grid data type
  // mismatch");
  PRECONDITION(!data || data->size() == d_numX * d_numY * d_numZ,
               "grid data size mismatch");

  d_spacing = spacing;
  d_offSet = offSet;
  if (!data) {
    d_storage.setLength(d_numX * d_numY * d_numZ);
    d_storage.setToVal(0.0);
  } else {
    d_storage = *data;
  }
}

UniformRealValueGrid3D::UniformRealValueGrid3D(const std::string &pkl) {
  initFromText(pkl.c_str(), pkl.size());
}

UniformRealValueGrid3D::UniformRealValueGrid3D(const char *pkl,
                                               unsigned int len) {
  initFromText(pkl, len);
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
  auto tPt = point;
  tPt -= d_offSet;  // d_origin;
  tPt /= d_spacing;
  constexpr double move = 0.5;
  auto xi = static_cast<int>(floor(tPt.x + move));
  auto yi = static_cast<int>(floor(tPt.y + move));
  auto zi = static_cast<int>(floor(tPt.z + move));

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
  auto id = getGridPointIndex(point);
  if (id < 0) {
    return -1;
  }
  return d_storage.getVal(static_cast<unsigned int>(id));
}

double UniformRealValueGrid3D::getVal(unsigned int pointId) const {
  return d_storage.getVal(pointId);
}

void UniformRealValueGrid3D::setVal(const RDGeom::Point3D &point, double val) {
  auto id = getGridPointIndex(point);
  if (id < 0) {
    return;
  }
  d_storage.setVal(static_cast<unsigned int>(id), val);
}

RDGeom::Point3D UniformRealValueGrid3D::getGridPointLoc(
    unsigned int pointId) const {
  if (pointId >= d_numX * d_numY * d_numZ) {
    throw IndexErrorException(pointId);
  }
  RDGeom::Point3D res((pointId % d_numX) * d_spacing,
                      ((pointId % (d_numX * d_numY)) / d_numX) * d_spacing,
                      (pointId / (d_numX * d_numY)) * d_spacing);
  res += d_offSet;  // d_origin;
  return res;
}

void UniformRealValueGrid3D::setVal(unsigned int pointId, double val) {
  d_storage.setVal(pointId, val);
}

bool UniformRealValueGrid3D::compareParams(
    const UniformRealValueGrid3D &other) const {
  if (d_numX != other.getNumX()) {
    return false;
  }
  if (d_numY != other.getNumY()) {
    return false;
  }
  if (d_numZ != other.getNumZ()) {
    return false;
  }
  if (fabs(d_spacing - other.getSpacing()) > SPACING_TOL) {
    return false;
  }
  auto dOffset = d_offSet;
  dOffset -= other.getOffset();
  return dOffset.lengthSq() <= OFFSET_TOL;
}

bool UniformRealValueGrid3D::compareVectors(
    const UniformRealValueGrid3D &other) const {
  return d_storage.compareVectors(other.d_storage);
}

bool UniformRealValueGrid3D::compareGrids(
    const UniformRealValueGrid3D &other) const {
  if (!compareParams(other)) {
    return false;
  }
  return d_storage.compareVectors(other.d_storage);
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator=(
    const UniformRealValueGrid3D &other) {
  if (this == &other) {
    return *this;
  }
  d_numX = other.d_numX;
  d_numY = other.d_numY;
  d_numZ = other.d_numZ;
  d_spacing = other.d_spacing;
  d_offSet = other.d_offSet;
  d_storage = other.d_storage;

  return *this;
}
UniformRealValueGrid3D &UniformRealValueGrid3D::operator|=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(compareParams(other), INCOMPATIBLE_GRIDS);

  d_storage |= other.d_storage;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator&=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(compareParams(other), INCOMPATIBLE_GRIDS);

  d_storage &= other.d_storage;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator+=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(compareParams(other), INCOMPATIBLE_GRIDS);

  d_storage += other.d_storage;
  return *this;
}

UniformRealValueGrid3D &UniformRealValueGrid3D::operator-=(
    const UniformRealValueGrid3D &other) {
  PRECONDITION(compareParams(other), INCOMPATIBLE_GRIDS);

  d_storage -= other.d_storage;
  return *this;
}

std::string UniformRealValueGrid3D::toString() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  std::int32_t tVers = ci_RealValueGrid3DPICKLE_VERSION * -1;
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
  std::string storePkl = d_storage.toString();
  std::uint32_t pklSz = storePkl.size();
  streamWrite(ss, pklSz);
  ss.write(storePkl.c_str(), pklSz * sizeof(char));

  return ss.str();
}

void UniformRealValueGrid3D::initFromText(const char *pkl,
                                          const unsigned int length) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pkl, length);
  std::int32_t tVers;
  streamRead(ss, tVers);
  tVers = -tVers;
  if (tVers != 0x1) {
    throw ValueErrorException("bad version in UniformRealValueGrid3D pickle");
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
  d_offSet = RDGeom::Point3D(oX, oY, oZ);

  std::uint32_t pklSz;
  streamRead(ss, pklSz);
  std::vector<char> buff(pklSz);
  ss.read(buff.data(), pklSz * sizeof(char));
  d_storage = RDKit::RealValueVect(buff.data(), pklSz);
}

UniformRealValueGrid3D operator&(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D ans(grd1);
  ans &= grd2;
  return ans;
};

UniformRealValueGrid3D operator|(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D ans(grd1);
  ans |= grd2;
  return ans;
};

UniformRealValueGrid3D operator+(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D ans(grd1);
  ans += grd2;
  return ans;
};

UniformRealValueGrid3D operator-(const UniformRealValueGrid3D &grd1,
                                 const UniformRealValueGrid3D &grd2) {
  UniformRealValueGrid3D ans(grd1);
  ans -= grd2;
  return ans;
};

}  // namespace RDGeom
