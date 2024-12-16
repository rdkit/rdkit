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
#include <RDGeneral/export.h>
#ifndef UNIFORMREALVALUEGRID3D_H_20140403
#define UNIFORMREALVALUEGRID3D_H_20140403

#include <DataStructs/RealValueVect.h>
#include <Geometry/point.h>
#include "Grid3D.h"

namespace RDGeom {
class RDKIT_RDGEOMETRYLIB_EXPORT UniformRealValueGrid3D
    : public Grid3D<RDKit::RealValueVect, double, double> {
 public:
  UniformRealValueGrid3D() { initGrid(1.0, 1.0, 1.0, 1.0, RDGeom::Point3D()); };
  //! \brief ctor
  /*
      \param dimX:    the x dimension of the grid, in Angstroms
      \param dimY:    the y dimension of the grid, in Angstroms
      \param dimZ:    the z dimension of the grid, in Angstroms
      \param spacing: the grid spacing, in Angstroms
      \param offset:  OPTIONAL: the offset of the grid from (0,0,0), in
     Angstroms.

      \b Note: the values of arguments such as \c dimX and \c spacing
      don't actually need to be in Angstroms, but they should be internally
      consistent.

   */
  UniformRealValueGrid3D(double dimX, double dimY, double dimZ,
                         double spacing = 0.5,
                         const RDGeom::Point3D *offset = nullptr,
                         const RDKit::RealValueVect *data = nullptr) {
    if (offset == nullptr) {
      initGrid(dimX, dimY, dimZ, spacing,
               RDGeom::Point3D(-0.5 * dimX, -0.5 * dimY, -0.5 * dimZ), data);
    } else {
      initGrid(dimX, dimY, dimZ, spacing, *offset, data);
    }
  }
  UniformRealValueGrid3D(const UniformRealValueGrid3D &other);
  UniformRealValueGrid3D &operator=(const UniformRealValueGrid3D &other);

  //! construct from a string pickle
  UniformRealValueGrid3D(const std::string &pkl);
  //! construct from a text pickle
  UniformRealValueGrid3D(const char *pkl, unsigned int);

  //! \brief Get the index of the grid point closest to point
  //!
  //! \return the integer index, -1 if the specified point is outside the grid
  int getGridPointIndex(const RDGeom::Point3D &point) const override;

  //! \brief Get the value at the grid point closest to the specified point
  //!
  //! \return the double value, -1 if the specified index is outside the grid
  double getVal(const RDGeom::Point3D &point) const override;

  //! \brief Get the value at a specified grid point
  //!
  //! \return the double value
  double getVal(unsigned int pointId) const override;

  //! \brief Set the value at the grid point closest to the specified point
  //!
  //! doesn't do anything if the point is outside the grid
  void setVal(const RDGeom::Point3D &point, double val) override;

  //! \brief get the location of the specified grid point
  RDGeom::Point3D getGridPointLoc(unsigned int pointId) const override;

  //! \brief Set the value at the specified grid point
  void setVal(unsigned int pointId, double val) override;

  //! \brief get the size of the grid (number of grid points)
  unsigned int getSize() const override { return d_numX * d_numY * d_numZ; };

  //! \brief get the index of the grid point given the x, y, z indices
  //!
  //! \return the integer value, -1 if the indices are outside the grid
  int getGridIndex(unsigned int xi, unsigned int yi, unsigned int zi) const;

  //! \brief get the x, y, and z indices of a grid-point index
  //!
  void getGridIndices(unsigned int idx, unsigned int &xi, unsigned int &yi,
                      unsigned int &zi) const;

  //! \brief get the number of grid points along x-axis
  unsigned int getNumX() const { return d_numX; };

  //! \brief get the number of grid points along y-axis
  unsigned int getNumY() const { return d_numY; };

  //! \brief get the number of grid points along z-axis
  unsigned int getNumZ() const { return d_numZ; };

  //! \brief get the grid's offset
  const RDGeom::Point3D &getOffset() const { return d_offSet; };

  //! \brief get the grid's spacing
  double getSpacing() const { return d_spacing; };

  //! \brief return a \b const pointer to our occupancy vector
  const RDKit::RealValueVect *getOccupancyVect() const override {
    return &d_storage;
  };

  //! brief returns raw vector
  const std::vector<double> &getData() const { return d_storage.getData(); }
  std::vector<double> &getData() { return d_storage.getData(); }

  //! \brief returns true if the grid \c other has parameters
  //!        compatible with ours.
  bool compareParams(const UniformRealValueGrid3D &other) const;

  //! \brief returns true if the grid \c other has the same values
  //!        as ours.
  bool compareVectors(const UniformRealValueGrid3D &other) const;

  //! \brief returns true if the grid \c other has the same parameters and
  //!        values as ours.
  bool compareGrids(const UniformRealValueGrid3D &other) const;

  //! \brief calculates the union between the data on this grid and
  //!  that on \c other.
  //!  This grid is modified.
  //!  NOTE that the grids must have the same parameters.
  UniformRealValueGrid3D &operator|=(const UniformRealValueGrid3D &other);
  //! \brief calculates the intersection between the data on this grid and
  //!  that on \c other.
  //!  This grid is modified.
  //!  NOTE that the grids must have the same parameters.
  UniformRealValueGrid3D &operator&=(const UniformRealValueGrid3D &other);
  //! \brief calculates the sum of the data on this grid and
  //!  that on \c other.
  //!  This grid is modified.
  //!  NOTE that the grids must have the same parameters.
  UniformRealValueGrid3D &operator+=(const UniformRealValueGrid3D &other);
  //! \brief calculates the difference between the data on this grid and
  //!  that on \c other.
  //!  This grid is modified.
  //!  NOTE that the grids must have the same parameters.
  UniformRealValueGrid3D &operator-=(const UniformRealValueGrid3D &other);

  //! \brief create and return a pickle
  std::string toString() const;

  /*
  UniformRealValueGrid3D operator& (const UniformRealValueGrid3D &other) const{
    PRECONDITION(dp_storage,"bad storage");
    PRECONDITION(compareParams(other),"mismatched params");
    UniformRealValueGrid3D
  res(d_numX*d_spacing,d_numY*d_spacing,d_numZ*d_spacing, d_spacing,&d_offSet);
    return res;
  };*/

 private:
  //! \brief internal initialization code
  /*
      \param dimX:    the x dimension of the grid, in Angstroms
      \param dimY:    the y dimension of the grid, in Angstroms
      \param dimZ:    the z dimension of the grid, in Angstroms
      \param spacing: the grid spacing, in Angstroms
      \param offset:  the offset of the grid from (0,0,0), in Angstroms.
      \param data:    (optional) a pointer to a DoubleVector to use, we take
                      ownership of the pointer.
   */
  void initGrid(double dimX, double dimY, double dimZ, double spacing,
                const RDGeom::Point3D &offSet,
                const RDKit::RealValueVect *data = nullptr);

  unsigned int d_numX, d_numY,
      d_numZ;        //! number of grid points along x, y, z axes
  double d_spacing;  //! grid spacing

  RDGeom::Point3D d_offSet;        //! the grid offset (from the origin)
  RDKit::RealValueVect d_storage;  //! storage for values at each grid point

  //! \brief construct from a pickle
  void initFromText(const char *pkl, const unsigned int length);
};

RDKIT_RDGEOMETRYLIB_EXPORT UniformRealValueGrid3D operator|(
    const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
RDKIT_RDGEOMETRYLIB_EXPORT UniformRealValueGrid3D operator&(
    const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
RDKIT_RDGEOMETRYLIB_EXPORT UniformRealValueGrid3D operator+(
    const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
RDKIT_RDGEOMETRYLIB_EXPORT UniformRealValueGrid3D operator-(
    const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
}  // namespace RDGeom

#endif
