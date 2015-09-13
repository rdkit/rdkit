//  UniformRealValueGrid3D.h
//  Created on: Apr 3, 2014
//  Author: hahnda6
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _UNIFORMREALVALUEGRID3D_H_20140403
#define _UNIFORMREALVALUEGRID3D_H_20140403

#include <DataStructs/RealValueVect.h>
#include <Geometry/point.h>
#include "Grid3D.h"

namespace RDGeom {
  class UniformRealValueGrid3D : public Grid3D<RDKit::RealValueVect, double, double> {

  public:
    UniformRealValueGrid3D() {
      initGrid(1.0,1.0,1.0,1.0,RDGeom::Point3D());
    };
    //! \brief ctor
    /*
        \param dimX:    the x dimension of the grid, in Angstroms
        \param dimY:    the y dimension of the grid, in Angstroms
        \param dimZ:    the z dimension of the grid, in Angstroms
        \param spacing: the grid spacing, in Angstroms
        \param offset:  OPTIONAL: the offset of the grid from (0,0,0), in Angstroms.

        \b Note: the values of arguments such as \c dimX and \c spacing
        don't actually need to be in Angstroms, but they should be internally
        consistent.

     */
    UniformRealValueGrid3D(double dimX, double dimY, double dimZ, double spacing=0.5,
                           const RDGeom::Point3D *offset=0, RDKit::RealValueVect *data=0){
      if (offset == 0) {
        initGrid(dimX, dimY, dimZ, spacing,
                 RDGeom::Point3D(-0.5*dimX, -0.5*dimY, -0.5*dimZ), data);
      } else {
        initGrid(dimX, dimY, dimZ, spacing, *offset, data);
      }
    }
    //! copy ctor
    UniformRealValueGrid3D(const UniformRealValueGrid3D &other);
    //! construct from a string pickle
    UniformRealValueGrid3D(const std::string &pkl);
    //! construct from a text pickle
    UniformRealValueGrid3D(const char *pkl,unsigned int);

    ~UniformRealValueGrid3D();

    //! \brief Get the index of the grid point closest to point
    //!
    //! \return the integer index, -1 if the specified point is outside the grid
    int getGridPointIndex(const RDGeom::Point3D &point) const;

    //! \brief Get the value at the grid point closest to the specified point
    //!
    //! \return the double value, -1 if the specified index is outside the grid
    double getVal(const RDGeom::Point3D &point) const;

    //! \brief Get the value at a specified grid point
    //!
    //! \return the double value
    double getVal(unsigned int pointId) const;

    //! \brief Set the value at the grid point closest to the specified point
    //!
    //! doesn't do anything if the point is outside the grid 
    void setVal(const RDGeom::Point3D &point, double val);

    //! \brief get the location of the specified grid point
    RDGeom::Point3D getGridPointLoc(unsigned int pointId) const;

    //! \brief Set the value at the specified grid point 
    void setVal(unsigned int pointId, double val);

    //! \brief get the size of the grid (number of grid points)
    unsigned int getSize() const { return d_numX*d_numY*d_numZ; };


    //! \brief get the index of the grid point given the x, y, z indices
    //!
    //! \return the integer value, -1 if the indices are outside the grid
    int getGridIndex(unsigned int xi, unsigned int yi, unsigned int zi) const;

    //! \brief get the x, y, and z indices of a grid-point index
    //!
    void getGridIndices(unsigned int idx,unsigned int &xi, unsigned int &yi, unsigned int &zi) const;


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
    const RDKit::RealValueVect *getOccupancyVect() const { return dp_storage;} ;

    //!brief returns shared pointer
    const boost::shared_array<double> &getDataPtr() {
      return dp_storage->getArray();
    }
    
    //! \brief returns true if the grid \c other has parameters
    //!        compatible with ours.
    bool compareParams(const UniformRealValueGrid3D &other) const;

    //! \brief returns true if the grid \c other has the same values
    //!        as ours.
    bool compareVectors(const UniformRealValueGrid3D &other) const;

    //! \brief returns true if the grid \c other has the same parameters and
    //!        values as ours.
    bool compareGrids(const UniformRealValueGrid3D &other) const;

    //! \brief copies a grid \c other into grid
    UniformRealValueGrid3D &operator=(const UniformRealValueGrid3D &other);
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
      UniformRealValueGrid3D res(d_numX*d_spacing,d_numY*d_spacing,d_numZ*d_spacing,
			d_spacing,&d_offSet);
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
                  const RDGeom::Point3D &offSet, RDKit::RealValueVect *data=0);


    unsigned int d_numX, d_numY, d_numZ; //! number of grid points along x, y, z axes
    double d_spacing; //! grid spacing

    RDGeom::Point3D d_offSet; //! the grid offset (from the origin)
    RDKit::RealValueVect *dp_storage; //! storage for values at each grid point

    //! \brief construct from a pickle
    void initFromText(const char *pkl,const unsigned int length);

  };

  UniformRealValueGrid3D operator| (const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
  UniformRealValueGrid3D operator& (const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
  UniformRealValueGrid3D operator+ (const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
  UniformRealValueGrid3D operator- (const UniformRealValueGrid3D &grd1, const UniformRealValueGrid3D &grd2);
}

#endif

