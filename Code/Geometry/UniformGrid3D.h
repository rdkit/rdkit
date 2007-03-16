// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _UNIFORMGRID3D_H_20050124_1703
#define _UNIFORMGRID3D_H_20050124_1703

#include "point.h"
#include <DataStructs/DiscreteValueVect.h>
#include "Grid3D.h"
#include <iostream>
namespace RDGeom {
  class UniformGrid3D : public Grid3D {
  
  public:
       
    //! \brief ctor
    /*
        \param dimX:    the x dimension of the grid, in Angstroms
        \param dimY:    the y dimension of the grid, in Angstroms
        \param dimZ:    the z dimension of the grid, in Angstroms
        \param spacing: the grid spacing, in Angstroms
        \param valType: the data type of the grid (determines the number of bits
                        per point)
        \param offset:  the offset of the grid from (0,0,0), in Angstroms.

        \b Note: the values of arguments such as \c dimX and \c spacing
        don't actually need to be in Angstroms, but they should be internally
        consistent.
        
    */
    UniformGrid3D(double dimX, double dimY, double dimZ, double spacing=0.5,
                  RDKit::DiscreteValueVect::DiscreteValueType valType=RDKit::DiscreteValueVect::TWOBITVALUE,
                  const RDGeom::Point3D *offset=0) {
      if (offset == 0) {
        initGrid(dimX, dimY, dimZ, spacing, valType,
                 Point3D(-0.5*dimX, -0.5*dimY, -0.5*dimZ));
      } else {
        initGrid(dimX, dimY, dimZ, spacing, valType, *offset);
      }
    }
    //! construct from a string pickle
    UniformGrid3D(const std::string pkl);
    //! construct from a text pickle
    UniformGrid3D(const char *pkl,unsigned int);

    ~UniformGrid3D();

    //! \brief Get the index of the grid point closest to point
    //!
    //! \return the integer index, -1 if the specified point is outside the grid
    int getGridPointIndex(const Point3D &point) const;

    //! \brief Get the value at the grid point closest to the specified point
    //!
    //! \return the integer value, -1 if the specified index is outside the grid
    int getVal(const Point3D &point) const;

    //! \brief Get the value at a specified grid point
    //!
    //! \return the unsigned integer value
    unsigned int getVal(unsigned int pointId) const;

    //! \brief Set the value at the grid point closest to the specified point
    //!
    //! doesn't do anything if the point is outside the grid 
    void setVal(const Point3D &point, unsigned int val);

    //! \brief get the location of the specified grid point
    Point3D getGridPointLoc(unsigned int pointId) const;

    //! \brief Set the value at the specified grid point 
    void setVal(unsigned int pointId, unsigned int val);

    //! \brief get the size of the grid (number of grid points)
    unsigned int getSize() const;

    //! \brief set the occupancy for a multi-layered sphere
    /*!
      This function encodes the occupancy for a sphere and multiple layers around it
      \param center              location of the sphere center
      \param radius              Radius of the base sphere
      \param stepSize            thickness of each layer on top of the base sphere
      \param maxNumLayers        maximum number of layers, if -1 this is determined by
                                 the number of bits used per grid points in the storage
      \param ignoreOutOfBound    if true, ignore if center is outside the grid, otherwise throw
                                 an exception
      
    */
    void setSphereOccupancy(const Point3D & center, double radius, 
                            double stepSize, int maxNumLayers=-1, 
                            bool ignoreOutOfBound=true);

    const RDKit::DiscreteValueVect *getStorage() const;

    //! \brief get the index of the grid point given the x, y, z indices
    //!
    //! \return the integer value, -1 if the indices are outside the grid
    int getGridIndex(unsigned int xi, unsigned int yi, unsigned int zi) const;

    //! \brief get the number of grid points along x-axis
    unsigned int getNumX() const;

    //! \brief get the number of grid points along y-axis
    unsigned int getNumY() const;

    //! \brief get the number of grid points along z-axis
    unsigned int getNumZ() const;

    //! \brief get the grid's offset
    const Point3D &getOffset() const;

    //! \brief get the grid's spacing
    double getSpacing() const;

    //! \brief return a \b const pointer to our occupancy vector
    const RDKit::DiscreteValueVect *getOccupancyVect() const;

    //! \brief returns true if the grid \c other has parameters
    //!        compatible with ours.
    virtual bool compareParams(const UniformGrid3D &other) const;

    //! \brief create and return a pickle
    std::string toString() const;

  private:
    //! \brief internal initialization code
    /*
        \param dimX:    the x dimension of the grid, in Angstroms
        \param dimY:    the y dimension of the grid, in Angstroms
        \param dimZ:    the z dimension of the grid, in Angstroms
        \param spacing: the grid spacing, in Angstroms
        \param valType: the data type of the grid (determines the number of bits
                        per point)
        \param offset:  the offset of the grid from (0,0,0), in Angstroms.
      
    */
    void initGrid(double dimX, double dimY, double dimZ, double spacing,
                  RDKit::DiscreteValueVect::DiscreteValueType valType,
                  const RDGeom::Point3D &offSet);
    unsigned int d_numX, d_numY, d_numZ; //! number of grid points along x, y, z axes
    double d_spacing; //! grid spacing
    Point3D d_offSet; //! the grid offset (from the origin)
    RDKit::DiscreteValueVect *dp_storage; //! storage for values at each grid point

    //! \brief construct from a pickle
    void initFromText(const char *pkl,const unsigned int length);

  };

  //! \brief writes the contents of the grid to a stream
  /*
    The grid is written in GRD format
  */
  void writeGridToStream(const UniformGrid3D &grid, std::ostream &outStrm);

  //! \brief writes the contents of the grid to a named file
  /*
    The grid is written in GRD format
  */
  void writeGridToFile(const UniformGrid3D &grid, std::string filename);

}

#endif

