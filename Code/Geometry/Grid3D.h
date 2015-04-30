// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _GRID3D_H_20050124_1113
#define _GRID3D_H_20050124_1113
#include <exception>
#include <string>

namespace RDKit {
  class DiscreteValueVect;
}
namespace RDGeom {
  class Point3D;
  
  class GridException : public std::exception {
  public:
    //! construct with an error message
    GridException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    GridException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~GridException () throw () {};
  private:
    std::string _msg;
  };
  
  //! Virtual base class for a grid object
  class Grid3D {
  public:
    virtual ~Grid3D() {};
    virtual int getGridPointIndex(const Point3D &point) const = 0;
    virtual int getVal(const Point3D &point) const = 0;
    virtual void setVal(const Point3D &point, unsigned int val) = 0;
  
    virtual Point3D getGridPointLoc(unsigned int pointId) const = 0;
    virtual unsigned int getVal(unsigned int pointId) const = 0;
    virtual void setVal(unsigned int pointId, unsigned int val) = 0;

    virtual unsigned int getSize() const = 0;

    virtual const RDKit::DiscreteValueVect *getOccupancyVect() const = 0;
  };
}

#endif 
