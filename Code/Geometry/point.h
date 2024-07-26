//
// Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef __RD_POINT_H__
#define __RD_POINT_H__
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <RDGeneral/Invariant.h>
#include <Numerics/Vector.h>
#include <boost/smart_ptr.hpp>

namespace RDGeom {

class RDKIT_RDGEOMETRYLIB_EXPORT Point {
  // this is the virtual base class, mandating certain functions
 public:
  virtual ~Point() {}

  virtual double operator[](unsigned int i) const = 0;
  virtual double &operator[](unsigned int i) = 0;

  virtual void normalize() = 0;
  virtual double length() const = 0;
  virtual double lengthSq() const = 0;
  virtual unsigned int dimension() const = 0;

  virtual Point *copy() const = 0;
};
#ifndef _MSC_VER
// g++ (at least as of v9.3.0) generates some spurious warnings from here.
// disable them
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// typedef class Point3D Point;
class RDKIT_RDGEOMETRYLIB_EXPORT Point3D : public Point {
 public:
  double x{0.0};
  double y{0.0};
  double z{0.0};

  Point3D() {}
  Point3D(double xv, double yv, double zv) : x(xv), y(yv), z(zv) {}

  ~Point3D() override = default;

  Point3D(const Point3D &other)
      : Point(other), x(other.x), y(other.y), z(other.z) {}

  Point *copy() const override { return new Point3D(*this); }

  inline unsigned int dimension() const override { return 3; }

  inline double operator[](unsigned int i) const override {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw ValueErrorException("Invalid index on Point3D");
        break;
    }
  }

  inline double &operator[](unsigned int i) override {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw ValueErrorException("Invalid index on Point3D");
        break;
    }
  }

  Point3D &operator=(const Point3D &other) {
    if (&other == this) {
      return *this;
    }
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  Point3D &operator+=(const Point3D &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Point3D &operator-=(const Point3D &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  Point3D &operator*=(double scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;
  }

  Point3D &operator/=(double scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
  }

  Point3D operator-() const {
    Point3D res(x, y, z);
    res.x *= -1.0;
    res.y *= -1.0;
    res.z *= -1.0;
    return res;
  }

  void normalize() override {
    double l = this->length();
    x /= l;
    y /= l;
    z /= l;
  }

  double length() const override {
    double res = x * x + y * y + z * z;
    return sqrt(res);
  }

  double lengthSq() const override {
    // double res = pow(x,2) + pow(y,2) + pow(z,2);
    double res = x * x + y * y + z * z;
    return res;
  }

  double dotProduct(const Point3D &other) const {
    double res = x * (other.x) + y * (other.y) + z * (other.z);
    return res;
  }

  /*! \brief determines the angle between a vector to this point
   *   from the origin and a vector to the other point.
   *
   *  The angle is unsigned: the results of this call will always
   *   be between 0 and M_PI
   */
  double angleTo(const Point3D &other) const {
    double lsq = lengthSq() * other.lengthSq();
    double dotProd = dotProduct(other);
    dotProd /= sqrt(lsq);

    // watch for roundoff error:
    if (dotProd <= -1.0) {
      return M_PI;
    }
    if (dotProd >= 1.0) {
      return 0.0;
    }

    return acos(dotProd);
  }

  /*! \brief determines the signed angle between a vector to this point
   *   from the origin and a vector to the other point.
   *
   *  The results of this call will be between 0 and M_2_PI
   */
  double signedAngleTo(const Point3D &other) const {
    double res = this->angleTo(other);
    // check the sign of the z component of the cross product:
    if ((this->x * other.y - this->y * other.x) < -1e-6) {
      res = 2.0 * M_PI - res;
    }
    return res;
  }

  /*! \brief Returns a normalized direction vector from this
   *   point to another.
   *
   */
  Point3D directionVector(const Point3D &other) const {
    Point3D res;
    res.x = other.x - x;
    res.y = other.y - y;
    res.z = other.z - z;
    res.normalize();
    return res;
  }

  /*! \brief Cross product of this point with the another point
   *
   * The order is important here
   *  The result is "this" cross with "other" not (other x this)
   */
  Point3D crossProduct(const Point3D &other) const {
    Point3D res;
    res.x = y * (other.z) - z * (other.y);
    res.y = -x * (other.z) + z * (other.x);
    res.z = x * (other.y) - y * (other.x);
    return res;
  }

  /*! \brief Get a unit perpendicular from this point (treating it as a vector):
   *
   */
  Point3D getPerpendicular() const {
    Point3D res(0.0, 0.0, 0.0);
    if (x) {
      if (y) {
        res.y = -1 * x;
        res.x = y;
      } else if (z) {
        res.z = -1 * x;
        res.x = z;
      } else {
        res.y = 1;
      }
    } else if (y) {
      if (z) {
        res.z = -1 * y;
        res.y = z;
      } else {
        res.x = 1;
      }
    } else if (z) {
      res.x = 1;
    }
    double l = res.length();
    POSTCONDITION(l > 0.0, "zero perpendicular");
    res /= l;
    return res;
  }
};

// given a  set of four pts in 3D compute the dihedral angle between the
// plane of the first three points (pt1, pt2, pt3) and the plane of the
// last three points (pt2, pt3, pt4)
// the computed angle is between 0 and PI
RDKIT_RDGEOMETRYLIB_EXPORT double computeDihedralAngle(const Point3D &pt1,
                                                       const Point3D &pt2,
                                                       const Point3D &pt3,
                                                       const Point3D &pt4);

// given a  set of four pts in 3D compute the signed dihedral angle between the
// plane of the first three points (pt1, pt2, pt3) and the plane of the
// last three points (pt2, pt3, pt4)
// the computed angle is between -PI and PI
RDKIT_RDGEOMETRYLIB_EXPORT double computeSignedDihedralAngle(
    const Point3D &pt1, const Point3D &pt2, const Point3D &pt3,
    const Point3D &pt4);

class RDKIT_RDGEOMETRYLIB_EXPORT Point2D : public Point {
 public:
  double x{0.0};
  double y{0.0};

  Point2D() {}
  Point2D(double xv, double yv) : x(xv), y(yv) {}
  ~Point2D() override = default;

  Point2D(const Point2D &other) : Point(other), x(other.x), y(other.y) {}
  //! construct from a Point3D (ignoring the z coordinate)
  Point2D(const Point3D &p3d) : Point(p3d), x(p3d.x), y(p3d.y) {}

  Point *copy() const override { return new Point2D(*this); }

  inline unsigned int dimension() const override { return 2; }

  inline double operator[](unsigned int i) const override {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      default:
        throw ValueErrorException("Invalid index on Point2D");
        break;
    }
  }

  inline double &operator[](unsigned int i) override {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      default:
        throw ValueErrorException("Invalid index on Point2D");
        break;
    }
  }

  Point2D &operator=(const Point2D &other) {
    x = other.x;
    y = other.y;
    return *this;
  }

  Point2D &operator+=(const Point2D &other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Point2D &operator-=(const Point2D &other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  Point2D &operator*=(double scale) {
    x *= scale;
    y *= scale;
    return *this;
  }

  Point2D &operator/=(double scale) {
    x /= scale;
    y /= scale;
    return *this;
  }

  Point2D operator-() const {
    Point2D res(x, y);
    res.x *= -1.0;
    res.y *= -1.0;
    return res;
  }

  void normalize() override {
    double ln = this->length();
    x /= ln;
    y /= ln;
  }

  void rotate90() {
    double temp = x;
    x = -y;
    y = temp;
  }

  double length() const override {
    // double res = pow(x,2) + pow(y,2);
    double res = x * x + y * y;
    return sqrt(res);
  }

  double lengthSq() const override {
    double res = x * x + y * y;
    return res;
  }

  double dotProduct(const Point2D &other) const {
    double res = x * (other.x) + y * (other.y);
    return res;
  }

  double angleTo(const Point2D &other) const {
    auto t1 = *this;
    auto t2 = other;
    t1.normalize();
    t2.normalize();
    double dotProd = t1.dotProduct(t2);
    // watch for roundoff error:
    if (dotProd < -1.0) {
      dotProd = -1.0;
    } else if (dotProd > 1.0) {
      dotProd = 1.0;
    }
    return acos(dotProd);
  }

  double signedAngleTo(const Point2D &other) const {
    double res = this->angleTo(other);
    if ((this->x * other.y - this->y * other.x) < -1e-6) {
      res = 2.0 * M_PI - res;
    }
    return res;
  }

  Point2D directionVector(const Point2D &other) const {
    Point2D res;
    res.x = other.x - x;
    res.y = other.y - y;
    res.normalize();
    return res;
  }
};

class RDKIT_RDGEOMETRYLIB_EXPORT PointND : public Point {
 public:
  typedef boost::shared_ptr<RDNumeric::Vector<double>> VECT_SH_PTR;

  PointND(unsigned int dim) {
    RDNumeric::Vector<double> *nvec = new RDNumeric::Vector<double>(dim, 0.0);
    dp_storage.reset(nvec);
  }

  PointND(const PointND &other) : Point(other) {
    RDNumeric::Vector<double> *nvec =
        new RDNumeric::Vector<double>(*other.getStorage());
    dp_storage.reset(nvec);
  }

  Point *copy() const override { return new PointND(*this); }

#if 0
	template <typename T>
    PointND(const T &vals){
      RDNumeric::Vector<double> *nvec = new RDNumeric::Vector<double>(vals.size(), 0.0);
      dp_storage.reset(nvec);
      unsigned int idx=0;
      typename T::const_iterator it;
      for(it=vals.begin();
          it!=vals.end();
          ++it){
        nvec->setVal(idx,*it);
        ++idx;
      };
    };
#endif

  ~PointND() override = default;

  inline double operator[](unsigned int i) const override {
    return dp_storage.get()->getVal(i);
  }

  inline double &operator[](unsigned int i) override {
    return (*dp_storage.get())[i];
  }

  inline void normalize() override { dp_storage.get()->normalize(); }

  inline double length() const override { return dp_storage.get()->normL2(); }

  inline double lengthSq() const override {
    return dp_storage.get()->normL2Sq();
  }

  unsigned int dimension() const override { return dp_storage.get()->size(); }

  PointND &operator=(const PointND &other) {
    if (this == &other) {
      return *this;
    }

    RDNumeric::Vector<double> *nvec =
        new RDNumeric::Vector<double>(*other.getStorage());
    dp_storage.reset(nvec);
    return *this;
  }

  PointND &operator+=(const PointND &other) {
    (*dp_storage.get()) += (*other.getStorage());
    return *this;
  }

  PointND &operator-=(const PointND &other) {
    (*dp_storage.get()) -= (*other.getStorage());
    return *this;
  }

  PointND &operator*=(double scale) {
    (*dp_storage.get()) *= scale;
    return *this;
  }

  PointND &operator/=(double scale) {
    (*dp_storage.get()) /= scale;
    return *this;
  }

  PointND directionVector(const PointND &other) {
    PRECONDITION(this->dimension() == other.dimension(),
                 "Point dimensions do not match");
    PointND np(other);
    np -= (*this);
    np.normalize();
    return np;
  }

  double dotProduct(const PointND &other) const {
    return dp_storage.get()->dotProduct(*other.getStorage());
  }

  double angleTo(const PointND &other) const {
    double dp = this->dotProduct(other);
    double n1 = this->length();
    double n2 = other.length();
    if ((n1 > 1.e-8) && (n2 > 1.e-8)) {
      dp /= (n1 * n2);
    }
    if (dp < -1.0) {
      dp = -1.0;
    } else if (dp > 1.0) {
      dp = 1.0;
    }
    return acos(dp);
  }

 private:
  VECT_SH_PTR dp_storage;
  inline const RDNumeric::Vector<double> *getStorage() const {
    return dp_storage.get();
  }
};
#ifndef _MSC_VER
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

typedef std::vector<RDGeom::Point *> PointPtrVect;
typedef PointPtrVect::iterator PointPtrVect_I;
typedef PointPtrVect::const_iterator PointPtrVect_CI;

typedef std::vector<RDGeom::Point3D *> Point3DPtrVect;
typedef std::vector<RDGeom::Point2D *> Point2DPtrVect;
typedef Point3DPtrVect::iterator Point3DPtrVect_I;
typedef Point3DPtrVect::const_iterator Point3DPtrVect_CI;
typedef Point2DPtrVect::iterator Point2DPtrVect_I;
typedef Point2DPtrVect::const_iterator Point2DPtrVect_CI;

typedef std::vector<const RDGeom::Point3D *> Point3DConstPtrVect;
typedef Point3DConstPtrVect::iterator Point3DConstPtrVect_I;
typedef Point3DConstPtrVect::const_iterator Point3DConstPtrVect_CI;

typedef std::vector<Point3D> POINT3D_VECT;
typedef std::vector<Point3D>::iterator POINT3D_VECT_I;
typedef std::vector<Point3D>::const_iterator POINT3D_VECT_CI;

typedef std::map<int, Point2D> INT_POINT2D_MAP;
typedef INT_POINT2D_MAP::iterator INT_POINT2D_MAP_I;
typedef INT_POINT2D_MAP::const_iterator INT_POINT2D_MAP_CI;

RDKIT_RDGEOMETRYLIB_EXPORT std::ostream &operator<<(std::ostream &target,
                                                    const RDGeom::Point &pt);

RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point3D operator+(const RDGeom::Point3D &p1,
                                                     const RDGeom::Point3D &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point3D operator-(const RDGeom::Point3D &p1,
                                                     const RDGeom::Point3D &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point3D operator*(const RDGeom::Point3D &p1,
                                                     double v);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point3D operator/(const RDGeom::Point3D &p1,
                                                     double v);

RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point2D operator+(const RDGeom::Point2D &p1,
                                                     const RDGeom::Point2D &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point2D operator-(const RDGeom::Point2D &p1,
                                                     const RDGeom::Point2D &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point2D operator*(const RDGeom::Point2D &p1,
                                                     double v);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point2D operator/(const RDGeom::Point2D &p1,
                                                     double v);

RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::PointND operator+(const RDGeom::PointND &p1,
                                                     const RDGeom::PointND &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::PointND operator-(const RDGeom::PointND &p1,
                                                     const RDGeom::PointND &p2);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::PointND operator*(const RDGeom::PointND &p1,
                                                     double v);
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::PointND operator/(const RDGeom::PointND &p1,
                                                     double v);
}  // namespace RDGeom

#endif
