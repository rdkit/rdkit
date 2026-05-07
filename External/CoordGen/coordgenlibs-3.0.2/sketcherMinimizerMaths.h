/*
 *  sketcherMinimizerMaths.h
 *
 *
 *  Created by Nicola Zonta on 15/3/2011.
 *  Copyright Schrodinger, LLC. All rights reserved
 *
 */

#ifndef sketcherMINIMIZERMATHS_H
#define sketcherMINIMIZERMATHS_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>

#define MACROCYCLE 9 // smallest MACROCYCLE

#define BONDLENGTH 50
#define FRACTION_OF_BONDLENGTH_FOR_CLASH 0.25

#define SKETCHER_EPSILON 0.0001f

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

inline float roundToTwoDecimalDigits(float f)
{
    return static_cast<float>(floor(f * 100 + 0.5) * 0.01);
}

inline float roundToPrecision(float f, int precision)
{
    return static_cast<float>(floor(f * pow(10.f, precision) + 0.5) *
                              pow(0.1, precision));
}

/* class to represent a point or vector in 2d */
class sketcherMinimizerPointF
{
  public:
    sketcherMinimizerPointF() = default;

    sketcherMinimizerPointF(const sketcherMinimizerPointF& p)
        : xp(p.x()), yp(p.y())
    {
    }

    sketcherMinimizerPointF(float xpos, float ypos) : xp(xpos), yp(ypos) {}

    sketcherMinimizerPointF& operator=(const sketcherMinimizerPointF& p)
    {
        if (this != &p) {
            xp = p.x();
            yp = p.y();
        }
        return *this;
    }

    inline float x() const { return xp; }
    inline float y() const { return yp; }
    inline float& rx() { return xp; }
    inline float& ry() { return yp; }
    void setX(float x) { xp = x; }
    void setY(float y) { yp = y; }
    float squareLength() const
    {
        float dd = x() * x() + y() * y();
        return dd;
    }

    /* return the length of the vector */
    float length() const
    {
        float dd = squareLength();
        if (dd > SKETCHER_EPSILON) {
            return sqrt(dd);
        } else {
            return 0;
        }
    }

    /* normalize the vector */
    void normalize()
    {
        float q = length();
        if (q > SKETCHER_EPSILON) {
            xp /= q;
            yp /= q;
        }
    }

    /* rotate the vector by the angle with given sine and cosine */
    void rotate(float s, float c)
    {
        float x = xp;
        float y = yp;
        xp = x * c + y * s;
        yp = -x * s + y * c;
    }

    /* parallel component of this along the give axis */
    sketcherMinimizerPointF
    parallelComponent(const sketcherMinimizerPointF& axis)
    {
        float dotProduct = x() * axis.x() + y() * axis.y();
        return axis * dotProduct / axis.squareLength();
    }

    /* round the coordinates to the given number of decimal figures */
    void round(int precision = 2)
    {
        if (precision == 2) {
            xp = roundToTwoDecimalDigits(xp);
            yp = roundToTwoDecimalDigits(yp);
        } else {
            xp = roundToPrecision(xp, precision);
            yp = roundToPrecision(yp, precision);
        }
    }

    sketcherMinimizerPointF& operator+=(const sketcherMinimizerPointF& p)
    {
        xp += p.xp;
        yp += p.yp;
        return *this;
    }
    sketcherMinimizerPointF& operator-=(const sketcherMinimizerPointF& p)
    {
        xp -= p.xp;
        yp -= p.yp;
        return *this;
    }
    template <typename T> sketcherMinimizerPointF& operator*=(T c)
    {
        xp *= static_cast<float>(c);
        yp *= static_cast<float>(c);
        return *this;
    }
    template <typename T> sketcherMinimizerPointF& operator/=(T c)
    {
        xp /= static_cast<float>(c);
        yp /= static_cast<float>(c);
        return *this;
    }

    // friend inline bool operator==(const sketcherMinimizerPointF &p1, const
    // sketcherMinimizerPointF &p2) ;
    //    friend inline bool operator!=(const sketcherMinimizerPointF &, const
    //    sketcherMinimizerPointF &);

    inline friend std::ostream& operator<<(std::ostream& out,
                                           sketcherMinimizerPointF& point)
    {
        out << "(" << point.xp << ", " << point.yp << ")";
        return out;
    }

    friend inline const sketcherMinimizerPointF
    operator+(const sketcherMinimizerPointF& p1,
              const sketcherMinimizerPointF& p2)
    {
        return sketcherMinimizerPointF(p1.xp + p2.xp, p1.yp + p2.yp);
    }
    friend inline const sketcherMinimizerPointF
    operator-(const sketcherMinimizerPointF& p1,
              const sketcherMinimizerPointF& p2)
    {
        return sketcherMinimizerPointF(p1.xp - p2.xp, p1.yp - p2.yp);
    }
    friend inline const sketcherMinimizerPointF
    operator*(float c, const sketcherMinimizerPointF& p1)
    {
        return sketcherMinimizerPointF(p1.xp * c, p1.yp * c);
    }
    template <typename T>
    friend inline const sketcherMinimizerPointF
    operator*(const sketcherMinimizerPointF& p1, T c)
    {
        auto cf = static_cast<float>(c);
        return sketcherMinimizerPointF(p1.xp * cf, p1.yp * cf);
    }
    template <typename T>
    friend inline const sketcherMinimizerPointF
    operator/(const sketcherMinimizerPointF& p1, T c)
    {
        auto cf = static_cast<float>(c);
        return sketcherMinimizerPointF(p1.xp / cf, p1.yp / cf);
    }

    friend inline const sketcherMinimizerPointF
    operator-(const sketcherMinimizerPointF& p1)
    {
        return sketcherMinimizerPointF(-p1.xp, -p1.yp);
    }
    //     friend inline const sketcherMinimizerPointF operator/(const
    //     sketcherMinimizerPointF &, float);

  private:
    float xp{0.f};
    float yp{0.f};
};

/* return true if the two segments intersect and if a result pointer was given,
 * set it to the intersection point */
struct sketcherMinimizerMaths {
    static bool
    intersectionOfSegments(const sketcherMinimizerPointF& s1p1,
                           const sketcherMinimizerPointF& s1p2,
                           const sketcherMinimizerPointF& s2p1,
                           const sketcherMinimizerPointF& s2p2,
                           sketcherMinimizerPointF* result = nullptr)
    {
        /*
                Suppose the two line segments run from p to p + r and from q to
           q + s. Then any point on the first line is representable as p + t r
           (for a scalar parameter t) and any point on the second line as q + u
           s (for a scalar parameter u).

                The two lines intersect if we can find t and u such that:

                    p + t r = q + u s

                        Cross both sides with s, getting

                        (p + t r) × s = (q + u s) × s

                        And since s × s = 0, this means

                        t (r × s) = (q − p) × s

                        And therefore, solving for t:

                            t = (q − p) × s / (r × s)

                            In the same way, we can solve for u:

                                (p + t r) × r = (q + u s) × r

                                u (s × r) = (p − q) × r

                                u = (p − q) × r / (s × r)

                                To reduce the number of computation steps, it's
           convenient to rewrite this as follows (remembering that s × r = − r ×
           s):

                                u = (q − p) × r / (r × s)

                                Now there are five cases:

                                If r × s = 0 and (q − p) × r = 0, then the two
           lines are collinear. If in addition, either 0 ≤ (q − p) · r ≤ r · r
           or 0 ≤ (p − q) · s ≤ s · s, then the two lines are overlapping.

                                If r × s = 0 and (q − p) × r = 0, but neither 0
           ≤ (q − p) · r ≤ r · r nor 0 ≤ (p − q) · s ≤ s · s, then the two lines
           are collinear but disjoint.

                                If r × s = 0 and (q − p) × r ≠ 0, then the two
           lines are parallel and non-intersecting.

                                If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the
           two line segments meet at the point p + t r = q + u s.

                                Otherwise, the two line segments are not
           parallel but do not intersect.
        */

        const sketcherMinimizerPointF& p = s1p1;
        sketcherMinimizerPointF r = s1p2 - s1p1;
        const sketcherMinimizerPointF& q = s2p1;
        sketcherMinimizerPointF s = s2p2 - s2p1;
        float rxs = crossProduct(r, s);
        if (rxs > -SKETCHER_EPSILON &&
            rxs < SKETCHER_EPSILON) { // parallel lines
            return false;
        }
        sketcherMinimizerPointF qminusp = q - p;
        float t = crossProduct(qminusp, s) / rxs;
        if (t < 0 || t > 1) {
            return false;
        }
        float u = crossProduct(qminusp, r) / rxs;
        if (u < 0 || u > 1) {
            return false;
        }
        if (result) {
            *result = p + t * r;
        }
        return true;
    }

    /* signed angle between p1p2 and p2p3 */
    static float signedAngle(const sketcherMinimizerPointF& p1,
                             const sketcherMinimizerPointF& p2,
                             const sketcherMinimizerPointF& p3)
    {
        sketcherMinimizerPointF v1 = p1 - p2;
        sketcherMinimizerPointF v2 = p3 - p2;
        return float(atan2(v1.x() * v2.y() - v1.y() * v2.x(),
                           v1.x() * v2.x() + v1.y() * v2.y()) *
                     180 / M_PI);
    }

    /* unsigned angle between p1p2 and p2p3 */
    static float unsignedAngle(const sketcherMinimizerPointF& p1,
                               const sketcherMinimizerPointF& p2,
                               const sketcherMinimizerPointF& p3)
    {
        float x1 = p1.x();
        float y1 = p1.y();
        float x2 = p2.x();
        float y2 = p2.y();
        float x3 = p3.x();
        float y3 = p3.y();
        float v1x = x1 - x2;
        float v1y = y1 - y2;
        float v2x = x3 - x2;
        float v2y = y3 - y2;

        float d = sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y);
        if (d < SKETCHER_EPSILON) {
            d = SKETCHER_EPSILON;
        }
        float cosine = (v1x * v2x + v1y * v2y) / d;

        if (cosine < -1) {
            cosine = -1;
        } else if (cosine > 1) {
            cosine = 1;
        }
        return float((acos(cosine)) * 180 / M_PI);
    }

    /* return true if the two points are very close in space */
    static bool pointsCoincide(const sketcherMinimizerPointF& p1,
                               const sketcherMinimizerPointF& p2)
    {
        return ((p1 - p2).squareLength() < SKETCHER_EPSILON * SKETCHER_EPSILON);
    }

    /* return true if p1 and p2 are in the same semiplane defined by the given
     * segment */
    static bool sameSide(const sketcherMinimizerPointF& p1,
                         const sketcherMinimizerPointF& p2,
                         const sketcherMinimizerPointF& lineP1,
                         const sketcherMinimizerPointF& lineP2)
    {
        float x = lineP2.x() - lineP1.x();
        float y = lineP2.y() - lineP1.y();
        //    ///cerr << "("<<p1.x()<<","<<p1.y()<<") ("<<p2.x
        //    ()<<","<<p2.y()<<")  ("<<lineP1.x()<<","<<lineP1.y()<<")
        //    ("<<lineP2.x ()<<","<<lineP2.y()<<")"<<endl;
        if (fabs(float(x)) > fabs((y))) { // what about q?
            float m = y / x;

            float d1 = p1.y() - lineP1.y() - m * (p1.x() - lineP1.x());
            float d2 = p2.y() - lineP1.y() - m * (p2.x() - lineP1.x());
            return (d2 * d1 > 0);
        } else {
            float m = x / y;
            float d1 = p1.x() - lineP1.x() - m * (p1.y() - lineP1.y());
            float d2 = p2.x() - lineP1.x() - m * (p2.y() - lineP1.y());
            return (d2 * d1 > 0);
        }
    }

    /* return the projection of p on the line defined by the given segment */
    static sketcherMinimizerPointF
    projectPointOnLine(const sketcherMinimizerPointF& p,
                       const sketcherMinimizerPointF& sp1,
                       const sketcherMinimizerPointF& sp2)
    {
        sketcherMinimizerPointF l1 = p - sp1;
        sketcherMinimizerPointF l3 = sp2 - sp1;

        float segmentl2 = l3.squareLength();
        if (segmentl2 < SKETCHER_EPSILON) {
            segmentl2 = SKETCHER_EPSILON;
        }
        float t = sketcherMinimizerMaths::dotProduct(l1, l3) / segmentl2;
        return sp1 + t * l3;
    }

    /* squared distance of the given point from the given segment */
    static float squaredDistancePointSegment(const sketcherMinimizerPointF& p,
                                             const sketcherMinimizerPointF& sp1,
                                             const sketcherMinimizerPointF& sp2,
                                             float* returnT = nullptr)
    {

        sketcherMinimizerPointF l1 = p - sp1;
        sketcherMinimizerPointF l2 = sp2 - p;
        sketcherMinimizerPointF l3 = sp2 - sp1;

        float segmentl2 = l3.x() * l3.x() + l3.y() * l3.y();
        // float l1l =  sqrt ( l1.x () * l1.x () + l1.y() * l1.y() );

        if (segmentl2 < SKETCHER_EPSILON) {
            segmentl2 = SKETCHER_EPSILON;
        }
        float t = (l1.x() * l3.x() + l1.y() * l3.y()) / segmentl2;
        if (returnT != nullptr) {
            if (t < 0) {
                *returnT = 0;
            } else if (t > 1) {
                *returnT = 1;
            } else {
                *returnT = t;
            }
        }
        float squaredistance = 0.f;
        if (t < 0.f) {
            squaredistance = l1.x() * l1.x() + l1.y() * l1.y();
        } else if (t > 1.f) {

            squaredistance = l2.x() * l2.x() + l2.y() * l2.y();

        } else {

            sketcherMinimizerPointF proj = sp1 + t * l3;
            sketcherMinimizerPointF l5 = p - proj;
            squaredistance = l5.x() * l5.x() + l5.y() * l5.y();
        }
        if (squaredistance < SKETCHER_EPSILON) {
            squaredistance = SKETCHER_EPSILON;
        }

        return squaredistance;
    }

    static float squaredDistance(const sketcherMinimizerPointF& p1,
                                 const sketcherMinimizerPointF& p2)
    {
        return (p1.x() - p2.x()) * (p1.x() - p2.x()) +
               (p1.y() - p2.y()) * (p1.y() - p2.y());
    }

    static std::vector<float> tridiagonalSolve(const std::vector<float>& a,
                                               const std::vector<float>& b,
                                               const std::vector<float>& c,
                                               const std::vector<float>& rhs)
    {
        assert(a.size() == b.size() && a.size() == c.size() &&
               a.size() == rhs.size());
        assert(b[0] != 0.f);

        auto n = static_cast<unsigned int>(rhs.size());
        std::vector<float> u(n);
        std::vector<float> gam(n);

        float bet = b[0];
        u[0] = rhs[0] / bet;
        for (unsigned int j = 1; j < n; j++) {
            gam[j] = c[j - 1] / bet;
            bet = b[j] - a[j] * gam[j];
            assert(bet != 0.f);
            u[j] = (rhs[j] - a[j] * u[j - 1]) / bet;
        }
        for (unsigned int j = 1; j < n; j++) {
            u[n - j - 1] -= gam[n - j] * u[n - j];
        }

        return u;
    }

    /* used by ClosedBezierControlPoints */
    static std::vector<float> cyclicSolve(const std::vector<float>& a,
                                          const std::vector<float>& b,
                                          const std::vector<float>& c,
                                          float alpha, float beta,
                                          const std::vector<float>& rhs)
    {
        assert(a.size() == b.size() && a.size() == c.size());
        auto n = static_cast<unsigned int>(b.size());
        assert(n > 2);
        float gamma = -b[0]; // Avoid subtraction error in forming bb[0].
        // Set up the diagonal of the modified tridiagonal system.
        std::vector<float> bb(n);
        bb[0] = b[0] - gamma;
        bb[n - 1] = b[n - 1] - alpha * beta / gamma;
        for (unsigned int i = 1; i < n - 1; i++) {
            bb[i] = b[i];
        }
        // Solve A · x = rhs.
        std::vector<float> solution = tridiagonalSolve(a, bb, c, rhs);
        std::vector<float> x = solution;

        // Set up the vector u.
        std::vector<float> u(n);
        u[0] = gamma;
        u[n - 1] = alpha;
        for (unsigned int i = 1; i < n - 1; i++) {
            u[i] = 0.f;
        }
        // Solve A · z = u.
        solution = tridiagonalSolve(a, bb, c, u);
        std::vector<float> z = solution;

        // Form v · x/(1 + v · z).
        double fact = (x[0] + beta * x[n - 1] / gamma) /
                      (1.f + z[0] + beta * z[n - 1] / gamma);

        // Now get the solution vector x.
        for (unsigned int i = 0; i < n; i++) {
            x[i] -= float(fact * z[i]);
        }
        return x;
    }

    static sketcherMinimizerPointF
    pointOnCubicBezier(const sketcherMinimizerPointF& p1,
                       const sketcherMinimizerPointF& cp1,
                       const sketcherMinimizerPointF& cp2,
                       const sketcherMinimizerPointF& p2, float t)
    {
        // using Casteljiau's algorithm
        auto v1 = (1 - t) * p1 + t * cp1;
        auto v2 = (1 - t) * cp1 + t * cp2;
        auto v3 = (1 - t) * cp2 + t * p2;
        auto v4 = (1 - t) * v1 + t * v2;
        auto v5 = (1 - t) * v2 + t * v3;
        return (1 - t) * v4 + t * v5;
    }

    /* find control points to a closed Bezier curve that passes through the
     * given points */
    static void ClosedBezierControlPoints(
        const std::vector<sketcherMinimizerPointF>& knots,
        std::vector<sketcherMinimizerPointF>& firstControlPoints,
        std::vector<sketcherMinimizerPointF>& secondControlPoints)
    {

        auto n = static_cast<unsigned int>(knots.size());
        if (n <= 2) {
            return;
        }

        // Calculate first Bezier control points

        std::vector<float> a(n), b(n), c(n);
        for (unsigned int i = 0; i < n; i++) {
            a[i] = 1;
            b[i] = 4;
            c[i] = 1;
        }

        std::vector<float> rhs(n);
        for (unsigned int i = 0; i < n; i++) {
            int j = i + 1;
            if (j > int(n - 1)) {
                j = 0;
            }
            rhs[i] = 4 * knots[i].x() + 2 * knots[j].x();
        }
        // Solve the system for X.
        std::vector<float> x = cyclicSolve(a, b, c, 1, 1, rhs);

        for (unsigned int i = 0; i < n; i++) {
            int j = i + 1;
            if (j > int(n - 1)) {
                j = 0;
            }
            rhs[i] = 4 * knots[i].y() + 2 * knots[j].y();
        }
        // Solve the system for Y.
        std::vector<float> y = cyclicSolve(a, b, c, 1, 1, rhs);

        // Fill output arrays.
        firstControlPoints.resize(n);
        secondControlPoints.resize(n);
        for (unsigned int i = 0; i < n; i++) {
            firstControlPoints[i] = sketcherMinimizerPointF(x[i], y[i]);
            secondControlPoints[i] = sketcherMinimizerPointF(
                2 * knots[i].x() - x[i], 2 * knots[i].y() - y[i]);
        }
    }

    /* return the mirror image of the given point wrt the given segment */
    static sketcherMinimizerPointF
    mirrorPoint(const sketcherMinimizerPointF& point,
                const sketcherMinimizerPointF& segmentPoint1,
                const sketcherMinimizerPointF& segmentPoint2)
    {
        sketcherMinimizerPointF segmentV = segmentPoint2 - segmentPoint1;
        sketcherMinimizerPointF v2 = point - segmentPoint1;
        sketcherMinimizerPointF parallelComponent =
            v2.parallelComponent(segmentV);
        sketcherMinimizerPointF normalComponent = v2 - parallelComponent;
        return segmentPoint1 + parallelComponent - normalComponent;
    }

    /* dot product of two vectors */
    static float dotProduct(const sketcherMinimizerPointF& a,
                            const sketcherMinimizerPointF& b)
    {

        return (a.x() * b.x() + a.y() * b.y());
    }

    /* cross product of two vectors */
    static float crossProduct(const sketcherMinimizerPointF& a,
                              const sketcherMinimizerPointF& b)
    {

        return (a.x() * b.y() - a.y() * b.x());
    }

    static float cannonBallDistance(float originX, float originY, float originZ,
                                    float directionX, float directionY,
                                    float directionZ, float targetX,
                                    float targetY, float targetZ, float ballR,
                                    float targetR, float cutOff = 4.f)
    // how far can a cannonball of radius ballR shot from origin travel befor
    // hitting a target ball of targetR radius
    {

        // assume that direction is normalized

        float targetdX = targetX - originX;
        float targetdY = targetY - originY;
        float targetdZ = targetZ - originZ;
        float rR = ballR + targetR;
        float d2 = (targetdX * targetdX) + (targetdY * targetdY) +
                   (targetdZ * targetdZ);

        if (d2 > (cutOff + rR) * (cutOff + rR)) {
            return cutOff;
        }

        if (d2 < rR * rR) {
            return 0;
        }

        float d = sqrt(d2);
        if (d > SKETCHER_EPSILON) {
            targetdX /= d;
            targetdY /= d;
            targetdZ /= d;
        }
        float cos = targetdX * directionX + targetdY * directionY +
                    targetdZ * directionZ;
        if (cos < 0) {
            return cutOff;
        }
        float sin = sqrt(1 - (cos * cos));
        float f = d * sin;
        if (f > rR) {
            return cutOff;
        }
        float result = sqrt(d2 - (f * f)) - sqrt((rR * rR) - (f * f));
        if (result > cutOff) {
            return cutOff;
        }
        return result;
    }

    /* length of a 3d vector */
    static float length3D(float x, float y, float z)
    {
        float m = x * x + y * y + z * z;
        if (m > SKETCHER_EPSILON) {
            m = sqrt(m);
        }
        return m;
    }

    /* dot product of two 3d vectors */
    static float dotProduct3D(float x1, float y1, float z1, float x2, float y2,
                              float z2)
    {
        return x1 * x2 + y1 * y2 + z1 * z2;
    }

    /* cross product of two 3d vectors */
    static void crossProduct3D(float x1, float y1, float z1, float x2, float y2,
                               float z2, float& xr, float& yr, float& zr)

    {
        xr = y1 * z2 - z1 * y2;
        yr = z1 * x2 - x1 * z2;
        zr = x1 * y2 - y1 * x2;
    }
    static float distance3D(float x1, float y1, float z1, float x2, float y2,
                            float z2)
    {
        return length3D(x2 - x1, y2 - y1, z2 - z1);
    }

    /* angle between two 3d vectors */
    static float angle3D(float x1, float y1, float z1, float x2, float y2,
                         float z2, float x3, float y3, float z3)
    {
        float xa = x1 - x2;
        float ya = y1 - y2;
        float za = z1 - z2;
        float xb = x3 - x2;
        float yb = y3 - y2;
        float zb = z3 - z2;
        float l1 = length3D(xa, ya, za);
        float l2 = length3D(xb, yb, zb);
        float dp = dotProduct3D(xa, ya, za, xb, yb, zb);
        return static_cast<float>(acos(dp / (l1 * l2)) * 180.f / M_PI);
    }

    /* dihedral angle defined by 4 3d points */
    static float dihedral3D(float x1, float y1, float z1, float x2, float y2,
                            float z2, float x3, float y3, float z3, float x4,
                            float y4, float z4)
    {
        float xa, ya, za;
        crossProduct3D(x1 - x2, y1 - y2, z1 - z2, x3 - x2, y3 - y2, z3 - z2, xa,
                       ya, za);
        float xb, yb, zb;
        crossProduct3D(x2 - x3, y2 - y3, z2 - z3, x4 - x3, y4 - y3, z4 - z3, xb,
                       yb, zb);
        return angle3D(xa, ya, za, 0, 0, 0, xb, yb, zb);
    }

}; // struct sketcherMinimizerMaths

#endif // sketcherMINIMIZERMATHS_H
