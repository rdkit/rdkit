//
//  Copyright (C) 2019-2023 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <vector>
#include <list>
#include <unordered_map>
#include <Geometry/point.h>
#include <cmath>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace conrec {
struct ConrecSegment {
  RDGeom::Point2D p1;
  RDGeom::Point2D p2;
  double isoVal;
  ConrecSegment(double x1, double y1, double x2, double y2, double isoVal)
      : p1(x1, y1), p2(x2, y2), isoVal(isoVal) {}
  ConrecSegment(const RDGeom::Point2D &p1, const RDGeom::Point2D &p2,
                double isoVal)
      : p1(p1), p2(p2), isoVal(isoVal) {}
};
// adapted from conrec.c by Paul Bourke:
// http://paulbourke.net/papers/conrec/conrec.c
/*
   Derivation from the fortran version of CONREC by Paul Bourke
   d               ! matrix of data to contour
   ilb,iub,jlb,jub ! index bounds of data matrix
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   nc              ! number of contour levels
   z               ! contour levels in increasing order
*/
inline void Contour(const double *d, size_t ilb, size_t iub, size_t jlb,
                    size_t jub, const double *x, const double *y, size_t nc,
                    double *z, std::vector<ConrecSegment> &res) {
  PRECONDITION(d, "no data");
  PRECONDITION(x, "no data");
  PRECONDITION(y, "no data");
  PRECONDITION(z, "no data");
  PRECONDITION(nc > 0, "no contours");
  PRECONDITION(iub > ilb, "bad bounds");
  PRECONDITION(jub > jlb, "bad bounds");

  int m1, m2, m3, case_value;
  double dmin, dmax, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  int i, j, m;
  size_t k;
  double h[5];
  int sh[5];
  double xh[5], yh[5];
  int im[4] = {0, 1, 1, 0}, jm[4] = {0, 0, 1, 1};
  int castab[3][3][3] = {{{0, 0, 8}, {0, 2, 5}, {7, 6, 9}},
                         {{0, 3, 4}, {1, 3, 1}, {4, 3, 0}},
                         {{9, 6, 7}, {5, 2, 0}, {8, 0, 0}}};
  double temp1, temp2;
  size_t ny = jub - jlb + 1;

  auto xsect = [&](int p1, int p2) {
    return (h[p2] * xh[p1] - h[p1] * xh[p2]) / (h[p2] - h[p1]);
  };
  auto ysect = [&](int p1, int p2) {
    return (h[p2] * yh[p1] - h[p1] * yh[p2]) / (h[p2] - h[p1]);
  };

  for (j = (jub - 1); j >= (int)jlb; j--) {
    for (i = ilb; i <= (int)iub - 1; i++) {
      temp1 = std::min(d[i * ny + j], d[i * ny + j + 1]);
      temp2 = std::min(d[(i + 1) * ny + j], d[(i + 1) * ny + j + 1]);
      dmin = std::min(temp1, temp2);
      temp1 = std::max(d[i * ny + j], d[i * ny + j + 1]);
      temp2 = std::max(d[(i + 1) * ny + j], d[(i + 1) * ny + j + 1]);
      dmax = std::max(temp1, temp2);
      if (dmax < z[0] || dmin > z[nc - 1]) {
        continue;
      }
      for (k = 0; k < nc; k++) {
        if (z[k] < dmin || z[k] > dmax) {
          continue;
        }
        for (m = 4; m >= 0; m--) {
          if (m > 0) {
            h[m] = d[(i + im[m - 1]) * ny + j + jm[m - 1]] - z[k];
            xh[m] = x[i + im[m - 1]];
            yh[m] = y[j + jm[m - 1]];
          } else {
            h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
            xh[0] = 0.50 * (x[i] + x[i + 1]);
            yh[0] = 0.50 * (y[j] + y[j + 1]);
          }
          if (h[m] > 0.0) {
            sh[m] = 1;
          } else if (h[m] < 0.0) {
            sh[m] = -1;
          } else {
            sh[m] = 0;
          }
        }

        /*
           Note: at this stage the relative heights of the corners and the
           centre are in the h array, and the corresponding coordinates are
           in the xh and yh arrays. The centre of the box is indexed by 0
           and the 4 corners by 1 to 4 as shown below.
           Each triangle is then indexed by the parameter m, and the 3
           vertices of each triangle are indexed by parameters m1,m2,and m3.
           It is assumed that the centre of the box is always vertex 2
           though this isimportant only when all 3 vertices lie exactly on
           the same contour level, in which case only the side of the box
           is drawn.
              vertex 4 +-------------------+ vertex 3
                       | \               / |
                       |   \    m-3    /   |
                       |     \       /     |
                       |       \   /       |
                       |  m=2    X   m=2   |       the centre is vertex 0
                       |       /   \       |
                       |     /       \     |
                       |   /    m=1    \   |
                       | /               \ |
              vertex 1 +-------------------+ vertex 2
        */
        /* Scan each triangle in the box */
        for (m = 1; m <= 4; m++) {
          m1 = m;
          m2 = 0;
          if (m != 4) {
            m3 = m + 1;
          } else {
            m3 = 1;
          }
          if ((case_value = castab[sh[m1] + 1][sh[m2] + 1][sh[m3] + 1]) == 0) {
            continue;
          }
          switch (case_value) {
            case 1: /* Line between vertices 1 and 2 */
              x1 = xh[m1];
              y1 = yh[m1];
              x2 = xh[m2];
              y2 = yh[m2];
              break;
            case 2: /* Line between vertices 2 and 3 */
              x1 = xh[m2];
              y1 = yh[m2];
              x2 = xh[m3];
              y2 = yh[m3];
              break;
            case 3: /* Line between vertices 3 and 1 */
              x1 = xh[m3];
              y1 = yh[m3];
              x2 = xh[m1];
              y2 = yh[m1];
              break;
            case 4: /* Line between vertex 1 and side 2-3 */
              x1 = xh[m1];
              y1 = yh[m1];
              x2 = xsect(m2, m3);
              y2 = ysect(m2, m3);
              break;
            case 5: /* Line between vertex 2 and side 3-1 */
              x1 = xh[m2];
              y1 = yh[m2];
              x2 = xsect(m3, m1);
              y2 = ysect(m3, m1);
              break;
            case 6: /* Line between vertex 3 and side 1-2 */
              x1 = xh[m3];
              y1 = yh[m3];
              x2 = xsect(m1, m2);
              y2 = ysect(m1, m2);
              break;
            case 7: /* Line between sides 1-2 and 2-3 */
              x1 = xsect(m1, m2);
              y1 = ysect(m1, m2);
              x2 = xsect(m2, m3);
              y2 = ysect(m2, m3);
              break;
            case 8: /* Line between sides 2-3 and 3-1 */
              x1 = xsect(m2, m3);
              y1 = ysect(m2, m3);
              x2 = xsect(m3, m1);
              y2 = ysect(m3, m1);
              break;
            case 9: /* Line between sides 3-1 and 1-2 */
              x1 = xsect(m3, m1);
              y1 = ysect(m3, m1);
              x2 = xsect(m1, m2);
              y2 = ysect(m1, m2);
              break;
            default:
              break;
          }

          /* Finally draw the line */
          res.emplace_back(RDGeom::Point2D(x1, y1), RDGeom::Point2D(x2, y2),
                           z[k]);
        } /* m */
      } /* k - contour */
    } /* i */
  } /* j */
}

struct tplHash {
  template <class T1, class T2, class T3>
  std::size_t operator()(const std::tuple<T1, T2, T3> &p) const {
    std::size_t res = 0;
    boost::hash_combine(res, std::get<0>(p));
    boost::hash_combine(res, std::get<1>(p));
    boost::hash_combine(res, std::get<2>(p));
    return res;
  }
};

inline std::vector<std::pair<std::vector<RDGeom::Point2D>, double>>
connectLineSegments(const std::vector<ConrecSegment> &segments,
                    double coordMultiplier = 1000,
                    double isoValMultiplier = 1e6) {
  std::vector<std::pair<std::vector<RDGeom::Point2D>, double>> res;
  std::unordered_map<std::tuple<int, int, long>, std::list<size_t>, tplHash>
      endPointHashes;

  auto makePointKey = [&coordMultiplier, &isoValMultiplier](const auto &pt,
                                                            double isoVal) {
    return std::make_tuple<int, int, long>(
        std::lround(coordMultiplier * pt.x),
        std::lround(coordMultiplier * pt.y),
        std::lround(isoValMultiplier * isoVal));
  };

  // first hash all of the endpoints
  for (auto i = 0u; i < segments.size(); ++i) {
    const auto &seg = segments[i];
    endPointHashes[makePointKey(seg.p1, seg.isoVal)].push_back(i);
    endPointHashes[makePointKey(seg.p2, seg.isoVal)].push_back(i);
  }

  boost::dynamic_bitset<> segmentsDone(segments.size());
  // a candidate end point hasn't been used already.
  auto isCandidate = [&segmentsDone](const auto &pr) {
    return !pr.second.empty() && !segmentsDone[pr.second.front()];
  };
  auto singlePoint =
      std::find_if(endPointHashes.begin(), endPointHashes.end(), isCandidate);
  while (singlePoint != endPointHashes.end()) {
    auto segId = singlePoint->second.front();
    auto currKey = singlePoint->first;
    auto currVal = segments[segId].isoVal;
    std::vector<RDGeom::Point2D> contour;
    while (1) {
      segmentsDone.set(segId, true);
      // move onto the next segment
      const auto seg = segments[segId];
      auto k1 = makePointKey(seg.p1, seg.isoVal);
      auto k2 = makePointKey(seg.p2, seg.isoVal);
      auto endPtKey = k2;
      if (k1 == currKey) {
        contour.push_back(seg.p1);
      } else if (k2 == currKey) {
        contour.push_back(seg.p2);
        endPtKey = k1;
      }
      // remove this segment from the two hash lists:
      auto &segs1 = endPointHashes[currKey];
      segs1.erase(std::find(segs1.begin(), segs1.end(), segId));

      auto &segs = endPointHashes[endPtKey];
      segs.erase(std::find(segs.begin(), segs.end(), segId));
      if (segs.empty()) {
        // we're at the end, push on the last point
        if (k1 == currKey) {
          contour.push_back(seg.p2);
        } else if (k2 == currKey) {
          contour.push_back(seg.p1);
        }
        break;
      }
      segId = segs.front();
      currKey = endPtKey;
    }
    res.push_back(std::make_pair(contour, currVal));
    singlePoint = std::find_if(singlePoint, endPointHashes.end(), isCandidate);
  }
  return res;
}

}  // namespace conrec
