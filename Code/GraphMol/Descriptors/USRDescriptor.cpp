//
//  Copyright (C) 2011-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//


#include <Geometry/point.h>
#include <Numerics/Vector.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <boost/foreach.hpp>
#include "USRDescriptor.h"


namespace RDKit{

  namespace {
    void calcDistances(const RDGeom::Point3DConstPtrVect &coords,
                       const RDGeom::Point3D &point,
                       RDNumeric::Vector<double> &distances) {
      PRECONDITION(distances.size() == coords.size(), "distances and coords must be the same size");
      RDGeom::Point3D tmpPt;
      unsigned int i = 0;
      // loop over coordinates
      BOOST_FOREACH(const RDGeom::Point3D *tpp, coords) {
        distances[i++] = (*tpp-point).length();
      }
    }

    void calcCentroid(const RDGeom::Point3DConstPtrVect &coords,
                      RDGeom::Point3D &pt) {
      PRECONDITION(coords.size() != 0,"no coordinates");
      // set pt to zero
      pt *= 0.0;
      // loop over coordinates
      BOOST_FOREACH(const RDGeom::Point3D *opt, coords) {
        pt += *opt;
      }
      pt /= coords.size();
    }

    void calcMoments(const RDNumeric::Vector<double> &dist,
                     RDNumeric::Vector<double> &moments) {
      PRECONDITION(moments.size() == 3, "moments must have 3 elements");
      PRECONDITION(dist.size() != 0,"no distances");
      // set all elements to zero
      moments *= 0.0;
      unsigned int numPts = dist.size();
      // 1. moment: mean
      for (unsigned int i = 0; i < numPts; ++i) {
        moments[0] += dist[i];
      }
      moments[0] /= numPts;
      // 2. moment: standard deviation
      // 3. moment: cubic root of skewness
      for (unsigned int i = 0; i < numPts; ++i) {
        double diff = dist[i] - moments[0];
        moments[1] += diff * diff;
        moments[2] += diff * diff * diff;
      }
      moments[1] = sqrt(moments[1] / numPts);
      moments[2] /= numPts;
      moments[2] = cbrt(moments[2] / (moments[1] * moments[1] * moments[1]));
    }

    void addMoments(std::vector<double> &d,
                    const RDNumeric::Vector<double> &m,
                    int idx) {
      for (unsigned int i = 0; i < m.size(); ++i) {
        d[i + idx] = m[i];
      }
    }

  } // end namespace

  namespace Descriptors {

    void USR(const ROMol &mol, std::vector<double> &descriptor, int confId) {
      PRECONDITION(descriptor.size() == 12, "descriptor must have 12 elements");
      unsigned int na = mol.getNumAtoms();
      // check that number of atoms > 3
      if (na < 3) {
        throw ValueErrorException("Number of atoms must be greater than 3");
      }
      // check that minimum a conformer exists
      if (mol.getNumConformers() == 0) {
        throw ConformerException("No conformations available on this molecule");
      }

      const Conformer &conf = mol.getConformer(confId);
      RDGeom::Point3DConstPtrVect coords(na);
      // loop over atoms
      for (unsigned int ai = 0; ai < na; ++ai) {
        coords[ai] = &conf.getAtomPos(ai);
      }
      calcUSRForPoints(coords, descriptor);
    }

    void calcUSRForPoints(const RDGeom::Point3DConstPtrVect &coords,
                          std::vector<double> &descriptor) {
      PRECONDITION(descriptor.size() == 12, "descriptor must have 12 elements");
      RDGeom::Point3D pt;
      RDGeom::Point3D pt2;
      RDNumeric::Vector<double> dist(coords.size());
      RDNumeric::Vector<double> moments(3, 0.0);

      // ctd = centroid
      calcCentroid(coords, pt);
      calcDistances(coords, pt, dist);
      calcMoments(dist, moments);
      addMoments(descriptor, moments, 0);

      // catc = closest atom to centroid
      pt = (*coords[dist.smallestValId()]); // catc
      pt2 = (*coords[dist.largestValId()]); // fatc
      calcDistances(coords, pt, dist);
      calcMoments(dist, moments);
      addMoments(descriptor, moments, 3);

      // fatc = farthest atom to centroid
      calcDistances(coords, pt2, dist);
      calcMoments(dist, moments);
      addMoments(descriptor, moments, 6);

      // fatf = farthest atom to fatc
      pt = (*coords[dist.largestValId()]);
      calcDistances(coords, pt, dist);
      calcMoments(dist, moments);
      addMoments(descriptor, moments, 9);
    }

    double calcUSRScore(const std::vector<double> &d1, const std::vector<double> &d2) {
      PRECONDITION(d1.size() == d2.size(), "descriptors must have the same size");
      double score = 0.0;
      unsigned int num = d1.size();
      for (unsigned int i = 0; i < num; ++i) {
        score += fabs(d1[i] - d2[i]);
      }
      score /= num;
      return 1.0 / (1.0 + score);
    }

  } // end of namespace Descriptors
} //end of namespace RDKit
