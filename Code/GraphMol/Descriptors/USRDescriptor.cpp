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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/foreach.hpp>
#include "USRDescriptor.h"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace RDKit{

  namespace {
    void calcDistances(const RDGeom::Point3DConstPtrVect &coords,
                       const RDGeom::Point3D &point,
                       std::vector<double> &distances) {
      distances.resize(coords.size());
      unsigned int i = 0;
      // loop over coordinates
      BOOST_FOREACH(const RDGeom::Point3D *tpp, coords) {
        distances[i++] = (*tpp-point).length();
      }
    }

    void calcCentroid(const RDGeom::Point3DConstPtrVect &coords,
                      RDGeom::Point3D &pt) {
      PRECONDITION(!coords.empty(),"no coordinates");
      // set pt to zero
      pt *= 0.0;
      // loop over coordinates
      BOOST_FOREACH(const RDGeom::Point3D *opt, coords) {
        pt += *opt;
      }
      pt /= coords.size();
    }

    unsigned int largestValId(const std::vector<double> &v) {
      PRECONDITION(!v.empty(),"no values");
      double res = v[0];
      unsigned int id = 0;
      for (unsigned int i = 1; i < v.size(); ++i) {
        if (v[i] > res) {
          res = v[i];
          id = i;
        }
      }
      return id;
    }

    unsigned int smallestValId(const std::vector<double> &v) {
      PRECONDITION(!v.empty(),"no values");
      double res = v[0];
      unsigned int id = 0;
      for (unsigned int i = 1; i < v.size(); ++i) {
        if (v[i] < res) {
          res = v[i];
          id = i;
        }
      }
      return id;
    }

    void calcMoments(const std::vector<double> &dist,
                     std::vector<double> &descriptor,
                     int idx) {
      std::vector<double> moments (3, 0.0);
      unsigned int numPts = dist.size();
      if (numPts > 0) {
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
      // add moments to descriptor
      std::copy(moments.begin(), moments.end(), descriptor.begin()+idx);
    }

    class ss_matcher {
      public:
        ss_matcher() {};
        ss_matcher(const std::string &pattern){
          RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
          TEST_ASSERT(p);
          m_matcher.reset(p);
        };
        const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };
      private:
        RDKit::ROMOL_SPTR m_matcher;
    };

    /*
    // Definitions for feature points adapted from:
    // Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
    const char *smartsPatterns[4] = {
      "[$([C;H2,H1](!=*)[C;H2,H1][C;H2,H1][$([C;H1,H2,H3]);!$(C=*)]),\
$(C([C;H2,H3])([C;H2,H3])[C;H2,H3])]", // Hydrophobic
      "[a]", //Aromatic
      "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]", // Donor
      "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]" // Acceptor
    };*/
    // Definitions for feature points from
    // http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd2?at=default
    const char *smartsPatterns[4] = {
      "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
      "[a]", // aromatic
      "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),\
$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
      "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]" // donor
    };
    std::vector<std::string> featureSmarts(smartsPatterns,smartsPatterns+4);
    typedef boost::flyweight<boost::flyweights::key_value<std::string,ss_matcher>,boost::flyweights::no_tracking > pattern_flyweight;

    void getAtomIdsForFeatures(const ROMol &mol, std::vector<std::vector<unsigned int> > &atomIds) {
      unsigned int numFeatures = featureSmarts.size();
      PRECONDITION(atomIds.size() == numFeatures, "atomIds must have 4 elements");
      std::vector<const ROMol *> featureMatchers;
      featureMatchers.reserve(numFeatures);
      BOOST_FOREACH(std::string feature, featureSmarts) {
        const ROMol *matcher = pattern_flyweight(feature).get().getMatcher();
        featureMatchers.push_back(matcher);
      }
      for (unsigned int i = 0; i < numFeatures; ++i) {
        std::vector<MatchVectType> matchVect;
        // to maintain thread safety, we have to copy the pattern molecules:
        SubstructMatch(mol,ROMol(*featureMatchers[i],true),matchVect);
        BOOST_FOREACH(MatchVectType mv, matchVect) {
          for (MatchVectType::const_iterator mIt = mv.begin(); mIt != mv.end(); ++mIt){
            atomIds[i].push_back(mIt->second);
          }
        }
      } // end loop over features
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
      // the four distances
      std::vector<std::vector<double> > dist(4);
      std::vector<RDGeom::Point3D> points(4);
      calcUSRDistributions(coords, dist, points);

      calcUSRFromDistributions(dist, descriptor);
    }

    void USRCAT(const ROMol &mol, std::vector<double> &descriptor,
                std::vector<std::vector<unsigned int> > &atomIds, int confId) {
      unsigned int na = mol.getNumAtoms();
      // check that number of atoms > 3
      if (na < 3) {
        throw ValueErrorException("Number of atoms must be greater than 3");
      }
      // check that minimum a conformer exists
      if (mol.getNumConformers() == 0) {
        throw ConformerException("No conformations available on this molecule");
      }

      // get atom selections
      unsigned int numClasses = atomIds.size();
      if (numClasses > 0) { // user provided atom selections
        PRECONDITION(descriptor.size() == (numClasses+1)*12, "descriptor must have (numClasses+1)*12 elements");
      } else { // use feature definitions of FeatureMorgan fingerprint
        numClasses = 4;
        PRECONDITION(descriptor.size() == 60, "descriptor must have 60 elements");
        atomIds.resize(numClasses);
        getAtomIdsForFeatures(mol, atomIds);
      }

      const Conformer &conf = mol.getConformer(confId);
      RDGeom::Point3DConstPtrVect coords(na);
      // loop over atoms
      for (unsigned int ai = 0; ai < na; ++ai) {
        coords[ai] = &conf.getAtomPos(ai);
      }
      // the original USR
      std::vector<std::vector<double> > dist(4);
      std::vector<RDGeom::Point3D> points(4);
      calcUSRDistributions(coords, dist, points);
      std::vector<double> tmpDescriptor(12);
      calcUSRFromDistributions(dist, tmpDescriptor);
      std::copy(tmpDescriptor.begin(), tmpDescriptor.end(), descriptor.begin());

      // loop over the atom selections
      unsigned int featIdx = 12;
      BOOST_FOREACH(std::vector<unsigned int> atoms, atomIds) {
        // reduce the coordinates to the atoms of interest
        RDGeom::Point3DConstPtrVect reducedCoords(atoms.size());
        unsigned int i = 0;
        BOOST_FOREACH(unsigned int idx, atoms) {
          reducedCoords[i++] = coords[idx];
        }
        calcUSRDistributionsFromPoints(reducedCoords, points, dist);
        calcUSRFromDistributions(dist, tmpDescriptor);
        std::copy(tmpDescriptor.begin(), tmpDescriptor.end(), descriptor.begin()+featIdx);
        featIdx += 12;
      }
    }

    void calcUSRDistributions(const RDGeom::Point3DConstPtrVect &coords,
                              std::vector<std::vector<double> > &dist,
                              std::vector<RDGeom::Point3D> &points) {
      PRECONDITION(dist.size() == 4, "dist must have 4 elements");
      PRECONDITION(points.size() == 4, "points must have 4 elements");
      // ctd = centroid
      calcCentroid(coords, points[0]);
      calcDistances(coords, points[0], dist[0]);
      // catc = closest atom to centroid
      points[1] = (*coords[smallestValId(dist[0])]);
      calcDistances(coords, points[1], dist[1]);
      // fatc = farthest atom to centroid
      points[2] = (*coords[largestValId(dist[0])]);
      calcDistances(coords, points[2], dist[2]);
      // fatf = farthest atom to fatc
      points[3] = (*coords[largestValId(dist[2])]);
      calcDistances(coords, points[3], dist[3]);
    }

    void calcUSRDistributionsFromPoints(const RDGeom::Point3DConstPtrVect &coords,
                                        const std::vector<RDGeom::Point3D> &points,
                                        std::vector<std::vector<double> > &dist) {
      PRECONDITION(points.size() == dist.size(), "points and dist must have the same size");
      for (unsigned int i = 0; i < points.size(); ++i) {
         calcDistances(coords, points[i], dist[i]);
      }
    }

    void calcUSRFromDistributions(const std::vector<std::vector<double> > &dist,
                                  std::vector<double> &descriptor) {
      PRECONDITION(descriptor.size() == 3*dist.size(), "descriptor must have 3 times more elements than dist");
      for (unsigned int i = 0; i < dist.size(); ++i) {
        calcMoments(dist[i], descriptor, 3*i);
      }
    }

    double calcUSRScore(const std::vector<double> &d1, const std::vector<double> &d2,
                           const std::vector<double> &weights) {
      unsigned int num = 12; // length of each subset
      PRECONDITION(d1.size() == d2.size(), "descriptors must have the same size");
      PRECONDITION(weights.size() == (d1.size()/num), "size of weights not correct");
      double score = 1.0;
      for (unsigned int w; w < weights.size(); ++w) {
        double tmpScore = 0.0;
        unsigned int offset = num*w;
        for (unsigned int i = 0; i < num; ++i) {
            tmpScore += fabs(d1[i+offset] - d2[i+offset]);
        }
        tmpScore /= num;
        score += weights[w]*tmpScore;
      }
      return 1.0 / score;
    }

  } // end of namespace Descriptors
} //end of namespace RDKit
