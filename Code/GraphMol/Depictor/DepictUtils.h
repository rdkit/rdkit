//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_DEPICT_UTILS_H_
#define _RD_DEPICT_UTILS_H_

// REVIEW: remove extra headers here
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <Geometry/Transform2D.h>
#include <Geometry/point.h>
#include <queue>

namespace RDDepict {
  extern double BOND_LEN;
  extern double COLLISION_THRES;
  extern double BOND_THRES;
  extern double ANGLE_OPEN;
  extern unsigned int MAX_COLL_ITERS;
  extern double HETEROATOM_COLL_SCALE;
  extern unsigned int NUM_BONDS_FLIPS;

  typedef std::vector<const RDGeom::Point2D *> VECT_C_POINT;
  
  typedef std::pair<int, int> PAIR_I_I;
  typedef std::vector<PAIR_I_I> VECT_PII;
  struct gtIIPair {
    bool operator() ( const PAIR_I_I &pd1, const PAIR_I_I &pd2) const {
      return pd1.first > pd2.first;
    }
  };

  typedef std::priority_queue<PAIR_I_I, VECT_PII, gtIIPair> PR_QUEUE;

  typedef std::pair<double, PAIR_I_I> PAIR_D_I_I;
  typedef std::list<PAIR_D_I_I> LIST_PAIR_DII;


  //! Some utility functions used in generating 2D coordinates
  
  //! Embed a ring as a convex polygon in 2D
  /*!
    The process here is very straightforward:

    We take the center of the ring to lie at the origin, so put the first 
    point at the origin and then sweep
    anti-clockwise by an angle A = 360/n for the next point.

    The length of the arm (l) we want to sweep is easy to compute given the 
    bond length (b) we want to use for each bond in the ring (for now
    we will assume that this bond legnth is the same for all bonds in the ring:

    l = b/sqrt(2*(1 - cos(A)) 

    the above formula derives from the triangle formula, where side 'c' is given
    in terms of sides 'a' and 'b' as:

    c = a^2 + b^2 - 2.a.b.cos(A) 

    where A is the angle between a and b
   */
  RDGeom::INT_POINT2D_MAP embedRing(const RDKit::INT_VECT &ring);

  void transformPoints(RDGeom::INT_POINT2D_MAP &nringCor, const RDGeom::Transform2D &trans);

  //! Find a point that bisects the angle at rcr
  /*!
    The new point lies between nb1 and nb2. The line (rcr, newPt) bisects the angle
    'ang' at rcr
  */
  RDGeom::Point2D computeBisectPoint(const RDGeom::Point2D &rcr, 
                                     double ang, const RDGeom::Point2D &nb1,
                                     const RDGeom::Point2D &nb2);

  //! Reflect a set of point through a the line joining two point
  /*!
    ARGUMENTS:
    \param coordMap       a map of <int, point2D> going from atom id to current
                          coordinates of the points that need to be reflected:
                          The coordinates are overwritten
    \param loc1           the first point of the line that is to be used as a mirror
    \param loc2           the second point of the line to be used as a mirror
   */
  void reflectPoints(RDGeom::INT_POINT2D_MAP &coordMap, const RDGeom::Point2D &loc1,
                   const RDGeom::Point2D &loc2);

  RDGeom::Point2D reflectPoint(const RDGeom::Point2D &point, const RDGeom::Point2D &loc1,
                               const RDGeom::Point2D &loc2);

  //! Set the neighbors yet to added to aid such that the atoms with the most subs fall on opposite sides
  /*!
    Ok this needs some explanation
    - Let A, B, C, D be the substituent on the central atom X (given
      by atom index aid)
    - also let A be the atom that is already embedded
    - Now we want the order in which the remaining neighbors B,C,D are
      added to X such that the atoms with atom with largest number of
      substituent fall on opposite sides of X so as to minimize atom
      clashes later in the depiction

    E.g. let say we have the following situation
<pre>    
            B
         |  |  
         A--X--C
         |  |
          --D--
            |
</pre>    
    In this case the the number substituent of A, B, C, D are 3, 1, 1,
    4 respectively so want to A and D to go opposite sides and so that
    we draw
<pre>    
            B
         |  |  |
         A--X--D--
         |  |  |
            C
</pre>    
    And the correct ordering of the neighbors is B,D,C
  */
  RDKit::INT_VECT setNbrOrder(unsigned int aid, const RDKit::INT_VECT &nbrs, const RDKit::ROMol &mol);

  //! \brief From a given set of rings find the ring the largest common elements with other rings
  /*
    Bit of a weird function - this is typically called once we have embedded some of the 
    rings in a fused system and we are looking for the ring that must be embedded (or merged)
    next. The heuristic used here is to pick the rings with the maximum number of atoms
    in common with the rings that are already embedded.
    
    \param doneRings    a vertor of ring IDs that have been embedded already
    \param fusedRings   list of all the rings in the the fused system
    \param nextId       this is where the ID for the next ring is written
    
    \return list of atom ids that are common
  */
  RDKit::INT_VECT findNextRingToEmbed(const RDKit::INT_VECT &doneRings, 
                                      const RDKit::VECT_INT_VECT &fusedRings, 
                                      int &nextId);

  typedef std::pair<int,int> INT_PAIR;
  typedef std::vector<INT_PAIR> INT_PAIR_VECT;
  typedef INT_PAIR_VECT::const_iterator INT_PAIR_VECT_CI;
  
  typedef std::pair<double, INT_PAIR> DOUBLE_INT_PAIR;
  
  //! Sort a list of atoms by their  CIP rank 
  /*!
    \param mol        molecule of interest
    \param commAtms   atoms that need to be ranked
    \param ascending  sort to an ascending order or a descending order
  */
  template<class T> T rankAtomsByRank(const RDKit::ROMol &mol, const T &commAtms,
                                      bool ascending=true);
  
   
  //! computes a subangle for an atom of given hybridization and degree
  /*!
    \param degree the degree of the atom (number of neighbors)
    \param htype  the atom's hybridization

    \return the subangle (in radians)
  */
  inline double computeSubAngle(unsigned int degree, RDKit::Atom::HybridizationType htype) {
    double angle = M_PI;
    switch (htype) {
    case RDKit::Atom::UNSPECIFIED :
    case RDKit::Atom::SP3 :
      if (degree == 4) {
        angle = M_PI/2;
      } else {
        angle = 2*M_PI/3;
      }
      break;
    case RDKit::Atom::SP2 :
      angle = 2*M_PI/3;
      break;
    default:
      angle = 2.*M_PI/degree;
    }
    return angle;
  }
  
  //! computes the rotation direction between two vectors
  /*!

    Let:

       v1 = loc1 - center

       v2 = loc2 - center

    If remaining angle(v1, v2) is < 180 and corss(v1, v2) > 0.0 then the rotation dir is +1.0

    else if remAngle(v1, v2) is > 180 and cross(v1, v2) < 0.0 then rotation dir is -1.0

    else if remAngle(v1, v2) is < 180 and cross(v1, v2) < 0.0 then rotation dir is -1.0

    finally if remAngle(v1, v2) is > 180 and cross(v1, v2) < 0.0 then rotation dir is +1.0
    
    \param center     the common point
    \param loc1       endpoint 1   
    \param loc2       endpoint 2
    \param remAngle   the remaining angle about center in radians

    \return the rotation direction (1 or -1)
  */
  inline int rotationDir(const RDGeom::Point2D &center, const RDGeom::Point2D &loc1, 
                          const RDGeom::Point2D &loc2, double remAngle) {
     RDGeom::Point2D pt1 = loc1 - center;
     RDGeom::Point2D pt2 = loc2 - center;
     double cross = pt1.x*pt2.y - pt1.y*pt2.x;
     double diffAngle = M_PI - remAngle;
     cross  *= diffAngle;
     if (cross >= 0.0) {
       return -1;
     } else {
       return 1;
     }
   }

  //! computes and return the normal of a vector between two points
  /*!
    \param center     the common point
    \param other      the endpoint

    \return the normal
  */
  inline RDGeom::Point2D computeNormal(const RDGeom::Point2D &center, 
                                       const RDGeom::Point2D &other) {
    RDGeom::Point2D res = other - center;
    res.normalize();
    double tmp = res.x;
    res.x = -res.y;
    res.y = tmp;
    return res;
  }

  //! computes the rotation angle between two vectors
  /*!
    \param center     the common point
    \param loc1       endpoint 1   
    \param loc2       endpoint 2

    \return the angle (in radians)
  */
  inline double computeAngle(const RDGeom::Point2D &center, 
                             const RDGeom::Point2D &loc1, 
                             const RDGeom::Point2D &loc2) {
    RDGeom::Point2D v1 = loc1 - center;
    RDGeom::Point2D v2 = loc2 - center;
    return v1.angleTo(v2);
  }

  //! \brief pick the ring to embed first in a fused system
  /*!
    \param mol        the molecule of interest
    \param fusedRings the collection of the molecule's fused rings

    \return the index of the ring with the least number of substitutions
  */
  int pickFirstRingToEmbed(const RDKit::ROMol &mol, const RDKit::VECT_INT_VECT &fusedRings);

  //! \brief find the rotatable bonds on the shortest path between two atoms
  //!   we will ignore ring atoms, and double bonds which are marked cis/trans
  /*!
    <b>Note</b> that rotatable in this context doesn't connect to the
      standard chemical definition of a rotatable bond; we're just talking
      about bonds than can be flipped in order to clean up the depiction.
  
    \param mol   the molecule of interest
    \param aid1  index of the first atom
    \param aid2  index of the second atom

    \return a set of the indices of the rotatable bonds
  */
  RDKit::INT_VECT getRotatableBonds(const RDKit::ROMol &mol, unsigned int aid1,
                                    unsigned int aid2);

  //! \brief find all the rotatable bonds in a molecule
  //!   we will ignore ring atoms, and double bonds which are marked cis/trans
  /*!
    <b>Note</b> that rotatable in this context doesn't connect to the
      standard chemical definition of a rotatable bond; we're just talking
      about bonds than can be flipped in order to clean up the depiction.
  
    \param mol   the molecule of interest
    
    \return a set of the indices of the rotatable bonds
  */
  RDKit::INT_VECT getAllRotatableBonds(const RDKit::ROMol &mol);

  //! Get the ids of the atoms and bonds that are connected to aid
  void getNbrAtomAndBondIds(unsigned int aid, const RDKit::ROMol *mol,
                            RDKit::INT_VECT &aids, RDKit::INT_VECT &bids);

  //! Find pairs of bonds that can be permuted at a non-ring degree 4 atom
  /*!
    This function will return only those pairs that cannot be 
    permuted by flipping a rotatble bond
    
           D
           |
           b3
           |
      A-b1-B-b2-C
           |
           b4
           |
           E
     For example in teh above situation on the pairs (b1, b3) and (b1, b4) will be returned
     All other permutations can be achieved via a rotatable bond flip.

     ARGUMENTS:
     \param center - location of the central atom
     \param nbrBids - a vector (of length 4) containing the ids of the bonds to the neighbors
     \param nbrLocs - locations of the neighbors
  */
  INT_PAIR_VECT findBondsPairsToPermuteDeg4(const RDGeom::Point2D &center, const RDKit::INT_VECT &nbrBids, 
                                            const VECT_C_POINT &nbrLocs);
}

#endif
