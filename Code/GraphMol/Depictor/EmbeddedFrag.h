//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_EMBEDDED_FRAG_H_
#define _RD_EMBEDDED_FRAG_H_

#include <RDGeneral/types.h>
#include <Geometry/Transform2D.h>
#include <Geometry/point.h>
#include "DepictUtils.h" 
#include <boost/smart_ptr.hpp>


namespace RDKit {
  class ROMol;
  class Bond;
}

namespace RDDepict {
  typedef boost::shared_array<double> DOUBLE_SMART_PTR;
  
  //! Class that contains the data for an atoms that has alredy been embedded
  class EmbeddedAtom {
  public:
    typedef enum {
      UNSPECIFIED=0,
      CISTRANS,
      RING} EAtomType;

    EmbeddedAtom() : aid(0), angle(-1.0), nbr1(-1), nbr2(-1), 
                     CisTransNbr(-1), ccw(true), rotDir(0), d_density(-1.0){
      neighs.clear();
    }

    EmbeddedAtom(unsigned int aid, const RDGeom::Point2D &pos) :
      aid(aid), angle(-1.0), nbr1(-1), nbr2(-1), 
      CisTransNbr(-1), ccw(true), rotDir(0), d_density(-1.0) {
      loc = pos;
    }
      
    EmbeddedAtom& operator=(const EmbeddedAtom& other) {
      loc = other.loc;
      angle = other.angle;
      nbr1 = other.nbr1;
      nbr2 = other.nbr2;
      CisTransNbr = other.CisTransNbr;
      rotDir = other.rotDir;
      normal = other.normal;
      ccw = other.ccw;
      neighs = other.neighs;
      d_density = other.d_density;
      return *this;
    }

    void Transform(const RDGeom::Transform2D &trans) {
      RDGeom::Point2D temp = loc + normal;
      trans.TransformPoint(loc);
      trans.TransformPoint(temp);
      normal = temp - loc;
    }

    void Reflect(const RDGeom::Point2D &loc1,
                 const RDGeom::Point2D &loc2) {
      RDGeom::Point2D temp = loc + normal; 
      loc = reflectPoint(loc, loc1, loc2);
      temp = reflectPoint(temp, loc1, loc2);
      normal = temp - loc;
      ccw = (!ccw);
    }
    
    unsigned int aid; // the id of the atom

    //! the angle that is already takes at this atom, so any new atom attaching to this
    //! atom with have to fall in the available part
    double angle; 
    
    //! the first neighbor of this atom that form the 'angle'
    int nbr1;
    
    //! the second neighbor of atom that from the 'angle'
    int nbr2;
    
    //! is this is a cis/trans atom the neighbor of this atom that is involved in the
    //! cis/trans system - defaults to -1
    int CisTransNbr;

    //! which direction do we rotate this normal to add the next bond 
    //! if ccw is true we rotate counter cloack wise, otherwise rotate clock wise, by an angle that is
    //! <= PI/2
    bool ccw;

    //! rotation direction around this atom when adding new atoms,
    //! we determine this for the first neighbor and stick to this direction after that
    //! useful only on atoms that are degree >= 4
    int rotDir;

    RDGeom::Point2D loc; // the current location of this atom
    //! this is a normal vector to one of the bonds that added this atom
    //! it provides the side on which we want to add a new bond to this atom,
    //! this is only relevant when we are dealing with non ring atoms. We would like to draw chains in
    //! a zig-zag manner
    RDGeom::Point2D normal;

    //! and these are the atom IDs of the neighbors that still need to be embedded
    RDKit::INT_VECT neighs;

    // density of the atoms around this atoms
    // - this is sum of inverse of the square of distances to other atoms from this atom
    // Used in the collision removal code - initialized to -1.0
    double d_density;
  };



  typedef std::map<unsigned int, EmbeddedAtom> INT_EATOM_MAP;
  typedef INT_EATOM_MAP::iterator INT_EATOM_MAP_I;
  typedef INT_EATOM_MAP::const_iterator INT_EATOM_MAP_CI;
  

  //! Class containing a fragment of a molecule that has already been embedded
  /*
    Here is how this class is designed to be used
    - find a set of fused rings and compte the coordinates for the atoms in those ring
    - them grow this sytem either by adding non ring neighbors
    - or by adding other embedded fragment
    - so at the end of the process  the whole molecule end up being one these embedded frag objects
  */
  class EmbeddedFrag {
    // REVIEW: think about moving member functions up to global level and just using
    // this class as a container

  public:
    //! Default constructor
    EmbeddedFrag() : d_done(false), dp_mol(0){
      d_eatoms.clear();
      d_attachPts.clear();
    };

    //! Intializer from a single atom id 
    /*!
      A single Embedded Atom with this atom ID is added and placed at the origin
    */
    EmbeddedFrag(unsigned int aid, const RDKit::ROMol *mol);

    //! Constructor when the coordinates have been specified for a set of atoms
    /*! 
       This simply initialized a set of EmbeddedAtom to have the same coordinates as the 
       one's specified. No testing is done to verify any kind of ocrrectlness. Also
       this fragment is less ready (to expand and add new neighbors) than when using other 
       constructors. This is because:
       - the user may have specified coords for only a part of the atoms in a fused ring systems
         in which case we need to find these atoms and merge these ring systems to this fragment
       - The atoms are not yet aware of their neighbor (what is left to add etc.) this again
         depends on  atoms properly so that new 
         neighbors can be added to them
    */
    EmbeddedFrag(const RDKit::ROMol *mol, const RDGeom::INT_POINT2D_MAP &coordMap);

    //! Initializer from a set of fused rings
    /*!
      ARGUMENTS:
      \param mol        the molecule of interest
      \param fusedRings a vector of rings, each ring is a list of atom ids 
    */
    EmbeddedFrag(const RDKit::ROMol *mol, const RDKit::VECT_INT_VECT &fusedRings);

    //! Initializer for a cis/trans system using the double bond
    /*!
      ARGUMENTS:
      \param dblBond   the double bond that is involed in the cis/trans configuration
    */
    explicit EmbeddedFrag(const RDKit::Bond* dblBond);

    //! Expand this embedded system by adding negihboring atoms or other embedded systems
    /*!

      Note that both nratms and efrags are modified in this function
      as we start merging them with the current fragment

    */
    void expandEfrag(RDKit::INT_LIST &nratms, std::list<EmbeddedFrag> &efrags);
                     
    //! Add a new non-ring atom to this object
    /*
      ARGUMENTS:
      \param aid     ID of the atom to be added
      \param toAid   ID of the atom that is already in this object to which this atom is added
    */
    void addNonRingAtom(unsigned int aid, unsigned int toAid);

    //! Merge this embedded object with another embedded fragment
    /*!

      The transformation (rotation + translation required to attached
      the passed in object will be computed and applied. The
      coordinates of the atoms in this object will remain fixed We
      will assume that there are no common atoms between the two
      fragments to start with
      
      ARGUMENTS:
      \param  embObj  another EmbeddedFrag object to be merged with this object
      \param  toAid   the atom in this embedded fragment to which the new object will be attached
      \param  nbrAid  the atom in the other fragment to attach to
    */
    void mergeNoCommon(EmbeddedFrag &embObj, unsigned int toAid, 
                       unsigned int nbrAid);

    //! Merge this embedded object with another embedded fragment
    /*!

      The transformation (rotation + translation required to attached
      the passed in object will be computed and applied. The
      coordinates of the atoms in this object will remain fixed This
      already know there are a atoms in common and we will use them to
      merge things
      
      ARGUMENTS:
      \param embObj    another EmbeddedFrag object to be merged with this object
      \param commAtms  a vector of ids of the common atoms

    */
    void mergeWithCommon(EmbeddedFrag &embObj, RDKit::INT_VECT &commAtms);

    void mergeFragsWithComm(std::list<EmbeddedFrag> &efrags);

    //! Mark this fragment to be done for final embedding
    void markDone() {
      d_done = true;
    }

    //! If this fragment done for the final embedding
    bool isDone() {
      return d_done;
    }

    //! Get the molecule that this embedded fragmetn blongs to
    const RDKit::ROMol *getMol() const { return dp_mol;}

    //! Find the common atom ids between this fragment and a second one
    RDKit::INT_VECT findCommonAtoms(const EmbeddedFrag &efrag2);

    //! Find a neighbor to a non-ring atom among the already embedded atoms
    /*!
      ARGUMENTS:
      \param aid  the atom id of interest
     
      RETURNS:
      \return the id of the atom if we found a neighbor
                   -1 otherwise
                   
      NOTE: by definition we can have only one neighbor in the embdded system. 
    */
    int findNeighbor(unsigned int aid);

    //! Tranform this object to a new coordinates system
    /*!
      ARGUMENTS:
      \param trans : the transformation that need to be applied to the atoms in this object
    */
    void Transform(const RDGeom::Transform2D &trans);

    void Reflect(const RDGeom::Point2D &loc1,
                 const RDGeom::Point2D &loc2);

    const INT_EATOM_MAP &GetEmbeddedAtoms() const {
      return d_eatoms;
    }

    void Translate(const RDGeom::Point2D &shift) {
      INT_EATOM_MAP_I eari;
      for (eari = d_eatoms.begin(); eari != d_eatoms.end(); eari++) {
        eari->second.loc += shift;
      }
    }

    EmbeddedAtom GetEmbeddedAtom(unsigned int aid) const {
      INT_EATOM_MAP_CI posi = d_eatoms.find(aid);
      if (posi == d_eatoms.end()) {
        PRECONDITION(0, "Embedded atom does not contain embedded atom specified");
      }
      return posi->second;
    }

    //! the number of atoms in the embedded system
    int Size() const {
      return d_eatoms.size();
    }
    
    //! \brief compute a box that encloses the fragment
    void computeBox();

    //! \brief Flip atoms on one side of a bond - used in removing colissions
    /*!
      ARGUMENTS:
      \param bondId - the bond used as the mirror to flip
      \param flipEnd - flip the atoms at the end of the bond

    */
    void flipAboutBond(unsigned int bondId,bool flipEnd=true);

    void openAngles(const double *dmat, unsigned int aid1, unsigned int aid2);

    std::vector<PAIR_I_I> findCollisions(const double *dmat, bool includeBonds=1);

    void computeDistMat(DOUBLE_SMART_PTR &dmat);

    double mimicDistMatAndDensityCostFunc(const DOUBLE_SMART_PTR *dmat, 
                                          double mimicDmatWt);

    void permuteBonds(unsigned int aid, unsigned int aid1, unsigned int aid2);

    void randomSampleFlipsAndPermutations(unsigned int nBondsPerSample=3,
                                          unsigned int nSamples=100, int seed=100,
                                          const DOUBLE_SMART_PTR *dmat=0, 
                                          double mimicDmatWt=0.0,
                                          bool permuteDeg4Nodes=false);

    //! Remove collisions in a structure by flipping rotable bonds 
    //! along the shortest path between two colliding atoms
    void removeCollisionsBondFlip();

    //! Remove collision by opening angles at the offending atoms
    void removeCollisionsOpenAngles();

    //! Remove collisions by shortening bonds along the shortest path between the atoms
    void removeCollisionsShortenBonds();

    //! helpers funtions to 
    

    //! \brief make list of neighbors for each atom in the embedded system that
    //!  still need to be embedded
    void setupNewNeighs(); 

    //! update the  unembedded neighbor atom list for a specified atom
    void updateNewNeighs(unsigned int aid); 

    //! \brief Find all atoms in this embedded system that are
    //!  within a specified distant of a point
    int findNumNeigh(const RDGeom::Point2D &pt, double radius);


    inline double getBoxPx() {return d_px;}
    inline double getBoxNx() {return d_nx;}
    inline double getBoxPy() {return d_py;}
    inline double getBoxNy() {return d_ny;}

    void canonicalizeOrientation();

    
  private:

    double totalDensity();

    void embedFusedRings(const RDKit::VECT_INT_VECT &fusedRings);

    //! \brief Find a transform to join a ring to the current embedded frag when we
    //! have only on common atom
    /*!
      So this is the state of affairs assumed here:
      - we already have some rings in the fused system embeded and the
        coordinates for the atoms
      - the coordinates for the atoms in the new ring (with the center
        of rings at the origin) are available nringCors. we want to
        translate and rotate this ring to join with the already
        embeded rings.
      - only one atom is common between this new ring and the atoms
        that are already embedded
      - so we need to compute a transform that includes a translation
        so that the common atom overlaps and the rotation to minimize
        overalp with other atoms.
        
      Here's what is done:
      - we bisect the remaining sweep angle at the common atom and
        attach the new ring such that the center of the new ring falls
        on this bisecting line
        
      NOTE: It is assumed here that the original coordinates for the
      new ring are such that the center is at the origin (this is the
      way rings come out of embedRing)
    */
    RDGeom::Transform2D computeOneAtomTrans(unsigned int commAid,
                                            const EmbeddedFrag &other);
    

    RDGeom::Transform2D computeTwoAtomTrans(unsigned int aid1, unsigned int aid2, 
                                            const RDGeom::INT_POINT2D_MAP &nringCor);
        
    //! Merge a ring with already embedded atoms
    /*!
      It is assumed that the new rings has already been oriented
      correctly etc.  This function just update all the relevant data,
      like the neighbor information and the sweep angle
    */
    void mergeRing(const EmbeddedFrag &embRing, unsigned int nCommon, 
                   const RDKit::INT_VECT &pinAtoms);

    //! Reflect a fragment if necessary through a line connecting two atoms
    /*!

      We want add the new fragment such that, most of its atoms fall
      on the side opoiste to where the atoms already embedded are aid1
      and aid2 give the atoms that were used to align the new ring to
      the embedded atoms and we will assume that that process has
      already taken place (i.e. transformRing has been called)

    */
    void reflectIfNecessaryDensity(EmbeddedFrag &embFrag, unsigned int aid1, 
                                   unsigned int aid2);

    //! Reflect a fragment if necessary based on the cis/trans specification
    /*!

      we want to add the new fragment such that the cis/trans
      specification on bond between aid1 and aid2 is not violated. We
      will assume that aid1 and aid2 from this fragments as well as
      embFrag are already aligned to each other.
      
      \param embFrag   the fragment that will be reflected if necessary
      \param ctCase    which fragment if the cis/trans dbl bond
                        - 1 means embFrag is the cis/trans fragment
                        - 2 mean "this" is the cis/trans fragment
      \param aid1      first atom that forms the plane (line) of reflection 
      \param aid2      seconf atom that forms the plane of reflection
    */
    void reflectIfNecessaryCisTrans(EmbeddedFrag &embFrag, unsigned int ctCase, 
                                    unsigned int aid1, unsigned int aid2);

    //! Reflect a fragment if necessary based on a thrid common point
    /*!

      we want add the new fragment such that the thrid point falls on
      the same side of aid1 and aid2. We will assume that aid1 and
      aid2 from this fragments as well as embFrag are already aligned
      to each other.

    */
    void reflectIfNecessaryThirdPt(EmbeddedFrag &embFrag, unsigned int aid1, unsigned int aid2, 
                                   unsigned int aid3);

    //! \brief Initialize this fragment from a ring and coordinates for its atoms
    /*!
      ARGUMENTS:
      /param ring     a vector of atom ids in the ring; it is assumed that there in 
                      clockwise or anti-clockwise order
      /param nringMap a map of atomId to coordinate map for the atoms in the ring
    */
    void initFromRingCoords(const RDKit::INT_VECT &ring, 
                            const RDGeom::INT_POINT2D_MAP &nringMap);

    //! Helper function to addNonRingAtom to a specified atoms in the fragment
    /*
      Add an atom to this embedded fragment when the fragment already
      has a atleast two previously added neighbors to 'toAid'. In this
      case we have to choose where the the new neighbor goes based on
      the angle that is already taken around the atom.
      
      ARGUMENTS:
      \param aid    ID of the atom to be added
      \param toAid  ID of the atom that is already in this object to which this atom is added
    */
    void addAtomToAtomWithAng(unsigned int aid, unsigned int toAid);


    //! Helper function to addNonRingAtom to a specified atoms in the fragment
    /*!

      Add an atom (aid) to an atom (toAid) in this embedded fragment
      when 'toAid' has one or no neighbors previously added. In this
      case where the new atom should fall is determined by the degree
      of 'toAid' and the congestion around it.
      
      ARGUMENTS:
      \param  aid     ID of the atom to be added
      \param  toAid   ID of the atom that is already in this object to which this atom is added
      \param  mol     the molecule we are dealing with
    */
    void addAtomToAtomWithNoAng(unsigned int aid, 
                                unsigned int toAid); //, const RDKit::ROMol *mol);

    //! Helper funtion to contructor that takes predefined coordinates
    /*! 

      Given an atom with more than 2 neighbors all embedded in this
      fragment this function tries to determine

      - how much of an angle if left for any new neighbors yet to be
        added
      - which atom should we rotate when we add a new neighbor and in
        which direction (clockwise or anticlockwise
        
      This is how it works 
      - find the pair of nbrs that have the largest angle
      - this will most likely be the angle that is available - unless
        we have fused rings and we found on of the ring angle !!!! -
        in this cae we find the next best
      - find the smallest anngle that contains one of these nbrs -
        this determined which
      - way we want to rotate
     
      ARGUMENTS:
      \param aid        the atom id where we are centered right now
      \param doneNbrs   list of neighbors that are already embedded around aid
    */
    void computeNbrsAndAng(unsigned int aid, const RDKit::INT_VECT &doneNbrs);
    //const RDKit::ROMol *mol);

    //! are we embedded with the final (molecule) coordinates
    bool d_done;
    double d_px, d_nx, d_py, d_ny;

    //! a map that takes one from teh atom id to the embeddedatom object for that atom.
    INT_EATOM_MAP d_eatoms;

    //RDKit::INT_DEQUE d_attachPts;
    RDKit::INT_LIST d_attachPts;

    // pointer to the owning molecule
    const RDKit::ROMol *dp_mol;

  };

  
}

#endif
