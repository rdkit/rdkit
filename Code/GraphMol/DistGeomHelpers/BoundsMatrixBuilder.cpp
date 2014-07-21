// $Id$
//
//  Copyright (C) 2004-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <DistGeom/BoundsMatrix.h>
#include "BoundsMatrixBuilder.h"
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <ForceField/UFF/BondStretch.h>
#include <Geometry/Utils.h>

#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>
#include <Numerics/SymmMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

const double DIST12_DELTA = 0.01;
const double ANGLE_DELTA = 0.0837;
const double RANGLE_DELTA = 0.0837; // tolerance for bond angles
const double TANGLE_DELTA = 0.0837; // tolerance for torsion angle
const double DIST13_TOL = 0.04;
const double GEN_DIST_TOL = 0.06; //  a general distance tolerance
const double DIST15_TOL = 0.08;
const double VDW_SCALE_15 = 0.7;
const double MAX_UPPER = 1000.0;
#include <map>

namespace RDKit {
  namespace DGeomHelpers {
    // forward declarations:
    typedef boost::shared_ptr<RDNumeric::IntSymmMatrix> SymmIntMatPtr;
    typedef boost::shared_ptr<RDNumeric::DoubleSymmMatrix> SymmDoubleMatPtr;
    

    typedef boost::dynamic_bitset<> BIT_SET;

    //! Bunch of functions to set distance bound based on topology

    typedef std::map<int, double> INT_DOUBLE_MAP;
    typedef INT_DOUBLE_MAP::const_iterator INT_DOUBLE_MAP_CI;
    
    typedef std::vector<long int> LINT_VECT;


    //! A structure used to store planar 14 paths - cis/trans
    struct Path14Configuration {
      unsigned int bid1, bid2, bid3;
      typedef enum {
        CIS = 0,
        TRANS,
        OTHER,
      } Path14Type;
      Path14Type type; 
    };
    
    typedef std::vector<Path14Configuration> PATH14_VECT;
    typedef PATH14_VECT::iterator PATH14_VECT_I;
    typedef PATH14_VECT::const_iterator PATH14_VECT_CI;

    class ComputedData {
    public: 
      ComputedData(unsigned int nAtoms, unsigned int nBonds) {
        bondLengths.resize(nBonds);
        RDNumeric::IntSymmMatrix *bAdj = new RDNumeric::IntSymmMatrix(nBonds, -1);
        bondAdj.reset(bAdj);
        RDNumeric::DoubleSymmMatrix *bAngles = new RDNumeric::DoubleSymmMatrix(nBonds, -1.0);
        bondAngles.reset(bAngles);
        cisPaths.resize(nBonds*nBonds*nBonds);
        transPaths.resize(nBonds*nBonds*nBonds);
        set15Atoms.resize(nAtoms*nAtoms);
      }
      
      ~ComputedData(){}

      DOUBLE_VECT bondLengths;
      SymmIntMatPtr bondAdj; // bond adjacency matrix
      SymmDoubleMatPtr bondAngles;
      PATH14_VECT paths14;
      BIT_SET cisPaths;
      BIT_SET transPaths;
      BIT_SET set15Atoms;
    };



    //! Set 1-2 distance bounds for between atoms in a molecule
    /*!   
      These are mostly bond lengths obtained from UFF parameters and then
      adjusted by a small tolerance to set the upper and lower limits
      
      \param mol          The molecule of interest
      \param mmat         Bounds matrix to which the bounds are written
      \param accumData    Used to store the data that have been calculated so far
                          about the molecule

    */
    void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, 
                     ComputedData &accumData);
    
    //! Set 1-3 distance bounds for atoms in a molecule
    /*!  
      These are computed using bond angles and bond lengths. There are special
      cases here, in particular for 3, 4, and 5 membered rings. Special attention
      is also paid to fused ring ring systems, when setting bounds on atoms that 
      have an atom between that is shared by multiple rings.
      
      \param mol          Molecule of interest
      \param mmat         Bounds matrix to which the bounds are written
      \param accumData    Used to store the data that have been calculated so far
                          about the molecule

      <b>Procedure</b>
      All 1-3 distances within all simple rings are first dealt with, while keeping track of 
      any atoms that are visited twice; these are the atoms that are part of multiple simple rings.
      Then all other 1-3 distance are set while treating 1-3 atoms that have a ring atom in 
      between differently from those that have a non-ring atom in between.
     */
    void set13Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, 
                     ComputedData &accumData);
    
    //! Set 1-4 distance bounds for atoms in a molecule
    /*!
      These are computed using the range of allowed torsion angles. There are several
      special casses and the rest are computed using 0 and 180 deg as the min. 
      and max. torsion angles. The special cases deal with ring systems, double bonds
      with cis-trans specifications, and a few special sub-structures
      
      \param mol          Molecule of interest
      \param mmat         Bounds matrix to which the bounds are written
      \param accumData    Used to store the data that have been calculated so far
                          about the molecule

      <b>Procedure</b>
      As in the case of 1-3 distances 1-4 distance that are part of simple rings are
      first dealt with. The remaining 1-4 cases are dealt with while paying attention
      to the special cases.
     */
    void set14Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, 
                     ComputedData &accumData);
    
    //! Set 1-5 distance bounds for atoms in a molecule
    /*!
      This is an optional call that recognizes a variety of special cases.
      
      \param mol          Molecule of interest
      \param mmat         Bounds matrix to which the bounds are written
      \param accumData    Used to store the data that have been calculated so far
                          about the molecule
    */
    void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                     ComputedData &accumData,double *distMatrix);

    //! Set lower distance bounds based on VDW radii for atoms that are not covered by 
    //! other bounds (1-2, 1-3, 1-4, or 1-5)
    /*!
      
      \param mol             Molecule of interest
      \param mmat            Bounds matrix to which the bounds are written
      \param useTopolScaling If true scale the sum of the vdW radii while setting lower bounds
                             so that a smaller value (0.7*(vdw1 + vdw2) ) is used for paths
	     		     that are less 5 bonds apart.
    */
    void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool useTopolScaling=true);

  }
}

namespace RDKit {
  namespace DGeomHelpers {
    void _checkAndSetBounds(unsigned int i, unsigned int j, 
                            double lb, double ub, 
                            DistGeom::BoundsMatPtr mmat) {
      // get the exisiting bounds
      double clb = mmat->getLowerBound(i, j);
      double cub = mmat->getUpperBound(i, j);

      CHECK_INVARIANT(ub>lb, "upper bound not greater than lower bound");
      CHECK_INVARIANT(lb>DIST12_DELTA || clb>DIST12_DELTA, "bad lower bound");

      if (clb <= DIST12_DELTA) {
        mmat->setLowerBound(i, j, lb);
      } else {
        if ((lb < clb) && (lb > DIST12_DELTA) ){
          mmat->setLowerBound(i, j, lb); // conservative bound setting
        }
      }
      
      if (cub >= MAX_UPPER) { //FIX this
        mmat->setUpperBound(i, j, ub);
      } else {
        if ((ub > cub) && (ub < MAX_UPPER)) {
          mmat->setUpperBound(i, j, ub);
        }
      }
    }
    
    void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, ComputedData &accumData) {
      unsigned int npt = mmat->numRows();
      CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
      CHECK_INVARIANT(accumData.bondLengths.size() >= mol.getNumBonds(), "Wrong size accumData");
      bool foundAll;
      UFF::AtomicParamVect atomParams;
      boost::tie(atomParams,foundAll)= UFF::getAtomTypes(mol);
      CHECK_INVARIANT(atomParams.size()==mol.getNumAtoms(),"parameter vector size mismatch");
      
      ROMol::ConstBondIterator bi;
      unsigned int begId, endId;
      double bl;
      for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
        begId = (*bi)->getBeginAtomIdx();
        endId = (*bi)->getEndAtomIdx();
        if(atomParams[begId] && atomParams[endId]){
          bl = ForceFields::UFF::Utils::calcBondRestLength((*bi)->getBondTypeAsDouble(),
                                                           atomParams[begId], atomParams[endId]);
          accumData.bondLengths[(*bi)->getIdx()] = bl;
          mmat->setUpperBound(begId, endId, bl + DIST12_DELTA);
          mmat->setLowerBound(begId, endId, bl - DIST12_DELTA);
        } else {
          // we don't have parameters for one of the atoms... so we're forced to use 
          // very crude bounds:
          double vw1 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(begId)->getAtomicNum());
          double vw2 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(endId)->getAtomicNum());
          double bl=(vw1+vw2)/2;
          accumData.bondLengths[(*bi)->getIdx()] = bl;
          mmat->setUpperBound(begId, endId, 1.5*bl);
          mmat->setLowerBound(begId, endId, .5*bl );
        }
      }
    }
    
    void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool useTopolScaling,
                          double *dmat) {
      unsigned int npt = mmat->numRows();
      PRECONDITION(npt == mol.getNumAtoms(), "Wrong size metric matrix");
      unsigned int i, j;

      double vw1, vw2;

      for (i = 1; i < npt; i++) {
        vw1 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(i)->getAtomicNum());
        for (j = 0; j < i; j++) {
          vw2 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(j)->getAtomicNum());
          if (mmat->getLowerBound(i,j) < DIST12_DELTA) {
            // ok this is what we are going to do
            // - for atoms that are 4 or 5 bonds apart (15 or 16 distances), we will scale
            //   the sum of the VDW radii so that the atoms can get closer
            //   For 15 we will use VDW_SCALE_15 and for 16 we will use 1 - 0.5*VDW_SCALE_15
            // - for all other pairs of atoms more than 5 bonds apart we use the sum of the VDW radii
            //    as the lower bound
            if (dmat[i*npt + j] == 4.0) {
              mmat->setLowerBound(i,j, VDW_SCALE_15*(vw1+vw2));
            } else if (dmat[i*npt + j] == 5.0) {
              mmat->setLowerBound(i,j, (VDW_SCALE_15 + 0.5*(1.0-VDW_SCALE_15))*(vw1+vw2));
            } else {
              mmat->setLowerBound(i,j, (vw1+vw2));
            }
          }
        }
      }
    }

    void _set13BoundsHelper(unsigned int aid1, unsigned int aid, unsigned int aid3, 
                            double angle, const ComputedData &accumData, DistGeom::BoundsMatPtr mmat, 
                            const ROMol &mol) {
      unsigned int bid1 = mol.getBondBetweenAtoms(aid1, aid)->getIdx();
      unsigned int bid2 = mol.getBondBetweenAtoms(aid, aid3)->getIdx();
      double dl = RDGeom::compute13Dist(accumData.bondLengths[bid1], accumData.bondLengths[bid2], 
                                        angle);
      double du = dl + DIST13_TOL; 
      dl -= DIST13_TOL;
      _checkAndSetBounds(aid1, aid3, dl, du, mmat);
    }
    
    void _setRingAngle(Atom::HybridizationType aHyb, unsigned int ringSize, 
                       double &angle) {
      // NOTE: this assumes that all angles in a ring are equal. This is
      // certainly not always the case, particular in aromatic rings with heteroatoms
      // like s1cncc1. This led to GitHub55, which was fixed elsewhere.
      
      if ((aHyb == Atom::SP2) || (ringSize==3) || (ringSize==4)) {
        angle = M_PI*(1 - 2.0/ringSize);
      } else if (aHyb == Atom::SP3) {
        if (ringSize == 5) {
          angle = 104*M_PI/180;
        } else {
          angle = 109.5*M_PI/180;
        }
      } else if (aHyb == Atom::SP3D) {
        angle = 105.0*M_PI/180;
      } else if (aHyb == Atom::SP3D2) {
        angle = 90.0*M_PI/180;
      } else {
        angle = 120*M_PI/180;
      }
    }
    

    struct lessVector : public std::binary_function<INT_VECT,INT_VECT,bool> {
      bool operator()(const INT_VECT &v1, const INT_VECT &v2) const {
        return v1.size() < v2.size();
      }
    };

  void set13Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, ComputedData &accumData) {
      unsigned int npt = mmat->numRows();
      CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
      CHECK_INVARIANT(accumData.bondAngles->numRows() == mol.getNumBonds(), "Wrong size bond angle matrix");
      CHECK_INVARIANT(accumData.bondAdj->numRows() == mol.getNumBonds(), "Wrong size bond adjacency matrix");

      // Since most of the special cases arise out of ring system, we will do the following here:
      // - Loop over all the rings and set the 13 distances between atoms in these rings. 
      //   While doing this keep track of the ring atoms that have already been used as the center atom.
      // - Set the 13 distance between atoms that have a ring atom in between; these can be either non-ring atoms, 
      //   or a ring atom and a non-ring atom, or ring atoms that belong to different simple rings
      // - finally set all other 13 distances
      const RingInfo *rinfo = mol.getRingInfo();
      CHECK_INVARIANT(rinfo, "");
      ROMol::OEDGE_ITER beg1, beg2, end1, end2;
      
      unsigned int aid2, aid1, aid3, bid1, bid2;
      double angle;

      VECT_INT_VECT atomRings = rinfo->atomRings();
      std::sort(atomRings.begin(), atomRings.end(), lessVector()); 
      // sort the rings based on the ring size
      VECT_INT_VECT_CI rii;
      INT_VECT visited(npt, 0);

      DOUBLE_VECT angleTaken(npt, 0.0);
      unsigned int i;
      unsigned int nb = mol.getNumBonds();
      BIT_SET donePaths(nb*nb);
      // first deal with all rings and atoms in them
      unsigned int id1, id2;
      for (rii = atomRings.begin(); rii != atomRings.end(); rii++) {
        unsigned int rSize = rii->size();
        aid1 = (*rii)[rSize-1];
        for (i = 0; i < rSize; i++) {
          aid2 = (*rii)[i];
          if (i == rSize-1) {
            aid3 = (*rii)[0];
          } else {
            aid3 = (*rii)[i+1];
          }
          const Bond *b1=mol.getBondBetweenAtoms(aid1, aid2);
          const Bond *b2=mol.getBondBetweenAtoms(aid2, aid3);
          CHECK_INVARIANT(b1,"no bond found");
          CHECK_INVARIANT(b2,"no bond found");
          bid1 = b1->getIdx();
          bid2 = b2->getIdx();
          id1 = nb*bid1 + bid2;
          id2 = nb*bid2 + bid1;

          if ((!donePaths[id1]) && (!donePaths[id2])) {
            // this invar stuff is to deal with bridged systems (Issue 215). In bridged 
            // systems we may be covering the same 13 (ring) paths multiple times and unnecssarily 
            // increasing the angleTaken at the central atom.
            _setRingAngle(mol.getAtomWithIdx(aid2)->getHybridization(), rSize, angle);
            _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
            accumData.bondAngles->setVal(bid1, bid2, angle);
            accumData.bondAdj->setVal(bid1, bid2, aid2);
            visited[aid2] += 1;
            angleTaken[aid2] += angle;
            donePaths[id1] = 1;
            donePaths[id2] = 1;
            //donePaths.push_back(invar);
          }
          aid1 = aid2;
        }
      }
      
      // now deal with the remaining atoms 
      for (aid2 = 0; aid2 < npt; aid2++) {
        const Atom *atom =  mol.getAtomWithIdx(aid2);
        unsigned int deg = atom->getDegree();
        unsigned int n13 = deg*(deg-1)/2;
        if (n13 == static_cast<unsigned int>(visited[aid2])) {
          // we are done with this atoms
          continue;
        }
        Atom::HybridizationType ahyb = atom->getHybridization();
        boost::tie(beg1,end1) = mol.getAtomBonds(atom);
        if (visited[aid2] >= 1) {
          // deal with atoms that we already visited; i.e. ring atoms. Set 13 distances for one of following cases:
          //  1) Non-ring atoms that have a ring atom in-between
          //  2) Non-ring atom and a ring atom that have a ring atom in between
          //  3) Ring atoms that belong to different rings (that are part of a fused system
          
          while (beg1 != end1) {
            const BOND_SPTR bnd1 = mol[*beg1];
            bid1 = bnd1->getIdx();
            aid1 = bnd1->getOtherAtomIdx(aid2);
            boost::tie(beg2,end2) = mol.getAtomBonds(atom);
            while (beg2 != beg1) {
              const BOND_SPTR bnd2 = mol[*beg2];
              bid2 = bnd2->getIdx();
              //invar = firstThousandPrimes[bid1]*firstThousandPrimes[bid2];
              if (accumData.bondAngles->getVal(bid1, bid2) < 0.0) {
              //if (bondAngles.find(invar) == bondAngles.end()) {
                // if we haven't dealt with these two bonds before
                
                // if we have a sp2 atom things are planar - we simply divide the remaining angle among the 
                // remaining 13 configurations (and there should only be one)
                if (ahyb == Atom::SP2) {
                  angle = (2*M_PI - angleTaken[aid2])/(n13-visited[aid2]);
                } else if (ahyb == Atom::SP3) {
                  // in the case of sp3 we will use the tetrahedral angle mostly - but 
                  // but with some special cases
                  angle = 109.5*M_PI/180;
                  // we will special-case a little bit here for 3, 4 members ring atoms that are sp3 hybirdized
                  // beyond that the angle reasonably close to the tetrahedral angle
                  if (rinfo->isAtomInRingOfSize(aid2, 3)) {
                    angle = 116.0*M_PI/180;
                  } else if (rinfo->isAtomInRingOfSize(aid2, 4)) {
                    angle = 112.0*M_PI/180;
                  }
                } else {
                  // other options we will simply based things on the number of substituent
                  if (deg == 5) {
                    angle = 105.0*M_PI/180;
                  } else if (deg == 6) {
                    angle = 135.0*M_PI/180;
                  } else { 
                    angle = 120.0*M_PI/180; // FIX: this default is probably not the best we can do here
                  }
                }
                aid3 = bnd2->getOtherAtomIdx(aid2);
                _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
                accumData.bondAngles->setVal(bid1, bid2, angle);
                accumData.bondAdj->setVal(bid1, bid2, aid2);
                angleTaken[aid2] += angle;
                visited[aid2] += 1;
              } 
              ++beg2;
            } // while loop over the second bond
            ++beg1;
          } // while loop over the first bond
        } else if (visited[aid2] == 0) { 
          // non-ring atoms - we will simply use angles based on hydridization
          while (beg1 != end1) {
            const BOND_SPTR bnd1 = mol[*beg1];
            bid1 = bnd1->getIdx();
            aid1 = bnd1->getOtherAtomIdx(aid2);
            boost::tie(beg2,end2) = mol.getAtomBonds(atom);
            while (beg2 != beg1) {
              const BOND_SPTR bnd2 = mol[*beg2];
              bid2 = bnd2->getIdx();
              if (ahyb == Atom::SP) {
                angle = M_PI;
              } else if (ahyb == Atom::SP2) {
                angle = 2*M_PI/3;
              } else if (ahyb == Atom::SP3) {
                angle = 109.5*M_PI/180;
              } else if (ahyb == Atom::SP3D){ 
                //FIX: this and the remaining two hybirdization states below should probably be special
                // cased. These defaults below are probably not the best we can do particularly when
                // stereo chemistry is know
                angle = 105.0*M_PI/180; 
              } else if (ahyb == Atom::SP3D2){
                angle = 135.0*M_PI/180;
              } else {
                angle = 120.0*M_PI/180;
              }
              aid3 = bnd2->getOtherAtomIdx(aid2);
              _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
              accumData.bondAngles->setVal(bid1, bid2, angle);
              accumData.bondAdj->setVal(bid1, bid2, aid2);
              angleTaken[aid2] += angle;
              visited[aid2] += 1;     
              ++beg2;
            } // while loop over second bond
            ++beg1;
          } // while loop over first bond
        } // done with non-ring atoms
      } // done with all atoms
    } // done with 13 distance setting

    Bond::BondStereo _getAtomStereo(const Bond *bnd, unsigned int aid1, 
                                          unsigned int aid4) {
      Bond::BondStereo stype = bnd->getStereo();
      if (stype > Bond::STEREOANY ) {
        const INT_VECT &stAtoms = bnd->getStereoAtoms();
        if ((static_cast<unsigned int>(stAtoms[0]) != aid1) ^ 
            (static_cast<unsigned int>(stAtoms[1]) != aid4)) {
          if (stype == Bond::STEREOZ) {
            stype = Bond::STEREOE;
          } else if (stype == Bond::STEREOE) {
            stype = Bond::STEREOZ;
          }
        }
      }
      return stype;

    }

    void _setInRing14Bounds(const ROMol &mol, const Bond* bnd1, const Bond* bnd2, const Bond* bnd3,
                            ComputedData &accumData, DistGeom::BoundsMatPtr mmat,double *dmat) {
      PRECONDITION(bnd1, "");
      PRECONDITION(bnd2, "");
      PRECONDITION(bnd3, "");
      unsigned int bid1, bid2, bid3;
      bid1 = bnd1->getIdx();
      bid2 = bnd2->getIdx();
      bid3 = bnd3->getIdx();
      const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2)); 
      PRECONDITION(atm2, "");
      Atom::HybridizationType ahyb2 = atm2->getHybridization(); 
      const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3)); 
      PRECONDITION(atm3, "");
      Atom::HybridizationType ahyb3 = atm3->getHybridization();
      
      unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
      unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());

      // check that this actually is a 1-4 contact:
      if(dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<2.9){
        //std::cerr<<"skip: "<<aid1<<"-"<<aid4<<" because d="<<dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<<std::endl;
        return;
      }
      
      double bl1 = accumData.bondLengths[bid1];
      double bl2 = accumData.bondLengths[bid2];
      double bl3 = accumData.bondLengths[bid3];

      double ba12 = accumData.bondAngles->getVal(bid1, bid2);
      double ba23 = accumData.bondAngles->getVal(bid2, bid3);

      CHECK_INVARIANT(ba12 > 0.0, "");
      CHECK_INVARIANT(ba23 > 0.0, "");
      double dl, du;
      unsigned int nb = mol.getNumBonds();
      //several special cases here
      Path14Configuration path14;
      path14.bid1 = bid1; path14.bid2 = bid2; path14.bid3 = bid3;
      Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
      if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2) 
          && (stype != Bond::STEREOE)) { 
        dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        du = dl + 2*GEN_DIST_TOL;
        path14.type = Path14Configuration::CIS;
        accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
        accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
      } else {
        // basically we will assume 0 to 180 allowed
        dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
        du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
        if(du<dl) std::swap(du,dl);
        if (fabs(du-dl) < DIST12_DELTA) {
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;
        }
        path14.type = Path14Configuration::OTHER;
      }
      
      //std::cerr<<"7: "<<aid1<<"-"<<aid4<<std::endl;
      _checkAndSetBounds(aid1, aid4, dl , du, mmat);
      accumData.paths14.push_back(path14);
      
    }
    
    void _setTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2, 
                                   const Bond *bnd3, ComputedData &accumData,
                                   DistGeom::BoundsMatPtr mmat,double *dmat) {
      PRECONDITION(bnd1, "");
      PRECONDITION(bnd2, "");
      PRECONDITION(bnd3, "");
      unsigned int bid1, bid2, bid3;
      bid1 = bnd1->getIdx();
      bid2 = bnd2->getIdx();
      bid3 = bnd3->getIdx();
      const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2)); 
      PRECONDITION(atm2, "");
      const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3)); 
      PRECONDITION(atm3, "");

      unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
      unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());

      // check that this actually is a 1-4 contact:
      if(dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<2.9){
        //std::cerr<<"skip: "<<aid1<<"-"<<aid4<<" because d="<<dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<<std::endl;
        return;
      }

      // when we have fused rings, it can happen that this isn't actually a 1-4 contact,
      // (this was the cause of sf.net bug 2835784) check that now:
      if(mol.getBondBetweenAtoms(aid1,atm3->getIdx()) ||
         mol.getBondBetweenAtoms(aid4,atm2->getIdx())) {
        return;
      }

      Atom::HybridizationType ahyb3 = atm3->getHybridization();
      Atom::HybridizationType ahyb2 = atm2->getHybridization(); 

      double bl1 = accumData.bondLengths[bid1];
      double bl2 = accumData.bondLengths[bid2];
      double bl3 = accumData.bondLengths[bid3];
      
      double ba12 = accumData.bondAngles->getVal(bid1, bid2);
      double ba23 = accumData.bondAngles->getVal(bid2, bid3);
      CHECK_INVARIANT(ba12 > 0.0, "");
      CHECK_INVARIANT(ba23 > 0.0, "");
      double dl, du;
      Path14Configuration path14;
      unsigned int nb = mol.getNumBonds();

      path14.bid1 = bid1; path14.bid2 = bid2; path14.bid3 = bid3;
      if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) { // FIX: check for trans
        // here we will assume 180 degrees: basically flat ring with an external substituent 
        dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
        du = dl; 
        dl -= GEN_DIST_TOL;
        du += GEN_DIST_TOL;
        path14.type = Path14Configuration::TRANS;
        accumData.transPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
        accumData.transPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
        
      } else {
        // here we will assume anything is possible
        dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
        du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);

        // in highly-strained situations these can get mixed up:
        if(du<dl){
          double tmpD=dl;
          dl=du;
          du=tmpD;
        }
        if (fabs(du-dl) < DIST12_DELTA) {
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;
        }
        path14.type = Path14Configuration::OTHER;
      }
      //std::cerr<<"1: "<<aid1<<"-"<<aid4<<": "<<dl<<" -> "<<du<<std::endl;
      _checkAndSetBounds(aid1, aid4, dl ,du, mmat);
      accumData.paths14.push_back(path14);
    }
    
    void _setTwoInDiffRing14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2, 
                                   const Bond *bnd3, ComputedData &accumData,
                                   DistGeom::BoundsMatPtr mmat,double *dmat) {
      // this turns out to be very similar to all bonds in the same ring situation.
      // There is probably some fine tuning that can be done when the atoms a2 and a3 are not sp2 hybridized,
      // but we will not worry about that now; simple use 0-180 deg for non-sp2 cases.
      _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat,dmat);
    }
    
    void _setShareRingBond14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2, 
                                   const Bond *bnd3, ComputedData &accumData,
                                   DistGeom::BoundsMatPtr mmat,double *dmat) {
      // once this turns out to be similar to bonds in the same ring
      _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat,dmat);
    }
    
    bool _checkH2NX3H1OX2(const Atom *atm) {
      if ((atm->getAtomicNum() == 6) && ( atm->getTotalNumHs() == 2) ) {
        // CH2
        return true;
      } else if ((atm->getAtomicNum() == 8) && (atm->getTotalNumHs() == 0) ) {
        // OX2
        return true;
      } else if ((atm->getAtomicNum() == 7) && (atm->getDegree() == 3) && 
                 (atm->getTotalNumHs() == 1)) { 
        // FIX: assumming hydrogen is not in the graph
        // this si NX3H1 situation
        return true;
      }
      return false;
    }
    
    bool _checkNhChChNh(const Atom *atm1, const Atom *atm2, const Atom *atm3, 
                        const Atom *atm4) {
      //checking for [!#1]~$ch!@$ch~[!#1], where ch = [CH2,NX3H1,OX2] situation
      if ((atm1->getAtomicNum() != 1) && (atm4->getAtomicNum() != 1) ) {
        // end atom not hydrogens
        if ( (_checkH2NX3H1OX2(atm2)) && (_checkH2NX3H1OX2(atm3)) ) {
          return true;
        }
      }
      return false;
    }
    
    // here we look for something like this:
    // It's an amide or ester:
    //
    //        4    <- 4 is the O
    //        |    <- That's the double bond
    //    1   3
    //     \ / \                                         T.S.I.Left Blank
    //      2   5  <- 2 is an oxygen/nitrogen
    bool _checkAmideEster14(const Bond *bnd1, const Bond *bnd3, const Atom *atm1,
			    const Atom *atm2, const Atom *atm3, const Atom *atm4) { 
      unsigned int a1Num = atm1->getAtomicNum();
      unsigned int a2Num = atm2->getAtomicNum();
      unsigned int a3Num = atm3->getAtomicNum();
      unsigned int a4Num = atm4->getAtomicNum();
      
      if ( a1Num != 1 && a3Num==6 &&
           bnd3->getBondType()==Bond::DOUBLE &&
           (a4Num==8 || a4Num==7) &&
           bnd1->getBondType()==Bond::SINGLE &&
           (a2Num==8 || (a2Num==7 && atm2->getTotalNumHs()==1)) ){
        return true;
      }
      return false;
    }

    bool _isCarbonyl(const ROMol &mol,const Atom *at){
      PRECONDITION(at,"bad atom");
      if(at->getAtomicNum()==6 && at->getDegree()>2 ){
	ROMol::ADJ_ITER nbrIdx,endNbrs;
	boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(at);
	while(nbrIdx!=endNbrs){
	  unsigned int atNum=mol.getAtomWithIdx(*nbrIdx)->getAtomicNum();
	  if( (atNum==8 || atNum==7) &&
	      mol.getBondBetweenAtoms(at->getIdx(),*nbrIdx)->getBondType()==Bond::DOUBLE) {
	    return true;
	  }
	  ++nbrIdx;
	}
      }
      return false;
    }

    bool _checkAmideEster15(const ROMol &mol,
			    const Bond *bnd1, const Bond *bnd3, const Atom *atm1,
			    const Atom *atm2, const Atom *atm3, const Atom *atm4) { 
      unsigned int a2Num = atm2->getAtomicNum();
      if ( (a2Num == 8) || 
           ((a2Num == 7) && (atm2->getTotalNumHs() == 1)) ) {
        if ((atm1->getAtomicNum() != 1) && (bnd1->getBondType() == Bond::SINGLE)) {
          if ((atm3->getAtomicNum() == 6) && 
              (bnd3->getBondType() == Bond::SINGLE) &&
	      _isCarbonyl(mol,atm3)
	      ) {
            return true;
          }
        }
      }
      return false;
    }


    void _setChain14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2, const Bond *bnd3, 
                           ComputedData &accumData,
                           DistGeom::BoundsMatPtr mmat,double *dmat){

      PRECONDITION(bnd1, "");
      PRECONDITION(bnd2, "");
      PRECONDITION(bnd3, "");
      unsigned int bid1, bid2, bid3;
      bid1 = bnd1->getIdx();
      bid2 = bnd2->getIdx();
      bid3 = bnd3->getIdx();
      const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
      PRECONDITION(atm2, "");
      const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3)); 
      PRECONDITION(atm3, "");
      
      unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
      unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
      const Atom *atm1 = mol.getAtomWithIdx(aid1);
      const Atom *atm4 = mol.getAtomWithIdx(aid4);
      
      double bl1 = accumData.bondLengths[bid1];
      double bl2 = accumData.bondLengths[bid2];
      double bl3 = accumData.bondLengths[bid3];
      
      double ba12 = accumData.bondAngles->getVal(bid1, bid2);
      double ba23 = accumData.bondAngles->getVal(bid2, bid3);
      CHECK_INVARIANT(ba12 > 0.0, "");
      CHECK_INVARIANT(ba23 > 0.0, "");
      bool setTheBound=true;
      double dl=0.0, du=0.0;
      
      // if the middle bond is double 
      Path14Configuration path14;
      path14.bid1 = bid1; path14.bid2 = bid2; path14.bid3 = bid3;
      unsigned int nb = mol.getNumBonds();
      switch(bnd2->getBondType()) {
      case Bond::DOUBLE :
        // if any of the other bonds are double - the torsion angle is zero
        // this is CC=C=C situation
        if ((bnd1->getBondType() == Bond::DOUBLE) || (bnd3->getBondType() == Bond::DOUBLE)) {
          dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
          du = dl + 2*GEN_DIST_TOL;
          path14.type = Path14Configuration::CIS;
          accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
          accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
          //BOOST_LOG(rdDebugLog) << "Special 5 " << aid1 << " " << aid4 << "\n";
        } else if (bnd2->getStereo() > Bond::STEREOANY) {
          Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
          if (stype == Bond::STEREOZ) {
            dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
            du = dl + 2*GEN_DIST_TOL;
            path14.type = Path14Configuration::CIS;
            //BOOST_LOG(rdDebugLog) << "Special 6 " <<  aid1 << " " << aid4 << "\n";
            accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
            accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
          } else {
            //BOOST_LOG(rdDebugLog) << "Special 7 " << aid1 << " " << aid4 << "\n";
            du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
            dl = du; 
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
            path14.type = Path14Configuration::TRANS;
            accumData.transPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
            accumData.transPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
          }
        } else {
          // double bond with no stereo setting can be 0 or 180
          dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
          du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
          if (fabs(du-dl) < DIST12_DELTA) {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
          }
          path14.type = Path14Configuration::OTHER;
        }
        break;
      case Bond::SINGLE :

        // Commenting out the following if block to fix issue 235, we may want to later provide 
        // the user with an option to invoke this special case  
#if 0
        if ( (_checkNhChChNh(atm1, atm2, atm3, atm4)) ||
             ((bnd1->getBondType() == Bond::DOUBLE) && (bnd3->getBondType() == Bond::DOUBLE) ) ) {
          // this is either 
          //  1. [!#1]~$ch!@$ch~[!#1] situation where ch = [CH2,NX3H1,OX2] or
          //  2. *=*-*=* situation
          // Both case cases we use 180 deg for torsion
          du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
          dl = du; 
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;
          path14.type = Path14Configuration::TRANS;
          transPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
          transPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
        } else 
#endif
        if ((atm2->getAtomicNum() == 16) && (atm3->getAtomicNum() == 16)) {
          // this is *S-S* situation
          //FIX: this cannot be right is sulfur has more than two coordinated
          // the torsion angle is 90 deg 
          dl = RDGeom::compute14Dist3D(bl1, bl2, bl3, ba12, ba23, M_PI/2) - GEN_DIST_TOL;
          du = dl + 2*GEN_DIST_TOL;
          path14.type = Path14Configuration::OTHER;
          //BOOST_LOG(rdDebugLog) << "Special 9 " << aid1 << " " << aid4 << "\n";
        } else if (( _checkAmideEster14(bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
                   ( _checkAmideEster14(bnd3, bnd1, atm4, atm3, atm2, atm1))) {

	  // It's an amide or ester:
	  //
	  //        4    <- 4 is the O
	  //        |    <- That's the double bond
	  //    1   3
	  //     \ / \                                         T.S.I.Left Blank
	  //      2   5  <- 2 is an oxygen/nitrogen
	  // 
	  // Here we set the distance between atoms 1 and 4,
	  //  we'll handle atoms 1 and 5 below.


          // fix for issue 251 - we were marking this as a cis configuration earlier
	  // -------------------------------------------------------
	  // Issue284:
	  //   As this code originally stood, we forced amide bonds to be trans. This is
	  //   convenient a lot of the time for generating nice-looking structures, but is
	  //   unfortunately totally bogus.  So here we'll allow the distance to
	  //   roam from cis to trans and hope that the force field planarizes things later.
	  //
	  //   What we'd really like to be able to do is specify multiple possible ranges
	  //   for the distances, but a single bounds matrix doesn't support this kind
	  //   of fanciness. 
	  //
#ifdef FORCE_TRANS_AMIDES	  
          dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23); 
          path14.type = Path14Configuration::TRANS;
          accumData.transPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
          accumData.transPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
#else
          if(atm2->getAtomicNum()==7 && atm2->getDegree()==3 &&
             atm1->getAtomicNum()==1 && atm2->getTotalNumHs(true)==1){
            // secondary amide, this is the H
            setTheBound=false;
          } else {
            // force the amide to be cis:
            dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
            path14.type = Path14Configuration::CIS;
            accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
            accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
          }
#endif
          du = dl;
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;

	  //BOOST_LOG(rdDebugLog) << "  amide: " << aid1 << " " << aid4 << ": " << dl << "->" << du << "\n";
        } else if (( _checkAmideEster15(mol, bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
                   ( _checkAmideEster15(mol, bnd3, bnd1, atm4, atm3, atm2, atm1))) {
	  // it's an amide or ester. 
	  //
	  //        4    <- 4 is the O
	  //        |    <- That's the double bond
	  //    1   3
	  //     \ / \                                          T.S.I.Left Blank
	  //      2   5  <- 2 is oxygen or nitrogen
	  // 
	  // we set the 1-4 contact above.

	  // If we're going to have a hope of getting good geometries
	  // out of here we need to set some reasonably smart bounds between 1
	  // and 5 (ref Issue355):

	  // NOTE THAT WE REVERSE THE ORDER HERE:

#ifdef FORCE_TRANS_AMIDES
	  // amide is trans, we're cis:
          dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
          path14.type = Path14Configuration::CIS;
          accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
          accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
#else
	  // amide is cis, we're trans:
          if(atm2->getAtomicNum()==7 && atm2->getDegree()==3 &&
             atm1->getAtomicNum()==1 && atm2->getTotalNumHs(true)==1){
            // secondary amide, this is the H
            setTheBound=false;
          } else {
            dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23); 
            path14.type = Path14Configuration::TRANS;
            accumData.transPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
            accumData.transPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
          }
#endif
          du = dl;
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;
	  //BOOST_LOG(rdDebugLog) << "    amide neighbor: " << aid1 << " " << aid4 << ": " << dl << "->" << du << "\n";
        } else {
          dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
          du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
          path14.type = Path14Configuration::OTHER;
        }
        break;
      default:
        //BOOST_LOG(rdDebugLog) << "Special 12 " << aid1 << " " << aid4 << "\n";
        dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23); 
        du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
        
        path14.type = Path14Configuration::OTHER;
      }
      if(setTheBound){
        if (fabs(du-dl) < DIST12_DELTA) {
          dl -= GEN_DIST_TOL;
          du += GEN_DIST_TOL;
        }
        //std::cerr<<"2: "<<aid1<<"-"<<aid4<<std::endl;
        _checkAndSetBounds(aid1, aid4, dl, du, mmat);
        accumData.paths14.push_back(path14);
      }
    }  
    
    void _record14Path(const ROMol &mol, unsigned int bid1, unsigned int bid2,
                        unsigned int bid3, ComputedData &accumData) {
      const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2)); 
      PRECONDITION(atm2, "");
      Atom::HybridizationType ahyb2 = atm2->getHybridization(); 
      const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3)); 
      PRECONDITION(atm3, "");
      Atom::HybridizationType ahyb3 = atm3->getHybridization();
      unsigned int nb = mol.getNumBonds();
      Path14Configuration path14;
      path14.bid1 = bid1; path14.bid2 = bid2; path14.bid3 = bid3;
      if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) { // FIX: check for trans
        path14.type = Path14Configuration::CIS;
        accumData.cisPaths[bid1*nb*nb + bid2*nb + bid3] = 1;
        accumData.cisPaths[bid3*nb*nb + bid2*nb + bid1] = 1;
      } else {
        path14.type = Path14Configuration::OTHER;
      }
      accumData.paths14.push_back(path14);
    }

    void set14Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, 
                     ComputedData &accumData,double *distMatrix) {
      unsigned int npt = mmat->numRows();
      CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
      
      const RingInfo *rinfo = mol.getRingInfo(); // FIX: make sure we have ring info
      CHECK_INVARIANT(rinfo, "");
      const VECT_INT_VECT &bondRings = rinfo->bondRings();
      VECT_INT_VECT_CI rii;
      unsigned int i, aid2, aid3;
      unsigned int bid1, bid2, bid3;
      ROMol::OEDGE_ITER beg1, beg2, end1, end2;
      ROMol::ConstBondIterator bi;
      
      unsigned int nb = mol.getNumBonds();
      BIT_SET ringBondPairs(nb*nb), donePaths(nb*nb*nb);
      unsigned int id1, id2, pid1, pid2, pid3, pid4;
      // first we will deal with 1-4 atoms that belong to the same ring    
      for (rii = bondRings.begin(); rii != bondRings.end(); rii++) {
        // we don't need deal with 3 membered rings
        unsigned int rSize = rii->size();

        bid1 = (*rii)[rSize-1];
        for (i  = 0; i < rSize; i++) {
          bid2 = (*rii)[i];
          if (i == rSize-1){
            bid3 = (*rii)[0];
          } else {
            bid3 = (*rii)[i+1];
          }
          pid1 = bid1*nb + bid2;
          pid2 = bid2*nb + bid1;
          id1 = bid1*nb*nb + bid2*nb + bid3;
          id2 = bid3*nb*nb + bid2*nb + bid1;
          
          ringBondPairs[pid1] = 1;
          ringBondPairs[pid2] = 1;
          donePaths[id1] = 1;
          donePaths[id2] = 1;
          
          if (rSize > 5) {
            _setInRing14Bounds(mol, mol.getBondWithIdx(bid1),
                               mol.getBondWithIdx(bid2), mol.getBondWithIdx(bid3),
                               accumData, mmat, distMatrix);
          } else {
            _record14Path(mol, bid1, bid2, bid3, accumData);
          }
         
          bid1 = bid2;
        } // loop over bonds in the ring
      } // end of all rings
      for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
        bid2 = (*bi)->getIdx();
        aid2 = (*bi)->getBeginAtomIdx();
        aid3 = (*bi)->getEndAtomIdx();
        boost::tie(beg1,end1) = mol.getAtomBonds(mol.getAtomWithIdx(aid2));
        while (beg1 != end1) {
          const Bond *bnd1 = mol[*beg1].get();
          bid1 = bnd1->getIdx();
          if (bid1 != bid2) {
            boost::tie(beg2,end2) = mol.getAtomBonds(mol.getAtomWithIdx(aid3));
            while (beg2 != end2) {
              const Bond *bnd3 = mol[*beg2].get();
              bid3 = bnd3->getIdx();
              if (bid3 != bid2) {
                id1 = nb*nb*bid1 + nb*bid2 + bid3;
                id2 = nb*nb*bid3 + nb*bid2 + bid1;
                if ((!donePaths[id1]) && (!donePaths[id2])) {
                  // we haven't dealt with this path before
                  pid1 = bid1*nb + bid2;
                  pid2 = bid2*nb + bid1;
                  pid3 = bid2*nb + bid3;
                  pid4 = bid3*nb + bid2;

                  if (ringBondPairs[pid1] || ringBondPairs[pid2] 
                      || ringBondPairs[pid3] || ringBondPairs[pid4]) {
                    // either (bid1, bid2) or (bid2, bid3) are in the
		    // same ring (note all three cannot be in the same
		    // ring; we dealt with that before)
                    _setTwoInSameRing14Bounds(mol, bnd1, (*bi), bnd3, accumData, mmat, distMatrix);
                  } else if ( ((rinfo->numBondRings(bid1) > 0) && (rinfo->numBondRings(bid2) > 0)) ||
                              ((rinfo->numBondRings(bid2) > 0) && (rinfo->numBondRings(bid3) > 0)) ) {
                    // (bid1, bid2) or (bid2, bid3) are ring bonds but
                    // belong to different rings.  Note that the third
                    // bond will not belong to either of these two
                    // rings (if it does, we would have taken care of
                    // it in the previous if block); i.e. if bid1 and
                    // bid2 are ring bonds that belong to ring r1 and
                    // r2, then bid3 is either an external bond or
                    // belongs to a third ring r3.
                    _setTwoInDiffRing14Bounds(mol, bnd1, (*bi), bnd3, accumData, mmat, distMatrix);
                  } else if (rinfo->numBondRings(bid2) > 0) {
                    // the middle bond is a ring bond and the other
                    // two do not belong to the same ring or are
                    // non-ring bonds
                    _setShareRingBond14Bounds(mol, bnd1, (*bi), bnd3, accumData, mmat, distMatrix);
                  } else { 
                    // middle bond not a ring
                    _setChain14Bounds(mol, bnd1, (*bi), bnd3, accumData, mmat, distMatrix);

                  }
                }
              }
              ++beg2;
            }
          }
          ++beg1;
        }
      }
    }  


    void initBoundsMat(DistGeom::BoundsMatrix *mmat,double defaultMin,
		       double defaultMax){
      unsigned int npt = mmat->numRows();
      
      for (unsigned int i = 1; i < npt; i++) {
	for (unsigned int j = 0; j < i; j++) {
	  mmat->setUpperBound(i,j,defaultMax);
	  mmat->setLowerBound(i,j,defaultMin);
	}
      }
    }
    void initBoundsMat(DistGeom::BoundsMatPtr mmat,double defaultMin,
		       double defaultMax){
      initBoundsMat(mmat.get(),defaultMin,defaultMax);
    };

    void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                        bool set15bounds, bool scaleVDW) {
      PRECONDITION(mmat.get(),"bad pointer");
      unsigned int nb = mol.getNumBonds();
      unsigned int na = mol.getNumAtoms();
      if(!na){
        throw ValueErrorException("molecule has no atoms");
      }
      ComputedData accumData(na, nb);
      double *distMatrix=0;
      distMatrix = MolOps::getDistanceMat(mol);

      set12Bounds(mol, mmat, accumData);
      set13Bounds(mol, mmat, accumData);

      set14Bounds(mol, mmat, accumData,distMatrix);

      if (set15bounds) {
        set15Bounds(mol, mmat, accumData,distMatrix);
      }

      setLowerBoundVDW(mol, mmat, scaleVDW,distMatrix);
    }


    // some helper functions to set 15 distances
    
    
    /*
     compute the lower and upper bounds for the distance between 15 atoms give than the first
     four atoms are in cis configuration. The 15 limits are computed assuming the following
     configuration
             5
              \  
         1     4
          \   / 
           2-3  
     ARGUMENTS:
       d1 - distance between 1 and 2
       d2 - distance between 2 and 3
       d3 - distance between 3 and 4
       d4 - distance between 4 and 5
       ang12 - angle(123)
       ang23 - angle(234)
       and34 - angle(345)
       dl - storage for lower 15 bound 
       du - storage for upper 15 bound
    */
    double _compute15DistsCisCis(double d1, double d2, double d3, double d4, 
                                 double ang12, double ang23, double ang34) {
      double dx14 = d2 - d3*cos(ang23) - d1*cos(ang12);
      double dy14 = d3*sin(ang23) - d1*sin(ang12);
      double d14 = sqrt(dx14*dx14 + dy14*dy14);
      double cval = (d3 - d2*cos(ang23) + d1*cos(ang12 + ang23))/d14;
      if (cval > 1.0) {
        cval = 1.0;
      } else if (cval < -1.0) {
        cval = -1.0;
      }

      double ang143 = acos(cval);
      double ang145 = ang34 - ang143;
      double res = RDGeom::compute13Dist(d14, d4, ang145);
      return res;
    }
    
    /*
     compute the lower and upper bounds for the distance between 15 atoms give than the first
     four atoms are in cis configuration. The 15 limits are computed assuming the following 
     configuration
                 
                 
      1     4-5  
       \   /     
        2-3      
    
     ARGUMENTS:
       d1 - distance between 1 and 2
       d2 - distance between 2 and 3
       d3 - distance between 3 and 4
       d4 - distance between 4 and 5
       ang12 - angle(123)
       ang23 - angle(234)
       and34 - angle(345)
       dl - storage for lower 15 bound 
       du - storage for upper 15 bound
    */
    double _compute15DistsCisTrans(double d1, double d2, double d3, double d4, 
                                   double ang12, double ang23, double ang34) {
      double dx14 = d2 - d3*cos(ang23) - d1*cos(ang12);
      double dy14 = d3*sin(ang23) - d1*sin(ang12);
      double d14 = sqrt(dx14*dx14 + dy14*dy14);
      double cval = (d3 - d2*cos(ang23) + d1*cos(ang12 + ang23))/d14;
      if (cval > 1.0) {
        cval = 1.0;
      } else if (cval < -1.0) {
        cval = -1.0;
      }

      double ang143 = acos(cval);
      double ang145 = ang34 + ang143;
      return RDGeom::compute13Dist(d14, d4, ang145);
    }

    /*
     compute the lower and upper bounds for the dsitance between 15 atoms given than the first
     four atoms are in trans configuration. The 15 limits are computed assuming the following 
     configuration
                            
      1         
       \        
        2-3     
           \    
            4-5 
                
                
     ARGUMENTS:
       d1 - distance between 1 and 2
       d2 - distance between 2 and 3
       d3 - distance between 3 and 4
       d4 - distance between 4 and 5
       ang12 - angle(123)
       ang23 - angle(234)
       and34 - angle(345)
       dl - storage for lower 15 bound 
       du - storage for upper 15 bound
    */
    double _compute15DistsTransTrans(double d1, double d2, double d3, double d4, double ang12, 
                                     double ang23, double ang34) {
      double dx14 = d2 - d3*cos(ang23) - d1*cos(ang12);
      double dy14 = d3*sin(ang23) + d1*sin(ang12);
      double d14 = sqrt(dx14*dx14 + dy14*dy14);
      double cval = (d3 - d2*cos(ang23) + d1*cos(ang12 - ang23))/d14;
      if (cval > 1.0) {
        cval = 1.0;
      } else if (cval < -1.0) {
        cval = -1.0;
      }

      double ang143 = acos(cval);
      double ang145 = ang34 + ang143;
      return RDGeom::compute13Dist(d14, d4, ang145);
    }
    
    /*
     compute the lower and upper bounds for the dsitance between 15 atoms given than the first
     four atoms are in trans configuration. The 15 limits are computed assuming the following 
     configuration
                            
                        1    
                         \   
                          2-3
                             \
                              4
                             /
                            5
     ARGUMENTS:
       d1 - distance between 1 and 2
       d2 - distance between 2 and 3
       d3 - distance between 3 and 4
       d4 - distance between 4 and 5
       ang12 - angle(123)
       ang23 - angle(234)
       and34 - angle(345)
       dl - storage for lower 15 bound 
       du - storage for upper 15 bound
    */
    double _compute15DistsTransCis(double d1, double d2, double d3, double d4, 
                                   double ang12, double ang23, double ang34) {
      double dx14 = d2 - d3*cos(ang23) - d1*cos(ang12);
      double dy14 = d3*sin(ang23) + d1*sin(ang12);
      double d14 = sqrt(dx14*dx14 + dy14*dy14);
      
      double cval = (d3 - d2*cos(ang23) + d1*cos(ang12 - ang23))/d14;
      if (cval > 1.0) {
        cval = 1.0;
      } else if (cval < -1.0) {
        cval = -1.0;
      }

      double ang143 = acos(cval);
      double ang145 = ang34 - ang143;
      return RDGeom::compute13Dist(d14, d4, ang145);
    }

    void _set15BoundsHelper(const ROMol &mol, unsigned int bid1, unsigned int bid2, unsigned int bid3,
			    unsigned int type, ComputedData &accumData,
                            DistGeom::BoundsMatPtr mmat, double *dmat) {
      unsigned int i, aid1, aid2, aid3, aid4, aid5;
      double d1, d2, d3, d4, ang12, ang23, ang34, du, dl, vw1, vw5;
      unsigned int nb = mol.getNumBonds();
      unsigned int na = mol.getNumAtoms();
      
      aid2 = accumData.bondAdj->getVal(bid1, bid2);
      aid1 = mol.getBondWithIdx(bid1)->getOtherAtomIdx(aid2);
      aid3 = accumData.bondAdj->getVal(bid2, bid3);
      aid4 = mol.getBondWithIdx(bid3)->getOtherAtomIdx(aid3);
      d1 = accumData.bondLengths[bid1];
      d2 = accumData.bondLengths[bid2];
      d3 = accumData.bondLengths[bid3];
      ang12 = accumData.bondAngles->getVal(bid1, bid2);
      ang23 = accumData.bondAngles->getVal(bid2, bid3);
      for (i = 0; i < nb; i++) {
        du = -1.0;
        dl = 0.0;
        if (accumData.bondAdj->getVal(bid3, i) == static_cast<int>(aid4)) {
          aid5 =  mol.getBondWithIdx(i)->getOtherAtomIdx(aid4);
          // make sure we did not com back to the first atom in the path - possible with 4 membered rings
          // this is a fix for Issue 244

          // check that this actually is a 1-5 contact:
          if(dmat[std::max(aid1,aid5)*mmat->numRows()+std::min(aid1,aid5)]<3.9){
            //std::cerr<<"skip: "<<aid1<<"-"<<aid5<<" because d="<<dmat[std::max(aid1,aid5)*mmat->numRows()+std::min(aid1,aid5)]<<std::endl;
            continue;
          }

          if (aid1 != aid5) { //FIX: do we need this
            unsigned int pid1 = aid1*na + aid5;
            unsigned int pid2 = aid5*na + aid1;
            if ((mmat->getLowerBound(aid1, aid5) < DIST12_DELTA) || 
                (accumData.set15Atoms[pid1]) || (accumData.set15Atoms[pid2])) {
              d4 = accumData.bondLengths[i];
              ang34 = accumData.bondAngles->getVal(bid3, i);
              unsigned int pathId = (bid2)*nb*nb + (bid3)*nb + i;
              if (type == 0) {
                if (accumData.cisPaths[pathId]) {
                  dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34);
                  du = dl + DIST15_TOL;
                  dl -= DIST15_TOL; 
                } else if (accumData.transPaths[pathId]) {
                  dl = _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34);
                  du = dl + DIST15_TOL;
                  dl -= DIST15_TOL;
                } else {
                  dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                  du = _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                } 
                
              } else if (type == 1) {
                if (accumData.cisPaths[pathId]) {
                  dl = _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34);
                  du = dl + DIST15_TOL;
                  dl -= DIST15_TOL;
                } else if (accumData.transPaths[pathId]) {
                  dl = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23, ang34);
                  du = dl + DIST15_TOL;
                  dl -= DIST15_TOL;
                } else {
                  dl = _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                  du = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                }
              } else {
                if (accumData.cisPaths[pathId]) {
                  dl = _compute15DistsCisCis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                  du = _compute15DistsCisTrans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                } else if (accumData.transPaths[pathId]) {
                  dl = _compute15DistsTransCis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                  du = _compute15DistsTransTrans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                } else {
                  vw1 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(aid1)->getAtomicNum());
                  vw5 = PeriodicTable::getTable()->getRvdw(mol.getAtomWithIdx(aid5)->getAtomicNum());
                  dl = VDW_SCALE_15*(vw1+vw5);
                }
              }
              if (du < 0.0) {
                du = MAX_UPPER;
              }
            
              //std::cerr<<"3: "<<aid1<<"-"<<aid5<<std::endl;
              _checkAndSetBounds(aid1, aid5, dl, du, mmat);
              accumData.set15Atoms[aid1*na + aid5] = 1;
              accumData.set15Atoms[aid5*na + aid1] = 1;
            }
          }
        }
      }
    }

    // set the 15 distance bounds
    void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat, 
                     ComputedData &accumData,double *distMatrix) {
      PATH14_VECT_CI pti;
      unsigned int bid1, bid2, bid3, type;
      for (pti = accumData.paths14.begin(); pti != accumData.paths14.end(); pti++) {
        bid1 = pti->bid1;
        bid2 = pti->bid2;
        bid3 = pti->bid3;
        type = pti->type;
        // 15 distances going one way with with 14 paths
        _set15BoundsHelper(mol, bid1, bid2, bid3, type, accumData, mmat, distMatrix);
        // goign the other way - reverse the 14 path
        _set15BoundsHelper(mol, bid3, bid2, bid1, type, accumData, mmat, distMatrix);
      }
    }
  }
}
