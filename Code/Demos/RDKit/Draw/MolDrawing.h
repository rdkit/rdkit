// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_MOLDRAWING_H_
#define _RD_MOLDRAWING_H_

#include <vector>
#include <boost/foreach.hpp>
#include <GraphMol/RDKitBase.h>
#include <Geometry/point.h>

/***********
  Return Format: vector of ints

  RESOLUTION dots_per_angstrom
  BOUNDS x1 y1 x2 y2
  LINE width dashed atom1_atnum atom2_atnum x1 y1 x2 y2
  WEDGE dashed atom1_atnum atom2_atnum x1 y1 x2 y2 x3 y3 
  ATOM idx atnum x y num_chars char1-char x 

*************/

namespace RDKit {
  namespace Drawing {
    typedef int ElementType;

    typedef enum {
      LINE=1,
      WEDGE,
      ATOM,
      BOUNDS,
      RESOLUTION
    } PrimType;

    std::vector<ElementType> DrawMol(const ROMol &mol,int confId=-1,unsigned int dotsPerAngstrom=100,
                                     double dblBondOffset=0.2,double dblBondLengthFrac=0.8){
      std::vector<ElementType> res;
      res.push_back(RESOLUTION);
      res.push_back(static_cast<ElementType>(dotsPerAngstrom));
      
      const Conformer &conf=mol.getConformer(confId);
      const RDGeom::POINT3D_VECT &locs=conf.getPositions();

      //------------
      // do the bounding box
      //------------
      double minx=1e6,miny=1e6,maxx=-1e6,maxy=-1e6;
      BOOST_FOREACH( const RDGeom::Point3D &pt,locs){
        minx=std::min(pt.x,minx);
        miny=std::min(pt.y,miny);
        maxx=std::max(pt.x,maxx);
        maxy=std::max(pt.y,maxy);
      }
      double dimx=(maxx-minx),dimy=(maxy-miny);
      res.push_back(BOUNDS);
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*0));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*0));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*dimx));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*dimy));

      // loop over atoms:
      ROMol::VERTEX_ITER bAts,eAts;
      boost::tie(bAts,eAts)=mol.getVertices();
      while(bAts!=eAts){
        int a1Idx=mol[*bAts]->getIdx();
        RDGeom::Point2D a1(locs[a1Idx].x-minx,locs[a1Idx].y-miny);
        ROMol::OEDGE_ITER nbr,endNbrs;
        boost::tie(nbr,endNbrs) = mol.getAtomBonds(mol[*bAts].get());
        while(nbr!=endNbrs){
          const BOND_SPTR bond=mol[*nbr];
          ++nbr;
          int a2Idx=bond->getOtherAtomIdx(a1Idx);
          if(a2Idx>a1Idx) continue;
          RDGeom::Point2D a2(locs[a2Idx].x-minx,locs[a2Idx].y-miny);
          
          res.push_back(LINE);
          res.push_back(1);
          res.push_back(0);
          res.push_back(mol[*bAts]->getAtomicNum());
          res.push_back(mol.getAtomWithIdx(a2Idx)->getAtomicNum());
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.x));
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.y));
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a2.x));
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a2.y));

          if(bond->getBondType()==Bond::DOUBLE || bond->getBondType()==Bond::TRIPLE ){
            RDGeom::Point2D obv=a2-a1;
            RDGeom::Point2D perp=obv;
            perp.rotate90();
            perp.normalize();

            if(bond->getBondType()==Bond::DOUBLE && mol.getRingInfo()->numBondRings(bond->getIdx())){
              // we're in a ring... we might need to flip sides:
              ROMol::OEDGE_ITER nbr2,endNbrs2;
              boost::tie(nbr2,endNbrs2) = mol.getAtomBonds(mol[*bAts].get());
              while(nbr2!=endNbrs2){
                const BOND_SPTR bond2=mol[*nbr2];
                ++nbr2;
                if(bond2->getIdx()==bond->getIdx() ||
                   !mol.getRingInfo()->numBondRings(bond2->getIdx())) continue;
                bool sharedRing=false;
                BOOST_FOREACH(const INT_VECT &ring,mol.getRingInfo()->bondRings()){
                  if(std::find(ring.begin(),ring.end(),bond->getIdx())!=ring.end() &&
                     std::find(ring.begin(),ring.end(),bond2->getIdx())!=ring.end()){
                    sharedRing=true;
                    break;
                  }
                }
                if(sharedRing){
                  // these two bonds share a ring.
                  int a3Idx=bond2->getOtherAtomIdx(a1Idx);
                  if(a3Idx!=a2Idx){
                    RDGeom::Point2D a3(locs[a3Idx].x-minx,locs[a3Idx].y-miny);
                    RDGeom::Point2D obv2=a3-a1;
                    if(obv2.dotProduct(perp)<0){
                      perp*=-1;
                    }
                  }
                }
              }
            }
            perp *= dblBondOffset;

            RDGeom::Point2D offsetStart=a1 + obv*(.5*(1.-dblBondLengthFrac));

            obv *= dblBondLengthFrac;

            res.push_back(LINE);
            res.push_back(1);
            res.push_back(0);
            res.push_back(mol[*bAts]->getAtomicNum());
            res.push_back(mol.getAtomWithIdx(a2Idx)->getAtomicNum());
            res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.x+perp.x)));
            res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.y+perp.y)));
            res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.x+obv.x+perp.x)));
            res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.y+obv.y+perp.y)));
            
            if(bond->getBondType()==Bond::TRIPLE){
              res.push_back(LINE);
              res.push_back(1);
              res.push_back(0);
              res.push_back(mol[*bAts]->getAtomicNum());
              res.push_back(mol.getAtomWithIdx(a2Idx)->getAtomicNum());
              res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.x-perp.x)));
              res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.y-perp.y)));
              res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.x+obv.x-perp.x)));
              res.push_back(static_cast<ElementType>(dotsPerAngstrom*(offsetStart.y+obv.y-perp.y)));
            }
          }
        }
        res.push_back(ATOM);
        res.push_back(mol[*bAts]->getAtomicNum());
        res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.x));
        res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.y));
        std::string symbol=mol[*bAts]->getSymbol();
        res.push_back(static_cast<ElementType>(symbol.length()));
        BOOST_FOREACH(char c, symbol){
          res.push_back(static_cast<ElementType>(c));
        }
        
        ++bAts;
      }
      
      return res;
    }
  }   // end of namespace Drawing
} // end of namespace RDKit

#endif
