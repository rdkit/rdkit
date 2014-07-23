// $Id$
//
//  Copyright (C) 2009-2012 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Includes contributions from Dave Cosgrove (davidacosgroveaz@gmail.com)
//
#ifndef _RD_MOLDRAWING_H_
#define _RD_MOLDRAWING_H_

#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <GraphMol/RDKitBase.h>
#include <Geometry/point.h>
#include <GraphMol/Depictor/RDDepictor.h>

/***********
  Return Format: vector of ints

  RESOLUTION dots_per_angstrom
  BOUNDS x1 y1 x2 y2
  LINE width dashed atom1_atnum atom2_atnum x1 y1 x2 y2
  WEDGE dashed atom1_atnum atom2_atnum x1 y1 x2 y2 x3 y3 
  ATOM idx atnum x y num_chars char1-charx orient



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
    typedef enum {
      C=0,
      N,
      E,
      S,
      W
    } OrientType;

    namespace detail {
      // **************************************************************************
      void drawLine( std::vector<ElementType> &res ,
                     int atnum1 , int atnum2 , int lineWidth , int dashed ,
                     double x1 , double y1 ,
                     double x2 , double y2 ) {

        res.push_back( LINE );
        res.push_back( static_cast<ElementType>(lineWidth) );
        res.push_back(dashed);
        res.push_back( static_cast<ElementType>(atnum1) );
        res.push_back( static_cast<ElementType>(atnum2) );
        res.push_back( static_cast<ElementType>(x1) );
        res.push_back( static_cast<ElementType>(y1) );
        res.push_back( static_cast<ElementType>(x2) );
        res.push_back( static_cast<ElementType>(y2) );

      }
      std::pair<std::string,OrientType> getAtomSymbolAndOrientation(const Atom &atom,RDGeom::Point2D nbrSum){
        std::string symbol="";
        OrientType orient=C;
        int isotope = atom.getIsotope();
        if(atom.getAtomicNum()!=6 ||
           atom.getFormalCharge()!=0 ||
           isotope ||
           atom.getNumRadicalElectrons()!=0 ||
           atom.hasProp("molAtomMapNumber") ||
           atom.getDegree()==0 ){
          symbol=atom.getSymbol();
          bool leftToRight=true;
          if(atom.getDegree()==1 && nbrSum.x>0){
            leftToRight=false;
          }
          if(isotope){
            symbol = boost::lexical_cast<std::string>(isotope)+symbol;
          }
          if(atom.hasProp("molAtomMapNumber")){
            std::string mapNum;
            atom.getProp("molAtomMapNumber",mapNum);
            symbol += ":" + mapNum;
          }
          int nHs=atom.getTotalNumHs();
          if(nHs>0){
            std::string h="H";
            if(nHs>1) {
              h += boost::lexical_cast<std::string>(nHs);
            }
            if(leftToRight) symbol += h;
            else symbol = h+symbol;
          }
          if( atom.getFormalCharge()!=0 ){
            int chg=atom.getFormalCharge();
            std::string sgn="+";
            if(chg<0){
              sgn="-";
            }
            chg=abs(chg);
            if(chg>1){
              sgn += boost::lexical_cast<std::string>(chg);
            } 
            if(leftToRight) symbol+=sgn;
            else symbol = sgn+symbol;
          }

          if(atom.getDegree()==1){
            double islope=0;
            if(fabs(nbrSum.y)>1){
              islope=nbrSum.x/fabs(nbrSum.y);
            } else {
              islope=nbrSum.x;
            }
            if(fabs(islope)>.85){
              if(islope>0){
                orient=W;
              } else {
                orient=E;
              }
            } else {
              if(nbrSum.y>0){
                orient=N;
              } else {
                orient=S;
              }
            }
          }
        }
        return std::make_pair(symbol,orient);
      }
    } // end of detail namespace
    // **************************************************************************
    std::vector<ElementType> DrawMol(const ROMol &mol,int confId=-1,
                                     const std::vector<int> *highlightAtoms=0,
                                     unsigned int dotsPerAngstrom=100,
                                     double dblBondOffset=0.3,
                                     double dblBondLengthFrac=0.8,
                                     double angstromsPerChar=0.20){
      if(!mol.getRingInfo()->isInitialized()){
        MolOps::findSSSR(mol);
      }
      std::vector<ElementType> res;
      res.push_back(RESOLUTION);
      res.push_back(static_cast<ElementType>(dotsPerAngstrom));
      
      const Conformer &conf=mol.getConformer(confId);
      const RDGeom::POINT3D_VECT &locs=conf.getPositions();

      // get atom symbols and orientations
      // (we need them for the bounding box calculation)
      std::vector< std::pair<std::string,OrientType> > atomSymbols;
      ROMol::VERTEX_ITER bAts,eAts;
      boost::tie(bAts,eAts)=mol.getVertices();
      while(bAts!=eAts){
        ROMol::OEDGE_ITER nbr,endNbrs;
        RDGeom::Point2D nbrSum(0,0);
        boost::tie(nbr,endNbrs) = mol.getAtomBonds(mol[*bAts].get());
        RDGeom::Point2D a1(locs[mol[*bAts]->getIdx()].x,
                           locs[mol[*bAts]->getIdx()].y);
        while(nbr!=endNbrs){
          const BOND_SPTR bond=mol[*nbr];
          ++nbr;
          int a2Idx=bond->getOtherAtomIdx(mol[*bAts]->getIdx());
          RDGeom::Point2D a2(locs[a2Idx].x,locs[a2Idx].y);
          nbrSum+=a2-a1;
        }
        atomSymbols.push_back(detail::getAtomSymbolAndOrientation(*mol[*bAts],nbrSum));
        ++bAts;
      }


      //------------
      // do the bounding box
      //------------
      double minx=1e6,miny=1e6,maxx=-1e6,maxy=-1e6;
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        RDGeom::Point3D pt=locs[i];
        std::string symbol;
        OrientType orient;
        boost::tie(symbol,orient)=atomSymbols[i];
        if(symbol!=""){
          // account for a possible expansion of the bounding box by the symbol
          if(pt.x<=minx){
            switch(orient){
            case C:
            case N:
            case S:
            case E:
              minx = pt.x-symbol.size()/2*angstromsPerChar;
              break;
            case W:
              minx = pt.x-symbol.size()*angstromsPerChar;
              break;
            }
          }
          if(pt.x>=maxx){
            switch(orient){
            case C:
            case N:
            case S:
            case W:
              maxx = pt.x+symbol.size()/2*angstromsPerChar;
              break;
            case E:
              maxx = pt.x+symbol.size()*angstromsPerChar;
              break;
            }
          }

          if(pt.y<=miny){
            miny = pt.y-1.5*angstromsPerChar;
          }
          if(pt.y>=maxy){
            maxy = pt.y+angstromsPerChar;
          }
        } else {
          minx=std::min(pt.x,minx);
          miny=std::min(pt.y,miny);
          maxx=std::max(pt.x,maxx);
          maxy=std::max(pt.y,maxy);
        }
      }
      double dimx=(maxx-minx),dimy=(maxy-miny);
      res.push_back(BOUNDS);
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*0));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*0));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*dimx));
      res.push_back(static_cast<ElementType>(dotsPerAngstrom*dimy));

      // loop over atoms:
      boost::tie(bAts,eAts)=mol.getVertices();
      while(bAts!=eAts){
        int a1Idx=mol[*bAts]->getIdx();
        RDGeom::Point2D a1(locs[a1Idx].x-minx,locs[a1Idx].y-miny);
        ROMol::OEDGE_ITER nbr,endNbrs;
        RDGeom::Point2D nbrSum(0,0);
        boost::tie(nbr,endNbrs) = mol.getAtomBonds(mol[*bAts].get());
        while(nbr!=endNbrs){
          const BOND_SPTR bond=mol[*nbr];
          ++nbr;
          int a2Idx=bond->getOtherAtomIdx(a1Idx);
          int lineWidth=1;
          if(highlightAtoms
             &&
             std::find(highlightAtoms->begin(),highlightAtoms->end(),a1Idx)
              != highlightAtoms->end()
             &&
             std::find(highlightAtoms->begin(),highlightAtoms->end(),a2Idx)
              != highlightAtoms->end() ){
            lineWidth=3;
          }
          RDGeom::Point2D a2(locs[a2Idx].x-minx,locs[a2Idx].y-miny);
          nbrSum+=a2-a1;
          if(a2Idx<a1Idx) continue;

          // draw bond from a1 to a2.
          int atnum1 = mol[*bAts]->getAtomicNum();
          int atnum2 = mol.getAtomWithIdx(a2Idx)->getAtomicNum();

          if( !mol.getRingInfo()->numBondRings(bond->getIdx()) && bond->getBondType()!=Bond::AROMATIC ) {
            // acyclic bonds
            RDGeom::Point2D obv=a2-a1;
            RDGeom::Point2D perp=obv;
            perp.rotate90();
            perp.normalize();

            if( bond->getBondType()==Bond::DOUBLE || bond->getBondType()==Bond::TRIPLE)  {
              RDGeom::Point2D startP=a1,endP=a2;              
              if( bond->getBondType()==Bond::TRIPLE){
                perp *= dblBondOffset;
                startP+=(obv*(1.-dblBondLengthFrac)/2);
              endP-=(obv*(1.-dblBondLengthFrac)/2);
              } else {
                perp *= 0.5 * dblBondOffset;
              }
              detail::drawLine( res , atnum1 , atnum2 , lineWidth , 0 ,
                                dotsPerAngstrom*(startP.x+perp.x) ,
                                dotsPerAngstrom*(startP.y+perp.y) ,
                                dotsPerAngstrom*(endP.x+perp.x) ,
                                dotsPerAngstrom*(endP.y+perp.y) );
              if(bond->getBondType() != Bond::AROMATIC){
                detail::drawLine( res , atnum1 , atnum2 , lineWidth , 0 ,
                                  dotsPerAngstrom*(startP.x-perp.x) ,
                                  dotsPerAngstrom*(startP.y-perp.y) ,
                                  dotsPerAngstrom*(endP.x-perp.x) ,
                                  dotsPerAngstrom*(endP.y-perp.y) );
              } else {
                detail::drawLine( res , atnum1 , atnum2 , lineWidth , 1 ,
                                  dotsPerAngstrom*(startP.x-perp.x) ,
                                  dotsPerAngstrom*(startP.y-perp.y) ,
                                  dotsPerAngstrom*(endP.x-perp.x) ,
                                  dotsPerAngstrom*(endP.y-perp.y) );

              }
            }
            if( bond->getBondType()==Bond::SINGLE || bond->getBondType()==Bond::TRIPLE ) {
              detail::drawLine( res , atnum1 , atnum2 , lineWidth , 0 ,
                                dotsPerAngstrom*(a1.x) ,
                                dotsPerAngstrom*(a1.y) ,
                                dotsPerAngstrom*(a2.x) ,
                                dotsPerAngstrom*(a2.y) );
            } else if( bond->getBondType()!=Bond::DOUBLE ) {
              detail::drawLine( res , atnum1 , atnum2 , lineWidth , 2 ,
                                dotsPerAngstrom*(a1.x) ,
                                dotsPerAngstrom*(a1.y) ,
                                dotsPerAngstrom*(a2.x) ,
                                dotsPerAngstrom*(a2.y) );
            }
          } else {
            // cyclic bonds
              detail::drawLine( res , atnum1 , atnum2 , lineWidth , 0 ,
                                dotsPerAngstrom*a1.x ,
                                dotsPerAngstrom*a1.y ,
                                dotsPerAngstrom*a2.x ,
                                dotsPerAngstrom*a2.y );

            if(bond->getBondType()==Bond::DOUBLE ||
               bond->getBondType()==Bond::AROMATIC ||
               bond->getBondType()==Bond::TRIPLE ){
              RDGeom::Point2D obv=a2-a1;
              RDGeom::Point2D perp=obv;
              perp.rotate90();
              perp.normalize();

              if( (bond->getBondType()==Bond::DOUBLE ||
                   bond->getBondType()==Bond::AROMATIC) &&
                  mol.getRingInfo()->numBondRings(bond->getIdx())){
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

              detail::drawLine( res , atnum1 , atnum2 , lineWidth , (bond->getBondType()==Bond::AROMATIC),
                                dotsPerAngstrom*(offsetStart.x+perp.x) ,
                                dotsPerAngstrom*(offsetStart.y+perp.y) ,
                                dotsPerAngstrom*(offsetStart.x+obv.x+perp.x) ,
                                dotsPerAngstrom*(offsetStart.y+obv.y+perp.y) );
            }
          }
        }
        std::string symbol;
        OrientType orient;
        boost::tie(symbol,orient)=atomSymbols[a1Idx];
        if(symbol!=""){
          res.push_back(ATOM);
          res.push_back(mol[*bAts]->getAtomicNum());
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.x));
          res.push_back(static_cast<ElementType>(dotsPerAngstrom*a1.y));
          res.push_back(static_cast<ElementType>(symbol.length()));
          BOOST_FOREACH(char c, symbol){
            res.push_back(static_cast<ElementType>(c));
          }
          res.push_back(static_cast<ElementType>(orient));
        }

        ++bAts;
      }

      return res;
    }

    std::vector<int> MolToDrawing(const RDKit::ROMol &mol,const std::vector<int> *highlightAtoms=0,
                                  bool kekulize=true){
      RDKit::RWMol *cp = new RDKit::RWMol(mol);
      if(kekulize){
        try{
          RDKit::MolOps::Kekulize(*cp);
        } catch (...) {
          delete cp;
          cp = new RDKit::RWMol(mol);
        }
      }
      if(!mol.getNumConformers()) {
        RDDepict::compute2DCoords(*cp);
      }
      std::vector<int> drawing=DrawMol(*cp,-1,highlightAtoms);
      delete cp;
      return drawing;
    }

  }   // end of namespace Drawing
} // end of namespace RDKit

#endif
