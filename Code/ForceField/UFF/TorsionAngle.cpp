// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TorsionAngle.h"
#include "Params.h"
#include <math.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace UFF {
    namespace Utils {
      double calculateCosTorsion(const RDGeom::Point3D &p1,const RDGeom::Point3D &p2,
                                 const RDGeom::Point3D &p3,const RDGeom::Point3D &p4){
        RDGeom::Point3D r1=p1-p2,r2=p3-p2,r3=p2-p3,r4=p4-p3;
        RDGeom::Point3D t1=r1.crossProduct(r2);
        RDGeom::Point3D t2=r3.crossProduct(r4);
        double d1=t1.length(),d2=t2.length();
        double cosPhi=t1.dotProduct(t2)/(d1*d2);
        return cosPhi;
      }

      // used locally
      bool isInGroup6(int num){
        return (num==8 || num==16 || num==34 || num==52 || num==84);
      }

      // used locally, implement equation 17 of the UFF paper.
      double equation17(double bondOrder23,
                        const AtomicParams *at2Params,
                        const AtomicParams *at3Params){
        return 5.*sqrt(at2Params->U1*at3Params->U1)*(1.+4.18*log(bondOrder23));
      }

    }


    TorsionAngleContrib::TorsionAngleContrib(ForceField *owner,
                                             unsigned int idx1,unsigned int idx2,
                                             unsigned int idx3,unsigned int idx4,
                                             double bondOrder23,
                                             int atNum2,int atNum3,
                                             RDKit::Atom::HybridizationType hyb2,
                                             RDKit::Atom::HybridizationType hyb3,
                                             const AtomicParams *at2Params,
                                             const AtomicParams *at3Params,
                                             bool endAtomIsSP2){
      PRECONDITION(owner,"bad owner");
      PRECONDITION(at2Params,"bad params pointer");
      PRECONDITION(at3Params,"bad params pointer");
      PRECONDITION((idx1!=idx2&&idx1!=idx3&&idx1!=idx4&&idx2!=idx3&&idx2!=idx4&&idx3!=idx4),
                   "degenerate points");
      RANGE_CHECK(0,idx1,owner->positions().size()-1);
      RANGE_CHECK(0,idx2,owner->positions().size()-1);
      RANGE_CHECK(0,idx3,owner->positions().size()-1);
      RANGE_CHECK(0,idx4,owner->positions().size()-1);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_at4Idx = idx4;
    
      this->calcTorsionParams(bondOrder23,atNum2,atNum3,hyb2,hyb3,
                              at2Params,at3Params,endAtomIsSP2);
    }

  
    void TorsionAngleContrib::calcTorsionParams(double bondOrder23,
                                                int atNum2,int atNum3,
                                                RDKit::Atom::HybridizationType hyb2,
                                                RDKit::Atom::HybridizationType hyb3,
                                                const AtomicParams *at2Params,
                                                const AtomicParams *at3Params,
                                                bool endAtomIsSP2){
      PRECONDITION((hyb2==RDKit::Atom::SP2||hyb2==RDKit::Atom::SP3) && (hyb3==RDKit::Atom::SP2||hyb3==RDKit::Atom::SP3),"bad hybridizations");

      if(hyb2==RDKit::Atom::SP3 && hyb3==RDKit::Atom::SP3){
        // general case:
        d_forceConstant = sqrt(at2Params->V1*at3Params->V1);
        d_order = 3;
        d_cosTerm=-1; // phi0=60
     
        // special case for single bonds between group 6 elements:
        if( bondOrder23==1.0 &&
            Utils::isInGroup6(atNum2) && Utils::isInGroup6(atNum3)){
          double V2=6.8,V3=6.8;
          if(atNum2 == 8) V2=2.0;
          if(atNum3 == 8) V3=2.0;
          d_forceConstant = sqrt(V2*V3);
          d_order = 2;
          d_cosTerm=-1; // phi0=90
        }
      } else if(hyb2==RDKit::Atom::SP2 && hyb3==RDKit::Atom::SP2){
        d_forceConstant = Utils::equation17(bondOrder23,at2Params,at3Params);
        d_order = 2;
        // FIX: is this angle term right?
        d_cosTerm = 1.0; // phi0= 180
      } else {
        // SP2 - SP3,  this is, by default, independent of atom type in UFF:
        d_forceConstant = 1.0;
        d_order = 6;
        d_cosTerm = 1.0; // phi0 = 0
        if(bondOrder23==1.0){
          // special case between group 6 sp3 and non-group 6 sp2:
          if( (hyb2==RDKit::Atom::SP3 && Utils::isInGroup6(atNum2) &&
               !Utils::isInGroup6(atNum3)) ||
              (hyb3==RDKit::Atom::SP3 && Utils::isInGroup6(atNum3) &&
               !Utils::isInGroup6(atNum2) ) ){
            d_forceConstant = Utils::equation17(bondOrder23,at2Params,at3Params);
            d_order=2;
            d_cosTerm=-1; // phi0 = 90;
          }

          // special case for sp3 - sp2 - sp2
          // (i.e. the sp2 has another sp2 neighbor, like propene)
          else if(endAtomIsSP2){
            d_forceConstant = 2.0;
            d_order=3;
            d_cosTerm=-1; // phi0 = 180;
          }
        }
      }
    }
    double TorsionAngleContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(this->d_order==2||this->d_order==3||this->d_order==6,"bad order");

      RDGeom::Point3D p1(pos[3*this->d_at1Idx],
                         pos[3*this->d_at1Idx+1],
                         pos[3*this->d_at1Idx+2]);
      RDGeom::Point3D p2(pos[3*this->d_at2Idx],
                         pos[3*this->d_at2Idx+1],
                         pos[3*this->d_at2Idx+2]);
      RDGeom::Point3D p3(pos[3*this->d_at3Idx],
                         pos[3*this->d_at3Idx+1],
                         pos[3*this->d_at3Idx+2]);
      RDGeom::Point3D p4(pos[3*this->d_at4Idx],
                         pos[3*this->d_at4Idx+1],
                         pos[3*this->d_at4Idx+2]);

      double cosPhi=Utils::calculateCosTorsion(p1,p2,p3,p4);
      double sinPhiSq=1-cosPhi*cosPhi;

      // E(phi) = V/2 * (1 - cos(n*phi_0)*cos(n*phi))
      double cosNPhi=0.0;
      switch(this->d_order){
      case 2:
        cosNPhi = cosPhi*cosPhi - sinPhiSq;
        break;
      case 3:
        // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x)
        cosNPhi = cosPhi*(cosPhi*cosPhi - 3*sinPhiSq);
        break;
      case 6:
        // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)
        cosNPhi = 1 + sinPhiSq*(-32*sinPhiSq*sinPhiSq + 48*sinPhiSq - 18);
        break;
      }
      double res=this->d_forceConstant/2.0 * (1 - this->d_cosTerm*cosNPhi);
      //std::cout << " torsion(" << this->d_at1Idx << "," << this->d_at2Idx << "," << this->d_at3Idx << "," << this->d_at4Idx << "): " << cosPhi << "(" << acos(cosPhi) << ")" << " -> " << res << std::endl;
      //if(this->d_at2Idx==5&&this->d_at3Idx==6) std::cerr << " torsion(" << this->d_at1Idx << "," << this->d_at2Idx << "," << this->d_at3Idx << "," << this->d_at4Idx << "): " << cosPhi << "(" << acos(cosPhi) << ")" << " -> " << res << std::endl;
      return res;
    }

    void TorsionAngleContrib::getGrad(double *pos,double *grad) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      RDGeom::Point3D p1(pos[3*this->d_at1Idx],
                         pos[3*this->d_at1Idx+1],
                         pos[3*this->d_at1Idx+2]);
      RDGeom::Point3D p2(pos[3*this->d_at2Idx],
                         pos[3*this->d_at2Idx+1],
                         pos[3*this->d_at2Idx+2]);
      RDGeom::Point3D p3(pos[3*this->d_at3Idx],
                         pos[3*this->d_at3Idx+1],
                         pos[3*this->d_at3Idx+2]);
      RDGeom::Point3D p4(pos[3*this->d_at4Idx],
                         pos[3*this->d_at4Idx+1],
                         pos[3*this->d_at4Idx+2]);
      double *g1=&(grad[3*this->d_at1Idx]);
      double *g2=&(grad[3*this->d_at2Idx]);
      double *g3=&(grad[3*this->d_at3Idx]);
      double *g4=&(grad[3*this->d_at4Idx]);

      RDGeom::Point3D r1=p1-p2,r2=p3-p2,r3=p2-p3,r4=p4-p3;
      RDGeom::Point3D t1=r1.crossProduct(r2);
      RDGeom::Point3D t2=r3.crossProduct(r4);
      double d1=t1.length(),d2=t2.length();
      if(d1==0.0 || d2==0.0){
        return;
      }
      
      double cosPhi=t1.dotProduct(t2)/(d1*d2);
      double sinPhi=1-cosPhi*cosPhi;
      if(sinPhi>=0.0) {
        sinPhi=sqrt(sinPhi);
      } else {
        sinPhi=0.0;
      }
      // dE/dPhi is independent of cartesians:
      double dE_dPhi=getThetaDeriv(cosPhi,sinPhi);
#if 0
      if(dE_dPhi!=dE_dPhi){
        std::cout << "\tNaN in Torsion("<<this->d_at1Idx<<","<<this->d_at2Idx<<","<<this->d_at3Idx<<","<<this->d_at4Idx<<")"<< std::endl;
        std::cout << "sin: " << sinPhi << std::endl;
        std::cout << "cos: " << cosPhi << std::endl;
      } 
      
#endif
      
      // -------
      // dTheta/dx is trickier:
      double dCos_dT1=1./d1 * (t2.x/d2 - cosPhi*t1.x/d1);
      double dCos_dT2=1./d1 * (t2.y/d2 - cosPhi*t1.y/d1);
      double dCos_dT3=1./d1 * (t2.z/d2 - cosPhi*t1.z/d1);
                                                    
      double dCos_dT4=1./d2 * (t1.x/d1 - cosPhi*t2.x/d2);
      double dCos_dT5=1./d2 * (t1.y/d1 - cosPhi*t2.y/d2);
      double dCos_dT6=1./d2 * (t1.z/d1 - cosPhi*t2.z/d2);
    
      double sinTerm;
      // FIX: use a tolerance here:
      if(sinPhi==0.0){
        // this is hacky, but it's per the
        // recommendation from Niketic and Rasmussen:
        sinTerm = 1/cosPhi; 
      } else {
        sinTerm = 1/sinPhi;
      }

      g1[0] += dE_dPhi*sinTerm*(dCos_dT3*r2.y - dCos_dT2*r2.z);
      g1[1] += dE_dPhi*sinTerm*(dCos_dT1*r2.z - dCos_dT3*r2.x);
      g1[2] += dE_dPhi*sinTerm*(dCos_dT2*r2.x - dCos_dT1*r2.y);

      g2[0] += dE_dPhi*sinTerm*(dCos_dT2*(r2.z-r1.z) + dCos_dT3*(r1.y-r2.y) +
                                dCos_dT5*(-1*r4.z) + dCos_dT6*(r4.y));
      g2[1] += dE_dPhi*sinTerm*(dCos_dT1*(r1.z-r2.z) + dCos_dT3*(r2.x-r1.x) +
                                dCos_dT4*(r4.z) + dCos_dT6*(-1*r4.x));
      g2[2] += dE_dPhi*sinTerm*(dCos_dT1*(r2.y-r1.y) + dCos_dT2*(r1.x-r2.x) +
                                dCos_dT4*(-1*r4.y) + dCos_dT5*(r4.x));
    
      g3[0] += dE_dPhi*sinTerm*(dCos_dT2*(r1.z) + dCos_dT3*(-1*r1.y) +
                                dCos_dT5*(r4.z-r3.z) + dCos_dT6*(r3.y-r4.y));
      g3[1] += dE_dPhi*sinTerm*(dCos_dT1*(-1*r1.z) + dCos_dT3*(r1.x) +
                                dCos_dT4*(r3.z-r4.z) + dCos_dT6*(r4.x-r3.x));
      g3[2] += dE_dPhi*sinTerm*(dCos_dT1*(r1.y) + dCos_dT2*(-1*r1.x) +
                                dCos_dT4*(r4.y-r3.y) + dCos_dT5*(r3.x-r4.x));

    
      g4[0] += dE_dPhi*sinTerm*(dCos_dT5*r3.z - dCos_dT6*r3.y);
      g4[1] += dE_dPhi*sinTerm*(dCos_dT6*r3.x - dCos_dT4*r3.z);
      g4[2] += dE_dPhi*sinTerm*(dCos_dT4*r3.y - dCos_dT5*r3.x);

    }


  
    double TorsionAngleContrib::getThetaDeriv(double cosTheta,double sinTheta) const {
      PRECONDITION(this->d_order==2||this->d_order==3||this->d_order==6,"bad order");
      double sinThetaSq=sinTheta*sinTheta;
      // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)

      double res=0.0;
      switch(this->d_order){
      case 2:
        res = 2*sinTheta*cosTheta;
        break;
      case 3:
        // sin(3*x) = 3*sin(x) - 4*sin^3(x)
        res = sinTheta*(3-4*sinThetaSq);
        break;
      case 6:
        // sin(6x) = cos(x) * [ 32*sin^5(x) - 32*sin^3(x) + 6*sin(x) ]
        res = cosTheta*sinTheta * (32*sinThetaSq * (sinThetaSq-1) + 6);
        break;
      }
      res *= this->d_forceConstant/2.0 * this->d_cosTerm * -1 * this->d_order;

      return res;
    }
  }
}
