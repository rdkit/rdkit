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
#include "AngleBend.h"
#include "BondStretch.h"
#include "Params.h"
#include <math.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
  namespace UFF {

    namespace Utils {
      double calcAngleForceConstant(double theta0,
				    double bondOrder12,double bondOrder23,
				    const AtomicParams *at1Params,
				    const AtomicParams *at2Params,
				    const AtomicParams *at3Params){
	double cosTheta0=cos(theta0);
	double r12 = calcBondRestLength(bondOrder12,at1Params,at2Params);
	double r23 = calcBondRestLength(bondOrder23,at2Params,at3Params);
	double r13 = sqrt(r12*r12 + r23*r23 - 2.*r12*r23*cosTheta0);
	double beta = 2.*Params::G/(r12*r23);

	double preFactor = beta*at1Params->Z1*at3Params->Z1 / int_pow<5>(r13);
	double rTerm = r12*r23;
	double innerBit = r12*r23*(1-cosTheta0*cosTheta0) - r13*r13*cosTheta0;
            
	double res=preFactor*rTerm*innerBit;
	return res;
      }
    } // end of namespace Utils
  
    AngleBendContrib::AngleBendContrib(ForceField *owner,
				       unsigned int idx1,unsigned int idx2,unsigned int idx3,
				       double bondOrder12,double bondOrder23,
				       const AtomicParams *at1Params,
				       const AtomicParams *at2Params,
				       const AtomicParams *at3Params,
				       unsigned int order){
      PRECONDITION(owner,"bad owner");
      PRECONDITION(at1Params,"bad params pointer");
      PRECONDITION(at2Params,"bad params pointer");
      PRECONDITION(at3Params,"bad params pointer");
      PRECONDITION((idx1!=idx2&&idx2!=idx3&&idx1!=idx3),"degenerate points");
      RANGE_CHECK(0,idx1,owner->positions().size()-1);
      RANGE_CHECK(0,idx2,owner->positions().size()-1);
      RANGE_CHECK(0,idx3,owner->positions().size()-1);
      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;
      d_at3Idx = idx3;
      d_order = order;
      this->d_forceConstant = Utils::calcAngleForceConstant(at2Params->theta0,
							    bondOrder12,bondOrder23,
							    at1Params,at2Params,at3Params);
      if(order==0){
	double sinTheta0=sin(at2Params->theta0);
	double cosTheta0=cos(at2Params->theta0);
	this->d_C2 = 1./(4.*std::max(sinTheta0*sinTheta0,1e-8));
	this->d_C1 = -4.*this->d_C2*cosTheta0;
	this->d_C0 = this->d_C2*(2.*cosTheta0*cosTheta0 + 1.);
      }
    }

    double AngleBendContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      double dist1=this->dp_forceField->distance(this->d_at1Idx,this->d_at2Idx,pos);
      double dist2=this->dp_forceField->distance(this->d_at2Idx,this->d_at3Idx,pos);

      RDGeom::Point3D p1(pos[3*this->d_at1Idx],
			 pos[3*this->d_at1Idx+1],
			 pos[3*this->d_at1Idx+2]);
      RDGeom::Point3D p2(pos[3*this->d_at2Idx],
			 pos[3*this->d_at2Idx+1],
			 pos[3*this->d_at2Idx+2]);
      RDGeom::Point3D p3(pos[3*this->d_at3Idx],
			 pos[3*this->d_at3Idx+1],
			 pos[3*this->d_at3Idx+2]);
      RDGeom::Point3D p12=p1-p2;
      RDGeom::Point3D p32=p3-p2;
      double cosTheta = p12.dotProduct(p32)/(dist1*dist2);
      // we need sin^2(theta) to get cos(2*theta), so compute that:
      double sinThetaSq = 1-cosTheta*cosTheta;
    
      double angleTerm = this->getEnergyTerm(cosTheta,sinThetaSq);
      double res = this->d_forceConstant*angleTerm;

      return res;
    }

    void AngleBendContrib::getGrad(double *pos,double *grad) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      double dist1=this->dp_forceField->distance(this->d_at1Idx,this->d_at2Idx,pos);
      double dist2=this->dp_forceField->distance(this->d_at2Idx,this->d_at3Idx,pos);

      //std::cout << "\tAngle("<<this->d_at1Idx<<","<<this->d_at2Idx<<","<<this->d_at3Idx<<") " << dist1 << " " << dist2 << std::endl;
      
      RDGeom::Point3D p1(pos[3*this->d_at1Idx],
			 pos[3*this->d_at1Idx+1],
			 pos[3*this->d_at1Idx+2]);
      RDGeom::Point3D p2(pos[3*this->d_at2Idx],
			 pos[3*this->d_at2Idx+1],
			 pos[3*this->d_at2Idx+2]);
      RDGeom::Point3D p3(pos[3*this->d_at3Idx],
			 pos[3*this->d_at3Idx+1],
			 pos[3*this->d_at3Idx+2]);
      double *g1=&(grad[3*this->d_at1Idx]);
      double *g2=&(grad[3*this->d_at2Idx]);
      double *g3=&(grad[3*this->d_at3Idx]);

      RDGeom::Point3D p12=p1-p2;
      RDGeom::Point3D p32=p3-p2;
      double cosTheta = p12.dotProduct(p32)/(dist1*dist2);
      double sinTheta = std::max(sqrt(1.0-cosTheta*cosTheta),1e-8);

      //std::cerr << "GRAD: " << cosTheta << " (" << acos(cosTheta)<< "), ";
      //std::cerr << sinTheta << " (" << asin(sinTheta)<< ")" << std::endl;
    
      // use the chain rule:
      // dE/dx = dE/dTheta * dTheta/dx

      // dE/dTheta is independent of cartesians:
      double dE_dTheta=getThetaDeriv(cosTheta,sinTheta);
    
      // -------
      // dTheta/dx is trickier:
      double dCos_dS1=1./dist1 * (p32.x/dist2 - cosTheta*p12.x/dist1);
      double dCos_dS2=1./dist1 * (p32.y/dist2 - cosTheta*p12.y/dist1);
      double dCos_dS3=1./dist1 * (p32.z/dist2 - cosTheta*p12.z/dist1);

      double dCos_dS4=1./dist2 * (p12.x/dist1 - cosTheta*p32.x/dist2);
      double dCos_dS5=1./dist2 * (p12.y/dist1 - cosTheta*p32.y/dist2);
      double dCos_dS6=1./dist2 * (p12.z/dist1 - cosTheta*p32.z/dist2);

    
      g1[0] += dE_dTheta*dCos_dS1/(-sinTheta);
      g1[1] += dE_dTheta*dCos_dS2/(-sinTheta);
      g1[2] += dE_dTheta*dCos_dS3/(-sinTheta);

      g2[0] += dE_dTheta*(-dCos_dS1 - dCos_dS4)/(-sinTheta);
      g2[1] += dE_dTheta*(-dCos_dS2 - dCos_dS5)/(-sinTheta);
      g2[2] += dE_dTheta*(-dCos_dS3 - dCos_dS6)/(-sinTheta);
    
      g3[0] += dE_dTheta*dCos_dS4/(-sinTheta);
      g3[1] += dE_dTheta*dCos_dS5/(-sinTheta);
      g3[2] += dE_dTheta*dCos_dS6/(-sinTheta);

    }


    double AngleBendContrib::getEnergyTerm(double cosTheta,double sinThetaSq) const {
      PRECONDITION(this->d_order==0||this->d_order==1||this->d_order==2||this->d_order==3||this->d_order==4,"bad order");
      // cos(2x) = cos^2(x) - sin^2(x);
      double cos2Theta = cosTheta*cosTheta - sinThetaSq;

      double res=0.0;
      if(this->d_order==0){
	res=this->d_C0 + this->d_C1*cosTheta + this->d_C2*cos2Theta;
      } else {
	switch(this->d_order){
	case 1:
	  res=cosTheta;
	  break;
	case 2:
	  res=cos2Theta;
	  break;
	case 3:
	  // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x)
	  res = cosTheta*(cosTheta*cosTheta-3*sinThetaSq);
	  break;
	case 4:
	  // cos(4x) = cos^4(x) - 6*cos^2(x)*sin^2(x)+sin^4(x)
	  res = int_pow<4>(cosTheta) - 6*cosTheta*cosTheta*sinThetaSq + sinThetaSq*sinThetaSq;
	  break;
	}
	res = 1-res;
	res /= (this->d_order*this->d_order);
      }
      return res;
    }


  
    double AngleBendContrib::getThetaDeriv(double cosTheta,double sinTheta) const {
      PRECONDITION(this->d_order==0||this->d_order==1||this->d_order==2||this->d_order==3||this->d_order==4,"bad order");

      double dE_dTheta=0.0;
      double sin2Theta = 2*sinTheta*cosTheta;

      if(this->d_order==0){
	dE_dTheta =  -1*this->d_forceConstant*(this->d_C1*sinTheta +
					       2.*this->d_C2*sin2Theta);
      } else {
	// E = k/n^2 [1-cos(n theta)]
	// dE = - k/n^2 * d cos(n theta)

	// these all use:
	// d cos(ax) = -a sin(ax)
      
	switch(this->d_order){
	case 1:
	  dE_dTheta = sinTheta;
	  break;
	case 2:
	  // sin(2*x) = 2*cos(x)*sin(x)
	  dE_dTheta = sin2Theta;
	  break;
	case 3:
	  // sin(3*x) = 3*sin(x) - 4*sin^3(x)
	  dE_dTheta = sinTheta*(3-4*sinTheta*sinTheta);
	  break;
	case 4:
	  // sin(4*x) = cos(x)*(4*sin(x) - 8*sin^3(x))
	  dE_dTheta = cosTheta*sinTheta*(4-8*sinTheta*sinTheta);
	  break;
	}
	dE_dTheta *= this->d_forceConstant/this->d_order;
      }
      return dE_dTheta;
    }
  
  }
}
