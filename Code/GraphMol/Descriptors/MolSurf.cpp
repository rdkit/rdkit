// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <vector>

namespace RDKit{
  namespace Descriptors {
    double getLabuteAtomContribs(const ROMol &mol,
				 std::vector<double> &Vi,
				 double &hContrib,
				 bool includeHs,
				 bool force){
      TEST_ASSERT(Vi.size()==mol.getNumAtoms());
      if(!force && mol.hasProp("_labuteAtomContribs")){
	mol.getProp("_labuteAtomContribs",Vi);
	mol.getProp("_labuteAtomHContrib",hContrib);
	double res;
	mol.getProp("_labuteASA",res);
	return res;
      }
      unsigned int nAtoms=mol.getNumAtoms();
      std::vector<double> rads(nAtoms);
      for(unsigned int i=0;i<nAtoms;++i){
	rads[i]=PeriodicTable::getTable()->getRb0(mol.getAtomWithIdx(i)->getAtomicNum());
	Vi[i]=0.0;
      }

      for(ROMol::ConstBondIterator bondIt=mol.beginBonds();
	  bondIt!=mol.endBonds();++bondIt){
	const double bondScaleFacts[4]={.1,0,.2,.3};
	double Ri=rads[(*bondIt)->getBeginAtomIdx()];
	double Rj=rads[(*bondIt)->getEndAtomIdx()];
	double bij=Ri+Rj;
	if(!(*bondIt)->getIsAromatic()){
	  if((*bondIt)->getBondType()<4){
	    bij -= bondScaleFacts[(*bondIt)->getBondType()];
	  }
	} else {
	  bij -= bondScaleFacts[0];
	}
	double dij=std::min(std::max(fabs(Ri-Rj),bij),Ri+Rj);
	Vi[(*bondIt)->getBeginAtomIdx()] += Rj*Rj-(Ri-dij)*(Ri-dij)/dij;
	Vi[(*bondIt)->getEndAtomIdx()] += Ri*Ri-(Rj-dij)*(Rj-dij)/dij;
      }
      hContrib=0.0;
      if(includeHs){
	double Rj=PeriodicTable::getTable()->getRb0(1);
	for(unsigned int i=0;i<nAtoms;++i){
	  double Ri=rads[i];
	  double bij=Ri+Rj;
	  double dij=std::min(std::max(fabs(Ri-Rj),bij),Ri+Rj);
	  Vi[i] += Rj*Rj-(Ri-dij)*(Ri-dij)/dij;
	  hContrib += Ri*Ri-(Rj-dij)*(Rj-dij)/dij;
	}
      }

      double res=0.0;
      for(unsigned int i=0;i<nAtoms;++i){
	double Ri=rads[i];
	Vi[i] = M_PI*Ri*(4.*Ri-Vi[i]);
	res+=Vi[i];
      }
      if(includeHs){
	double Rj=PeriodicTable::getTable()->getRb0(1);
	hContrib = M_PI*Rj*(4.*Rj-hContrib);
	res+=hContrib;
      }
      mol.setProp("_labuteAtomContribs",Vi,true);
      mol.setProp("_labuteAtomHContrib",hContrib,true);
      mol.setProp("_labuteASA",res,true);

      return res;
    }
    double calcLabuteASA(const ROMol &mol,bool includeHs,bool force){
      if(!force && mol.hasProp("_labuteASA")){
	double res;
	mol.getProp("_labuteASA",res);
	return res;
      }
      std::vector<double> contribs;
      contribs.resize(mol.getNumAtoms());
      double hContrib;
      double res;
      res=getLabuteAtomContribs(mol,contribs,hContrib,includeHs,force);
      return res;
    }

  } // end of namespace Descriptors
} // end of namespace RDKit
