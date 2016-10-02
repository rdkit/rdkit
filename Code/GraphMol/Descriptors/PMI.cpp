//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include "PMI.h"
#include <Eigen/Dense>

namespace RDKit {
namespace Descriptors{
namespace {

bool getMoments(const ROMol& mol,int confId, bool useAtomicMasses,
                double &pm1, double &pm2, double &pm3){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  const char *pn1 = useAtomicMasses ? "PMI1_mass" : "PMI1";
  const char *pn2 = useAtomicMasses ? "PMI2_mass" : "PMI2";
  const char *pn3 = useAtomicMasses ? "PMI3_mass" : "PMI3";

  if(mol.hasProp(pn1) && mol.hasProp(pn2) && mol.hasProp(pn3) ) {
    mol.getProp(pn1,pm1);
    mol.getProp(pn2,pm2);
    mol.getProp(pn3,pm3);
    return true;
  }

  const Conformer &conf=mol.getConformer(confId);

  Eigen::Matrix3d axes;
  Eigen::Vector3d moments;
  bool res;
  bool ignoreHs=false;
  if(useAtomicMasses){
    std::vector<double> weights;
    weights.resize(mol.getNumAtoms());
    for(ROMol::ConstAtomIterator cai=mol.beginAtoms();
        cai!=mol.endAtoms();++cai){
      weights[(*cai)->getIdx()] = (*cai)->getMass();
    }
    res = MolTransforms::computePrincipalAxesAndMoments(conf,axes,moments,ignoreHs,false,&weights);
  } else {
    res = MolTransforms::computePrincipalAxesAndMoments(conf,axes,moments,ignoreHs);
  }
  if(res){
    pm1 = moments(0);
    pm2 = moments(1);
    pm3 = moments(2);
    // std::cerr<<"   moments: "<<pm1<<" "<<pm2<<" "<<pm3<<std::endl;
    mol.setProp(pn1,pm1,true);
    mol.setProp(pn2,pm2,true);
    mol.setProp(pn3,pm3,true);
  }
  return res;
}

} // end of anonymous namespace

double NPR1(const ROMol& mol, int confId, bool useAtomicMasses){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  double pm1,pm2,pm3;
  if(!getMoments(mol,confId,useAtomicMasses,pm1,pm2,pm3)){
    // the eigenvector calculation failed
    return 0.0; // FIX: throw an exception here?
  }
  if(pm3<1e-8) return 0.0;
  return pm1/pm3;
}
double NPR2(const ROMol& mol, int confId, bool useAtomicMasses){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  double pm1,pm2,pm3;
  if(!getMoments(mol,confId,useAtomicMasses,pm1,pm2,pm3)){
    // the eigenvector calculation failed
    return 0.0; // FIX: throw an exception here?
  }
  if(pm3<1e-8) return 0.0;
  return pm2/pm3;
}
double PMI1(const ROMol& mol, int confId, bool useAtomicMasses){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  double pm1,pm2,pm3;
  if(!getMoments(mol,confId,useAtomicMasses,pm1,pm2,pm3)){
    // the eigenvector calculation failed
    return 0.0; // FIX: throw an exception here?
  }
  return pm1;
}
double PMI2(const ROMol& mol, int confId, bool useAtomicMasses){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  double pm1,pm2,pm3;
  if(!getMoments(mol,confId,useAtomicMasses,pm1,pm2,pm3)){
    // the eigenvector calculation failed
    return 0.0; // FIX: throw an exception here?
  }
  return pm2;
}
double PMI3(const ROMol& mol, int confId, bool useAtomicMasses){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers");
  double pm1,pm2,pm3;
  if(!getMoments(mol,confId,useAtomicMasses,pm1,pm2,pm3)){
    // the eigenvector calculation failed
    return 0.0; // FIX: throw an exception here?
  }
  return pm3;
}

} // end of Descriptors namespace
} // end of RDKit namespace
