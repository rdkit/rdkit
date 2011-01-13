// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolOps.h"
#include "RDKitBase.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include <iomanip>


namespace RDKit {
  namespace MolOps {
    double computeBalabanJ(double *distMat, int nb, int nAts){
      // NOTE that the distance matrix is modified here for the sake of 
      // efficiency
      PRECONDITION(distMat,"bogus distance matrix")
        double sum = 0.0;
      int nActive=nAts;
      int mu = nb - nActive + 1;
    
      if(mu==-1) return 0.0;
    
      for(int i=0;i<nAts;i++){
        int iTab=i*nAts;
        sum = 0.0;
        for(int j=0;j<nAts;j++){
          if(j!=i){
            sum += distMat[iTab+j];
          }
        }
        distMat[iTab+i]*=sum;
      }
      double accum = 0.0;
      for (int i = 0; i < nAts; i++) {
        int iTab = i*nAts+i;
        for (int j = i+1; j < nAts ; j++) {
          // NOTE: this isn't strictly the Balaban J value, because we
          // aren't only adding in adjacent atoms.  Since we're doing a
          // discriminator, that shouldn't be a problem.
          if(j!=i){
            accum += (1.0/sqrt(distMat[iTab]*distMat[j*nAts+j]));
          }
        }
      }
      return nActive/((mu+1)*accum);
    }

    double computeBalabanJ(const ROMol &mol, 
                           bool useBO,
                           bool force,
                           const std::vector<int> *bondPath,
                           bool cacheIt) {
      double res=0.0;
      if (!force && mol.hasProp("BalanbanJ")) {
        mol.getProp("BalabanJ", res);
      }
      else {
        double *dMat;
        int nb=0,nAts=0;
        if(bondPath){
          boost::dynamic_bitset<> atomsUsed(mol.getNumAtoms());
          boost::dynamic_bitset<> bondsUsed(mol.getNumBonds());
          for(std::vector<int>::const_iterator ci=bondPath->begin();
              ci!=bondPath->end();ci++){
            bondsUsed[*ci]=1;
          }
          std::vector<const Bond *> bonds;
          bonds.reserve(bondPath->size());
          std::vector<int> atomsInPath;
          atomsInPath.reserve(bondPath->size()+1);

          ROMol::EDGE_ITER beg,end;
          boost::tie(beg,end)=mol.getEdges();
          while(beg!=end){
            const Bond *bond=mol[*beg].get();
            if(bondsUsed[bond->getIdx()]){
              int begIdx=bond->getBeginAtomIdx();
              int endIdx=bond->getEndAtomIdx();
              bonds.push_back(bond);
              if(!atomsUsed[begIdx]){
                atomsInPath.push_back(begIdx);
                atomsUsed[begIdx]=1;
              }
              if(!atomsUsed[endIdx]){
                atomsInPath.push_back(endIdx);
                atomsUsed[endIdx]=1;
              }
            }
            beg++;
          }
          nb = bondPath->size();
          nAts = atomsInPath.size();
          dMat = MolOps::getDistanceMat(mol,atomsInPath,bonds,true,true);
          res = computeBalabanJ(dMat,nb,nAts);
          delete [] dMat;
        } else {
          nb = mol.getNumBonds();
          nAts = mol.getNumAtoms();
          dMat = MolOps::getDistanceMat(mol,true,true,true,0);
          res = computeBalabanJ(dMat,nb,nAts);
          delete [] dMat;
        }


        if(cacheIt) mol.setProp("BalabanJ", res, true);
      }
      return res;
    }
  } // end of namespace MolOps
} // end of namespace RDKit
