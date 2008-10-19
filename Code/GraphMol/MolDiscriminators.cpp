// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "MolOps.h"
#include "RDKitBase.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include <iomanip>

#ifdef RDK_USELAPACKPP
//lapack ++ includes
#include <lafnames.h>
#include <lapack.h>
#include <symd.h>
#include <lavd.h>
#include <laslv.h>
#else
// uBLAS and boost.bindings includes
#include <boost/numeric/bindings/traits/ublas_matrix.hpp> 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/ublas/io.hpp> 
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
#endif

namespace RDKit {

  // local utility namespace
  namespace {
    double LocalBalaban(double *distMat, int nb, int nAts){
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
  }
  
  namespace MolOps {
    DiscrimTuple computeDiscriminators(double *distMat, unsigned int nb, unsigned int na) {
      PRECONDITION(distMat,"bogus distance matrix");
      double ev1 = 0.0;
      double ev2 = 0.0;
#ifdef RDK_USELAPACKPP
      LaSymmMatDouble A(distMat, na, na);
      LaVectorDouble eigs(na);
      LaEigSolve(A, eigs);
#else
      ublas::matrix<double> A(na,na);
      ublas::vector<double> eigs(na);
      for(unsigned int i=0;i<na;++i){
        for(unsigned int j=i;j<na;++j){
          A(i,j)=distMat[i*na+j];
        }
      }
      lapack::syev('N','L',A,eigs);
#endif
      if (na > 1) {
        ev1 = eigs(0);
      }
      if (na > 2) {
        ev2 = eigs(na-1); 
      }
      double J = LocalBalaban(distMat, nb, na);
    
      return boost::make_tuple(J,ev1,ev2);
    }

    DiscrimTuple computeDiscriminators(const ROMol &mol, 
                                       bool useBO,
                                       bool force) {
      DiscrimTuple res;
      if ((mol.hasProp("Discrims")) && (!force)) {
        mol.getProp("Discrims", res);
      }
      else {
        unsigned int nAts = mol.getNumAtoms();
        double *dMat;
        unsigned int nb = mol.getNumBonds();
        dMat = MolOps::getDistanceMat(mol,useBO,true);
        //  Our discriminators (based upon eigenvalues of the distance matrix
        //  and Balaban indices) are good, but is a case they don't properly
        //  handle by default.  These two fragments:
        //    C-ccc and c-ccc
        //  give the same discriminators because their distance matrices are
        //  identical.  We'll work around this by adding 0.5 to the diagonal
        //  elements of the distance matrix corresponding to aromatic atoms:
        ROMol::ConstAromaticAtomIterator atomIt;
        for(atomIt=mol.beginAromaticAtoms();
            atomIt!=mol.endAromaticAtoms();
            atomIt++){
          unsigned int idx=(*atomIt)->getIdx();
          dMat[idx*nAts+idx] += 0.5;
        }
#if 0
        BOOST_LOG(rdDebugLog)<< "--------------------" << std::endl;
        for(int i=0;i<nAts;i++){
          for(int j=0;j<nAts;j++){
            BOOST_LOG(rdDebugLog)<< "\t" << std::setprecision(4) << dMat[i*nAts+j];
          }
          BOOST_LOG(rdDebugLog)<< std::endl;
        }
#endif
      
        res = MolOps::computeDiscriminators(dMat, nb, nAts);
        mol.setProp("Discrims", res, true);
      }
      return res;
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

          ROMol::GRAPH_MOL_BOND_PMAP::const_type bondMap = mol.getBondPMap();
          ROMol::EDGE_ITER beg,end;
          boost::tie(beg,end)=mol.getEdges();
          while(beg!=end){
            const Bond *bond=bondMap[*beg];
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
          res = LocalBalaban(dMat,nb,nAts);
          delete [] dMat;
        } else {
          nb = mol.getNumBonds();
          nAts = mol.getNumAtoms();
          dMat = MolOps::getDistanceMat(mol,true,true,true,0);
          res = LocalBalaban(dMat,nb,nAts);
          delete [] dMat;
        }


        if(cacheIt) mol.setProp("BalabanJ", res, true);
      }
      return res;
    }
  } // end of namespace MolOps
} // end of namespace RDKit
