// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/MolOps.h>
#include <memory>
#include <boost/shared_array.hpp>
#include <algorithm>
#include <cstring>

namespace RDKit{

  const int LOCAL_INF=(int)1e8;

  // local utility namespace 
  namespace {
    /* ----------------------------------------------

    This implements the Floyd-Warshall all-pairs-shortest-paths
    algorithm as described on pg 564 of the White Book (Cormen, Leiserson, Rivest)

    Arguments:
    adjMat: the adjacency matrix.  This is overwritten with the shortest
    paths, should be dim x dim
    dim: the size of adjMat
    pathMat: the path matrix, should be dim x dim
   
    -----------------------------------------------*/
    template<class T> void
    FloydWarshall(int dim,T *adjMat,int *pathMat)
    {
      int k,i,j;
      T *currD,*lastD,*tTemp;
      int *currP,*lastP,*iTemp;

      currD = new T[dim*dim];
      currP = new int[dim*dim];
      lastD = new T[dim*dim];
      lastP = new int[dim*dim];

      memcpy(static_cast<void *>(lastD),static_cast<void *>(adjMat),dim*dim*sizeof(T));
  
      // initialize the paths
      for(i=0;i<dim;i++){
	int itab=i*dim;
	for(j=0;j<dim;j++){
	  if(i==j || adjMat[itab+j]==LOCAL_INF){
	    pathMat[itab+j] = -1;
	  }else{
	    pathMat[itab+j] = i;
	  }
	}
      }
      memcpy(static_cast<void *>(lastP),static_cast<void *>(pathMat),dim*dim*sizeof(int));
  
      for(k=0;k<dim;k++){
	int ktab=k*dim;
	for(i=0;i<dim;i++){
	  int itab=i*dim;
	  for(j=0;j<dim;j++){
	    T v1=lastD[itab+j];
	    T v2=lastD[itab+k]+lastD[ktab+j];
	    if(v1<=v2){
	      currD[itab+j] = v1;
	      currP[itab+j] = lastP[itab+j];
	    }else{
	      currD[itab+j] = v2;
	      currP[itab+j] = lastP[ktab+j];
	    }
	  }
	}
	tTemp = currD;
	currD = lastD;
	lastD = tTemp;

	iTemp = currP;
	currP = lastP;
	lastP = iTemp;
      }
      memcpy(static_cast<void *>(adjMat),static_cast<void *>(lastD),dim*dim*sizeof(T));
      memcpy(static_cast<void *>(pathMat),static_cast<void *>(lastP),dim*dim*sizeof(int));

      delete [] currD;
      delete [] currP;
      delete [] lastD;
      delete [] lastP;
    }


    template<class T> void
    FloydWarshall(int dim,T *adjMat,int *pathMat,const std::vector<int> &activeAtoms)
    {
      T *currD,*lastD,*tTemp;
      int *currP,*lastP,*iTemp;

      currD = new T[dim*dim];
      currP = new int[dim*dim];
      lastD = new T[dim*dim];
      lastP = new int[dim*dim];

      memcpy(static_cast<void *>(lastD),static_cast<void *>(adjMat),dim*dim*sizeof(T));
  
      // initialize the paths
      for(std::vector<int>::const_iterator i=activeAtoms.begin();
	  i!=activeAtoms.end();
	  i++){
	int itab=(*i)*dim;
	for(std::vector<int>::const_iterator j=activeAtoms.begin();
	    j!=activeAtoms.end();
	    j++){
	  if((*i)==(*j) || adjMat[itab+(*j)]==LOCAL_INF){
	    pathMat[itab+(*j)] = -1;
	  }else{
	    pathMat[itab+(*j)] = (*i);
	  }
	}
      }
      memcpy(static_cast<void *>(lastP),static_cast<void *>(pathMat),dim*dim*sizeof(int));
  
      for(std::vector<int>::const_iterator k=activeAtoms.begin();
	  k!=activeAtoms.end();
	  k++){
	int ktab=(*k)*dim;
	for(std::vector<int>::const_iterator i=activeAtoms.begin();
	    i!=activeAtoms.end();
	    i++){
	  int itab=(*i)*dim;
	  for(std::vector<int>::const_iterator j=activeAtoms.begin();
	      j!=activeAtoms.end();
	      j++){
	    T v1=lastD[itab+(*j)];
	    T v2=lastD[itab+(*k)]+lastD[ktab+(*j)];
	    if(v1<=v2){
	      currD[itab+(*j)] = v1;
	      currP[itab+(*j)] = lastP[itab+(*j)];
	    }else{
	      currD[itab+(*j)] = v2;
	      currP[itab+(*j)] = lastP[ktab+(*j)];
	    }
	  }
	}
	tTemp = currD;
	currD = lastD;
	lastD = tTemp;

	iTemp = currP;
	currP = lastP;
	lastP = iTemp;
      }
      memcpy(static_cast<void *>(adjMat),static_cast<void *>(lastD),dim*dim*sizeof(T));
      memcpy(static_cast<void *>(pathMat),static_cast<void *>(lastP),dim*dim*sizeof(int));

      delete [] currD;
      delete [] currP;
      delete [] lastD;
      delete [] lastP;
    }
  } // end of local utility namespace
  
  namespace MolOps {
    double *getDistanceMat(const ROMol &mol,bool useBO,
				   bool useAtomWts,
				   bool force,
				   const char *propNamePrefix){
      std::string propName;
      boost::shared_array<double> sptr;
      if(propNamePrefix){
	propName = propNamePrefix;
      } else {
	propName = "";
      }
      propName+="DistanceMatrix";
      // make sure we don't use the nonBO cache for the BO matrix and vice versa:
      if(useBO) propName+="BO";
      if(!force && mol.hasProp(propName)){
	mol.getProp(propName,sptr);
	return sptr.get();
      }    
      int nAts=mol.getNumAtoms();
      double *dMat = new double[nAts*nAts];
      int i,j;
      // initialize off diagonals to LOCAL_INF and diagonals to 0
      for(i=0;i<nAts*nAts;i++) dMat[i] = LOCAL_INF;
      for(i=0;i<nAts;i++) dMat[i*nAts+i]=0.0;

      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
	const BOND_SPTR bond = mol[*firstB];
	i = bond->getBeginAtomIdx();
	j = bond->getEndAtomIdx();
	double contrib;
	if(useBO){
	  if(!bond->getIsAromatic()){
	    contrib = 1./bond->getBondTypeAsDouble();
	  } else {
	    contrib = 2. / 3.;
	  }
	} else {
	  contrib = 1.0;
	}
	dMat[i*nAts+j]=contrib;
	dMat[j*nAts+i]=contrib;
	++firstB;
      }

      int *pathMat = new int[nAts*nAts];
      memset(static_cast<void *>(pathMat),0,nAts*nAts*sizeof(int));
      FloydWarshall(nAts,dMat,pathMat);
    
      if(useAtomWts){
	for (i = 0; i < nAts; i++) {
	  int anum = mol.getAtomWithIdx(i)->getAtomicNum();
	  dMat[i*nAts+i] = 6.0/anum;
	}
      }
      sptr.reset(dMat);
      mol.setProp(propName,sptr,true);
      boost::shared_array<int> iSptr(pathMat);
      mol.setProp(propName+"_Paths",iSptr,true);
  
      return dMat;
    };



    double *getDistanceMat(const ROMol &mol,
				   const std::vector<int> &activeAtoms,
				   const std::vector<const Bond *> &bonds,
				   bool useBO,
				   bool useAtomWts){
  
      const int nAts=activeAtoms.size();
  
      double *dMat = new double[nAts*nAts];
      int i,j;
      // initialize off diagonals to LOCAL_INF and diagonals to 0
      for(i=0;i<nAts*nAts;i++) dMat[i] = LOCAL_INF;
      for(i=0;i<nAts;i++) dMat[i*nAts+i]=0.0;

      for(std::vector<const Bond *>::const_iterator bi=bonds.begin();
	  bi!=bonds.end();bi++){
	const Bond *bond=*bi;
	i=std::find(activeAtoms.begin(),activeAtoms.end(),
		    static_cast<int>(bond->getBeginAtomIdx())) - activeAtoms.begin();
	j=std::find(activeAtoms.begin(),activeAtoms.end(),
		    static_cast<int>(bond->getEndAtomIdx())) - activeAtoms.begin();
	double contrib;
	if(useBO){
	  if(!bond->getIsAromatic()){
	    contrib = 1./bond->getBondTypeAsDouble();
	  } else {
	    contrib = 2. / 3.;
	  }
	} else {
	  contrib = 1.0;
	}
	dMat[i*nAts+j]=contrib;
	dMat[j*nAts+i]=contrib;
      }

      int *pathMat = new int[nAts*nAts];
      memset(static_cast<void *>(pathMat),0,nAts*nAts*sizeof(int));
      FloydWarshall(nAts,dMat,pathMat);
      delete [] pathMat;

      if(useAtomWts){
	for(i=0;i<nAts;i++){
	  int anum = mol.getAtomWithIdx(activeAtoms[i])->getAtomicNum();
	  dMat[i*nAts+i] = 6.0/anum;
	}
      }
      return dMat;
    };


    // NOTE: do *not* delete results
    double *getAdjacencyMatrix(const ROMol &mol,
                               bool useBO,
                               int emptyVal,
                               bool force,
                               const char *propNamePrefix,
                               const boost::dynamic_bitset<> *bondsToUse
                               )
    {
      std::string propName;
      boost::shared_array<double> sptr;
      if(propNamePrefix){
	propName = propNamePrefix;
      } else {
	propName = "";
      }
      propName+="AdjacencyMatrix";
      if(!force && mol.hasProp(propName)){
	mol.getProp(propName,sptr);
	return sptr.get();
      }    

      int nAts = mol.getNumAtoms();
      double *res = new double[nAts*nAts];
      memset(static_cast<void *>(res),emptyVal,nAts*nAts*sizeof(double));

      for(ROMol::ConstBondIterator bondIt=mol.beginBonds();
          bondIt!=mol.endBonds();bondIt++){
        if(bondsToUse && !(*bondsToUse)[(*bondIt)->getIdx()]){
          continue;
        }
        if(!useBO){
          int beg=(*bondIt)->getBeginAtomIdx();
          int end=(*bondIt)->getEndAtomIdx();
          res[beg*nAts+end] = 1;
          res[end*nAts+beg] = 1;
        }
        else {
          int begIdx=(*bondIt)->getBeginAtomIdx();
          int endIdx=(*bondIt)->getEndAtomIdx();
          Atom const *beg=mol.getAtomWithIdx(begIdx);
          Atom const *end=mol.getAtomWithIdx(endIdx);
          res[begIdx*nAts+endIdx] = (*bondIt)->getValenceContrib(beg);
          res[endIdx*nAts+begIdx] = (*bondIt)->getValenceContrib(end);
        }
      }
      sptr.reset(res);
      mol.setProp(propName,sptr,true);
  
      return res;
    };
  
    INT_LIST getShortestPath(const ROMol &mol, int aid1, int aid2) {

      int nats = mol.getNumAtoms();
      RANGE_CHECK(0,aid1,nats-1);
      RANGE_CHECK(0,aid2,nats-1);
      CHECK_INVARIANT(aid1 != aid2, "");

      INT_VECT pred, doneAtms;
  
      //pred.reserve(nats);
      //doneAtms.reserve(nats);
      pred.resize(nats);
      doneAtms.resize(nats);
      int ai;
      for (ai = 0; ai < nats; ++ai) {
	doneAtms[ai] = 0;
      }

      std::deque<int> bfsQ;
  
      bfsQ.push_back(aid1);
      bool done = false;
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      while ((!done) && (bfsQ.size() > 0)) {
	int curAid = bfsQ.front();
	boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(curAid));
	while (nbrIdx != endNbrs) {
	  if (doneAtms[*nbrIdx] == 0) {
	    pred[*nbrIdx] = curAid;
	    if (static_cast<int>(*nbrIdx) == aid2) {
	      done = true;
	      break;
	    }
	    bfsQ.push_back(*nbrIdx);
	  }
	  nbrIdx++;
	}
	doneAtms[curAid] = 1;
	bfsQ.pop_front();
      }

      INT_LIST res;
      if(done){
        done = false;
        int prev = aid2;
        res.push_back(aid2);
        while (!done) {
          prev = pred[prev];
          if (prev != aid1) {
            res.push_front(prev);
          } else {
            done = true;
          }
        }
        res.push_front(aid1);
      }  
      return res;
 
    }

    double *get3DDistanceMat(const ROMol &mol,
                             int confId,
                             bool useAtomWts,
                             bool force,
                             const char *propNamePrefix){
      const Conformer &conf=mol.getConformer(confId);
      std::string propName;
      boost::shared_array<double> sptr;
      if(propNamePrefix){
	propName = propNamePrefix;
      } else {
	propName = "_";
      }
      propName+="3DDistanceMatrix_Conf"+boost::lexical_cast<std::string>(conf.getId());
      if(!force && mol.hasProp(propName)){
	mol.getProp(propName,sptr);
	return sptr.get();
      }    

      unsigned int nAts=mol.getNumAtoms();
      double *dMat = new double[nAts*nAts];
      sptr.reset(dMat);

      for(unsigned int i=0;i<nAts;++i){
        if(useAtomWts){
          dMat[i*nAts+i]=6.0/mol.getAtomWithIdx(i)->getAtomicNum();
        } else {
          dMat[i*nAts+i]=0.0;
        }
        for(unsigned int j=i+1;j<nAts;++j){
          double dist = (conf.getAtomPos(i)-conf.getAtomPos(j)).length();
          dMat[i*nAts+j]=dist;
          dMat[j*nAts+i]=dist;          
        }
      }
      
      mol.setProp(propName,sptr,true);
      return dMat;
    }
  } // end of namespace MolOps
} // end of namespace RDKit
