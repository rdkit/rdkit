//
//  Copyright (C) 2014 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/roger_canon.h>
#include <RDGeneral/RDLog.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>
#include <vector>
#include <boost/random.hpp>
#include <cstdlib>

using namespace RDKit;

int pcmp(const void *a,const void *b){
  if((*(int *)a)<(*(int *)b)){
    return -1;
  } else if((*(int *)a)>(*(int *)b)){
    return 1;
  }
  return 0;
}
int icmp(int a,int b){
  if(a<b){
    return -1;
  } else if(a>b){
    return 1;
  }
  return 0;
}

void qs1(  const std::vector< std::vector<int> > &vects){
  BOOST_LOG(rdInfoLog)<<"sorting (qsort) vectors"<<std::endl;
  for(unsigned int i=0;i<vects.size();++i){
    std::vector<int> tv=vects[i];
    int *data=&tv.front();
    qsort(data,tv.size(),sizeof(int),pcmp);
    for(unsigned int j=1;j<tv.size();++j){
      assert(tv[j]>=tv[j-1]);
    }
  }
  BOOST_LOG(rdInfoLog)<< "done: " << vects.size()<<std::endl;
}

void hs1(  const std::vector< std::vector<int> > &vects){
  BOOST_LOG(rdInfoLog)<<"sorting (hanoi sort) vectors"<<std::endl;
  for(unsigned int i=0;i<vects.size();++i){
    std::vector<int> tv=vects[i];
    int *data=&tv.front();
    int *count=(int *)malloc(tv.size()*sizeof(int));
    RDKit::Canon::hanoisort(data,tv.size(),count,icmp);
    for(unsigned int j=1;j<tv.size();++j){
      assert(tv[j]>=tv[j-1]);
    }
    free(count);
  }
  BOOST_LOG(rdInfoLog)<< "done: " << vects.size()<<std::endl;
}

void test1(){
  typedef boost::random::mersenne_twister<boost::uint32_t,32,4,2,31,0x9908b0df,11,7,0x9d2c5680,15,0xefc60000,18, 3346425566U>  rng_type;
  typedef boost::uniform_int<> distrib_type;
  typedef boost::variate_generator<rng_type &,distrib_type> source_type;
  rng_type generator(42u);

  const unsigned int nVects=500000;
  const unsigned int vectSize=50;
  const unsigned int nClasses=15;

  distrib_type dist(0,nClasses);
  source_type randomSource(generator,dist);

  BOOST_LOG(rdInfoLog)<<"populating vectors"<<std::endl;
  std::vector< std::vector<int> > vects(nVects);
  for(unsigned int i=0;i<nVects;++i){
    vects[i] = std::vector<int>(vectSize);
    for(unsigned int j=0;j<nClasses;++j){
      vects[i][j] = j;
    }
    for(unsigned int j=nClasses+1;j<vectSize;++j){
      vects[i][j] = randomSource();
    }
  }

  qs1(vects);
  hs1(vects);
};


class atomcomparefunctor {
  Canon::canon_atom *d_atoms;
  atomcomparefunctor() : d_atoms(NULL) {};
  atomcomparefunctor(Canon::canon_atom *atoms) : d_atoms(atoms) {};
  int operator()(int i,int j) const {
    assert(d_atoms);
    unsigned int ivi= d_atoms[i].invar;
    unsigned int ivj= d_atoms[j].invar;
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;
    else
      return 0;
  }
};
void test2(){
  {
    std::string smi="FC1C(Cl)C1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
    }
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      std::cerr<<i<<": "<<atoms[i].invar<<" "<<atoms[i].index<<std::endl;
    }

  }
};


int main(){
  RDLog::InitLogs();
  test1();
  test2();
  return 0;
}

