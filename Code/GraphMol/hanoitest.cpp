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
#include <RDGeneral/Invariant.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RankAtoms.h>

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
      TEST_ASSERT(tv[j]>=tv[j-1]);
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
      TEST_ASSERT(tv[j]>=tv[j-1]);
    }
    free(count);
  }
  BOOST_LOG(rdInfoLog)<< "done: " << vects.size()<<std::endl;
}

void test1(){
  BOOST_LOG(rdInfoLog) << "Testing the hanoi sort" << std::endl;

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
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


class atomcomparefunctor {
  Canon::canon_atom *d_atoms;
public:
  atomcomparefunctor() : d_atoms(NULL) {};
  atomcomparefunctor(Canon::canon_atom *atoms) : d_atoms(atoms) {};
  int operator()(int i,int j) const {
    PRECONDITION(d_atoms,"no atoms");
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
  BOOST_LOG(rdInfoLog) << "Testing hanoi with a functor." << std::endl;
  // make sure that hanoi works with a functor and "molecule data"
  {
    std::string smi="FC1C(Cl)C1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    std::vector<int> indices(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
      indices[i]=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    int *data=&indices.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::hanoisort(data,atoms.size(),count,ftor);

    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[indices[i]].invar >= atoms[indices[i-1]].invar);
        if(atoms[indices[i]].invar != atoms[indices[i-1]].invar){
          TEST_ASSERT(count[indices[i]]!=0);
        } else {
          TEST_ASSERT(count[indices[i]]==0);
        }
      } else {
        TEST_ASSERT(count[indices[i]]!=0);
      }
    }
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};

void test3(){
  BOOST_LOG(rdInfoLog) << "Testing basic partition refinement." << std::endl;
  // basic partition refinement
  {
    std::string smi="FC1C(Cl)CCC1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);

    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
  }
  {
    // this time with smarter invariants
    std::string smi="FC1C(Cl)CCC1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<boost::uint64_t> invars(m->getNumAtoms());
    RankAtoms::buildAtomInvariants(*m,invars,true,true);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=invars[i];
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);

    //std::cerr<<"!!!!! ----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


class atomcomparefunctor2 {
  Canon::canon_atom *dp_atoms;
  const ROMol *dp_mol;
  unsigned int getAtomNeighborhood(unsigned int i) const{
    unsigned int res=0;
    const Atom *at=dp_mol->getAtomWithIdx(i);
    std::vector<unsigned int> nbrs(at->getDegree());
    unsigned int nbridx=0;
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = dp_mol->getAtomBonds(at);
    while(beg!=end){
      const BOND_SPTR bond=(*dp_mol)[*beg];
      nbrs[nbridx]=static_cast<unsigned int>(100*bond->getBondTypeAsDouble())+dp_atoms[bond->getOtherAtomIdx(i)].index;
      ++beg;
      ++nbridx;
    }
    std::sort(nbrs.begin(),nbrs.end());
    for(nbridx=0;nbridx<at->getDegree();++nbridx){
      res+=(nbridx*1)*1000+nbrs[nbridx];
    }
    return res;
  }
public:
  atomcomparefunctor2() : dp_atoms(NULL), dp_mol(NULL) {};
  atomcomparefunctor2(Canon::canon_atom *atoms, const ROMol &m) : dp_atoms(atoms), dp_mol(&m) {};
  int operator()(int i,int j) const {
    PRECONDITION(dp_atoms,"no atoms");
    PRECONDITION(dp_mol,"no molecule");
    unsigned int ivi= dp_atoms[i].invar;
    unsigned int ivj= dp_atoms[j].invar;
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;
    else{
      ivi=dp_atoms[i].index+1+getAtomNeighborhood(i);
      ivj=dp_atoms[j].index+1+getAtomNeighborhood(j);
      if(ivi<ivj)
        return -1;
      else if(ivi>ivj)
        return 1;
      else
        return 0;
    }
  }
};


void test4(){
  BOOST_LOG(rdInfoLog) << "Testing partition refinement with neighbors." << std::endl;
  // partition refinement with neighbors
  {
    std::string smi="FC1C(Cl)CCC1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }


    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    atomcomparefunctor2 ftor2(&atoms.front(),*m);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor2,true,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
  }

  {
    std::string smi="FC1C(CO)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }

    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    atomcomparefunctor2 ftor2(&atoms.front(),*m);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor2,true,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
      }
      TEST_ASSERT(count[order[i]]==1); // here we should manage to uniquify everything
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
  }

  {
    std::string smi="FC1C(CC)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }


    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    atomcomparefunctor2 ftor2(&atoms.front(),*m);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor2,true,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0]==4 && count[4]==2); 
    TEST_ASSERT(order[1]==9 && count[9]==0); 
    TEST_ASSERT(order[2]==3 && count[3]==2); 
    TEST_ASSERT(order[3]==8 && count[8]==0); 
    TEST_ASSERT(order[4]==5 && count[5]==2); 
    TEST_ASSERT(order[5]==6 && count[6]==0); 
    TEST_ASSERT(order[6]==2 && count[2]==2); 
    TEST_ASSERT(order[7]==7 && count[7]==0); 
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


void test5(){
  BOOST_LOG(rdInfoLog) << "testing canonicalization via tie breaking." << std::endl;
  // canonicalization via tie breaking
  {
    std::string smi="FC1C(CC)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].invar=m->getAtomWithIdx(i)->getAtomicNum();
      atoms[i].index=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<std::endl;
    // }

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
        if(atoms[order[i]].invar != atoms[order[i-1]].invar){
          TEST_ASSERT(count[order[i]]!=0);
        } else {
          TEST_ASSERT(count[order[i]]==0);
        }
      } else {
        TEST_ASSERT(count[order[i]]!=0);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }


    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    atomcomparefunctor2 ftor2(&atoms.front(),*m);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor2,true,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;
    
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
      }
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    }
    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0]==4 && count[4]==2); 
    TEST_ASSERT(order[1]==9 && count[9]==0); 
    TEST_ASSERT(order[2]==3 && count[3]==2); 
    TEST_ASSERT(order[3]==8 && count[8]==0); 
    TEST_ASSERT(order[4]==5 && count[5]==2); 
    TEST_ASSERT(order[5]==6 && count[6]==0); 
    TEST_ASSERT(order[6]==2 && count[2]==2); 
    TEST_ASSERT(order[7]==7 && count[7]==0); 

    RDKit::Canon::BreakTies(*m,data,ftor2,true,order,count,activeset,next);
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      if(i>0){
        TEST_ASSERT(atoms[order[i]].invar >= atoms[order[i-1]].invar);
      }
      TEST_ASSERT(count[order[i]]==1);
    }
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


void test6(){
  BOOST_LOG(rdInfoLog) << "testing canonicalization using the wrapper." << std::endl;
  // canonicalization using the wrapper
  {
    std::string smi="FC1C(CC)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::RankMolAtoms(*m,atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      TEST_ASSERT(!seen[i]);
      seen.set(i,1);
    }
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};



int main(){
  RDLog::InitLogs();
  //test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  return 0;
}

