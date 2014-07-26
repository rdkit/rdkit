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

class int_compare_ftor {
  const int *dp_ints;
public:
  int_compare_ftor() : dp_ints(NULL) {};
  int_compare_ftor(const int *ints) : dp_ints(ints) {};
  int operator()(int i,int j) const {
    PRECONDITION(dp_ints,"no ints");
    unsigned int ivi= dp_ints[i];
    unsigned int ivj= dp_ints[j];
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;
    else
      return 0;
  }
};

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
    const int *data=&vects[i].front();
    int_compare_ftor icmp(data);
    int *indices=(int *)malloc(vects[i].size()*sizeof(int));
    for(unsigned int j=0;j<vects[i].size();++j) indices[j]=j;
    int *count=(int *)malloc(vects[i].size()*sizeof(int));
    RDKit::Canon::hanoisort(indices,vects[i].size(),count,icmp);
    for(unsigned int j=1;j<vects[i].size();++j){
      TEST_ASSERT(data[indices[j]]>=data[indices[j-1]]);
    }
    free(count);
    free(indices);
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
    for(unsigned int j=0;j<vectSize;++j){
      vects[i][j] = randomSource();
    }
  }

  //qs1(vects);
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
    unsigned int ivi,ivj;

    // always start with the current class:
    ivi= d_atoms[i].index;
    ivj= d_atoms[j].index;
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    ivi= d_atoms[i].atom->getAtomicNum();
    ivj= d_atoms[j].atom->getAtomicNum();
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    return 0;
  }
};
class atomcomparefunctor2 {
  Canon::canon_atom *d_atoms;
public:
  atomcomparefunctor2() : d_atoms(NULL) {};
  atomcomparefunctor2(Canon::canon_atom *atoms) : d_atoms(atoms) {};
  int operator()(int i,int j) const {
    PRECONDITION(d_atoms,"no atoms");
    unsigned int ivi,ivj;

    // always start with the current class:
    ivi= d_atoms[i].index;
    ivj= d_atoms[j].index;
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    // start by comparing degree
    ivi= d_atoms[i].atom->getDegree();
    ivj= d_atoms[j].atom->getDegree();
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    // move onto atomic number
    ivi= d_atoms[i].atom->getAtomicNum();
    ivj= d_atoms[j].atom->getAtomicNum();
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

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
      atoms[i].atom = m->getAtomWithIdx(i);
      atoms[i].index=0;
      indices[i]=i;
    }
    atomcomparefunctor ftor(&atoms.front());

    int *data=&indices.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::hanoisort(data,atoms.size(),count,ftor);

    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<indices[i]<<" "<<" index: "<<atoms[indices[i]].index<<" count: "<<count[indices[i]]<<std::endl;
      if(i>0){
        TEST_ASSERT(atoms[indices[i]].atom->getAtomicNum() >= atoms[indices[i-1]].atom->getAtomicNum());
        if(atoms[indices[i]].atom->getAtomicNum() != atoms[indices[i-1]].atom->getAtomicNum()){
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
      atoms[i].atom = m->getAtomWithIdx(i);
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

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<atoms[order[i]].invar<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    // }

    TEST_ASSERT(order[0]==1);
    TEST_ASSERT(order[6]==0);
    TEST_ASSERT(order[7]==3);
    TEST_ASSERT(count[order[0]]==6);
    TEST_ASSERT(count[order[1]]==0);
    TEST_ASSERT(count[order[6]]==1);
    TEST_ASSERT(count[order[7]]==1);
  }
  {
    // this time with smarter invariants
    std::string smi="FC1C(Cl)CCC1C";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].atom=m->getAtomWithIdx(i);
      atoms[i].index=i;
    }
    atomcomparefunctor2 ftor(&atoms.front());

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    // }

    TEST_ASSERT(order[0]==7);
    TEST_ASSERT(order[1]==0);
    TEST_ASSERT(order[2]==3);
    TEST_ASSERT(order[3]==4);
    TEST_ASSERT(order[5]==1);
    TEST_ASSERT(count[order[0]]==1);
    TEST_ASSERT(count[order[1]]==1);
    TEST_ASSERT(count[order[2]]==1);
    TEST_ASSERT(count[order[3]]==2);
    TEST_ASSERT(count[order[4]]==0);
    TEST_ASSERT(count[order[5]]==3);
    TEST_ASSERT(count[order[6]]==0);

  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


class atomcomparefunctor3 {
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
      res+=(nbridx+1)*1000+nbrs[nbridx];
    }
    return res;
  }
  int basecomp(int i,int j) const {
    PRECONDITION(dp_atoms,"no atoms");
    unsigned int ivi,ivj;

    // always start with the current class:
    ivi= dp_atoms[i].index;
    ivj= dp_atoms[j].index;
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;
    
    // start by comparing degree
    ivi= dp_atoms[i].atom->getDegree();
    ivj= dp_atoms[j].atom->getDegree();
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    // move onto atomic number
    ivi= dp_atoms[i].atom->getAtomicNum();
    ivj= dp_atoms[j].atom->getAtomicNum();
    if(ivi<ivj)
      return -1;
    else if(ivi>ivj)
      return 1;

    return 0;
  }
public:
  bool df_useNbrs;
  atomcomparefunctor3() : dp_atoms(NULL), dp_mol(NULL), df_useNbrs(false) {};
  atomcomparefunctor3(Canon::canon_atom *atoms, const ROMol &m) : dp_atoms(atoms), dp_mol(&m),
                                                                  df_useNbrs(false) {};
  int operator()(int i,int j) const {
    PRECONDITION(dp_atoms,"no atoms");
    PRECONDITION(dp_mol,"no molecule");
    int v=basecomp(i,j);
    if(v) return v;
    unsigned int ivi,ivj;
    if(df_useNbrs){
      ivi=dp_atoms[i].index+1+getAtomNeighborhood(i);
      ivj=dp_atoms[j].index+1+getAtomNeighborhood(j);
      //std::cerr<<"               "<<i<<"-"<<j<<": "<<ivi<<" "<<ivj<<std::endl;
      if(ivi<ivj)
        return -1;
      else if(ivi>ivj)
        return 1;
    }
    return 0;
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
      atoms[i].atom=m->getAtomWithIdx(i);
      atoms[i].index=i;
    }
    atomcomparefunctor3 ftor(&atoms.front(),*m);
    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);

    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs=true;
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor,true,order,count,activeset,next);

    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]]==1);
      if(i>0){
        TEST_ASSERT(ftor(order[i],order[i-1])>=0);
      }
    }
    delete m;
  }

  {
    std::string smi="FC1C(CO)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].atom=m->getAtomWithIdx(i);
      atoms[i].index=i;
    }
    atomcomparefunctor3 ftor(&atoms.front(),*m);

    RDKit::Canon::canon_atom *data=&atoms.front();
    int *count=(int *)malloc(atoms.size()*sizeof(int));
    int *order=(int *)malloc(atoms.size()*sizeof(int));
    int activeset;
    int *next=(int *)malloc(atoms.size()*sizeof(int));
    RDKit::Canon::CreateSinglePartition(atoms.size(),order,count,data);
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);

    RDKit::Canon::RefinePartitions(*m,data,ftor,false,order,count,activeset,next);
    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    // }

    ftor.df_useNbrs=true;
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor,true,order,count,activeset,next);
    
    //std::cerr<<"----------------------------------"<<std::endl;
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]]==1);
      if(i>0){
        //std::cerr<<"  ftor: "<<ftor(order[i],order[i-1])<<std::endl;
        TEST_ASSERT(ftor(order[i],order[i-1])>=0);
      }
    }
    delete m;
  }

  {
    std::string smi="FC1C(CC)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);
    std::vector<Canon::canon_atom> atoms(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      atoms[i].atom=m->getAtomWithIdx(i);
      atoms[i].index=i;
    }
    atomcomparefunctor3 ftor(&atoms.front(),*m);

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

    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs=true;

    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor,true,order,count,activeset,next);
    //std::cerr<<"----------------------------------"<<std::endl;


    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      if(i>0){
        //std::cerr<<"  ftor: "<<ftor(order[i],order[i-1])<<std::endl;
        TEST_ASSERT(ftor(order[i],order[i-1])>=0);
      }
    }
    
    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0]==4 && count[4]==2); 
    TEST_ASSERT(order[1]==9 && count[9]==0); 
    TEST_ASSERT(order[2]==0 && count[0]==1); 
    TEST_ASSERT(order[3]==3 && count[3]==2); 
    TEST_ASSERT(order[4]==8 && count[8]==0); 
    TEST_ASSERT(order[5]==5 && count[5]==2); 
    TEST_ASSERT(order[6]==6 && count[6]==0); 
    TEST_ASSERT(order[7]==2 && count[2]==2); 
    TEST_ASSERT(order[8]==7 && count[7]==0); 
    TEST_ASSERT(order[9]==1 && count[1]==1); 

    delete m;
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
      atoms[i].atom=m->getAtomWithIdx(i);
      atoms[i].index=i;
    }
    atomcomparefunctor3 ftor(&atoms.front(),*m);

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

    //std::cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    ftor.df_useNbrs=true;
    RDKit::Canon::ActivatePartitions(atoms.size(),order,count,activeset,next);
    RDKit::Canon::RefinePartitions(*m,data,ftor,true,order,count,activeset,next);

    // std::cerr<<"----------------------------------"<<std::endl;
    // for(unsigned int i=0;i<m->getNumAtoms();++i){
    //   std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
    // }

    // here we can't manage to get everything  unique
    TEST_ASSERT(order[0]==4 && count[4]==2); 
    TEST_ASSERT(order[1]==9 && count[9]==0); 
    TEST_ASSERT(order[2]==0 && count[0]==1); 
    TEST_ASSERT(order[3]==3 && count[3]==2); 
    TEST_ASSERT(order[4]==8 && count[8]==0); 
    TEST_ASSERT(order[5]==5 && count[5]==2); 
    TEST_ASSERT(order[6]==6 && count[6]==0); 
    TEST_ASSERT(order[7]==2 && count[2]==2); 
    TEST_ASSERT(order[8]==7 && count[7]==0); 
    TEST_ASSERT(order[9]==1 && count[1]==1); 

    RDKit::Canon::BreakTies(*m,data,ftor,true,order,count,activeset,next);
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<order[i]<<" "<<" index: "<<atoms[order[i]].index<<" count: "<<count[order[i]]<<std::endl;
      TEST_ASSERT(count[order[i]]==1);
    }
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};


void test6(){
  BOOST_LOG(rdInfoLog) << "testing canonicalization using the wrapper." << std::endl;
  // canonicalization using the wrapper
#if 1
  {
    std::string smi="FC1C(CC)CCC1CC";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::RankMolAtoms(*m,atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i],1);
    }
    //std::copy(atomRanks.begin(),atomRanks.end(),std::ostream_iterator<unsigned int>(std::cerr," "));
    //std::cerr<<std::endl;
    TEST_ASSERT(atomRanks[0]==2);
    TEST_ASSERT(atomRanks[1]==7);
    TEST_ASSERT(atomRanks[2]==8);
    TEST_ASSERT(atomRanks[3]==3);
    TEST_ASSERT(atomRanks[4]==0);
    TEST_ASSERT(atomRanks[5]==5);
    TEST_ASSERT(atomRanks[6]==6);
    TEST_ASSERT(atomRanks[7]==9);
    TEST_ASSERT(atomRanks[8]==4);
    TEST_ASSERT(atomRanks[9]==1);
    delete m;
  }

  {
    std::string smi="CC[C@@H]1CCC[C@@H](O1)C(=O)O";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::RankMolAtoms(*m,atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i],1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //   std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    TEST_ASSERT(atomRanks[0]==0);
    TEST_ASSERT(atomRanks[1]==3);
    TEST_ASSERT(atomRanks[2]==9);
    TEST_ASSERT(atomRanks[3]==5);
    TEST_ASSERT(atomRanks[4]==4);
    TEST_ASSERT(atomRanks[5]==6);
    TEST_ASSERT(atomRanks[6]==10);
    TEST_ASSERT(atomRanks[7]==7);
    TEST_ASSERT(atomRanks[8]==8);
    TEST_ASSERT(atomRanks[9]==2);
    TEST_ASSERT(atomRanks[10]==1);

    delete m;
  }

  {
    std::string smi="N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](Cc3c[nH]c4ccccc34)C(=O)OCc5ccccc5";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::RankMolAtoms(*m,atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i],1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //   std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    delete m;
  }
#endif  
  {
    std::string smi="BrC=C1CCC(C(=O)O1)c2cccc3ccccc23";
    RWMol *m =SmilesToMol(smi);
    TEST_ASSERT(m);

    std::vector<unsigned int> atomRanks;
    RDKit::Canon::RankMolAtoms(*m,atomRanks);
    boost::dynamic_bitset<> seen(m->getNumAtoms());
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      //std::cerr<<i<<" "<<atomRanks[i]<<std::endl;
      TEST_ASSERT(!seen[atomRanks[i]]);
      seen.set(atomRanks[i],1);
    }

    // for(unsigned int ii=0;ii<atomRanks.size();++ii){
    //   std::cerr<<ii<<":"<<atomRanks[ii]<<std::endl;
    // }
    delete m;
  }


  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
};



int main(){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  return 0;
}

