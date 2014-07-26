//
//  Copyright (C) 2014 Greg Landrum
//  Adapted from pseudo-code from Roger Sayle
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/ROMol.h>
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include <iostream>
#include <cassert>


namespace RDKit {
  namespace Canon{
    struct canon_atom{
      const Atom *atom;
      int index;
    };


#if 0    
    bool hanoi( int* base, int nel, int *temp, int *count, func compar ) {
      //std::cerr<<"  hanoi: "<<nel<<std::endl;
      register int *b1,*b2;
      register int *t1,*t2;
      register int *s1,*s2;
      register int n1,n2;
      register bool result;
      register int *ptr;

      if( nel == 1 ) {
        count[base[0]] = 1;
        return false;
      } else if( nel == 2 ) {
        n1 = base[0];
        n2 = base[1];
        int stat = (*compar)(n1,n2);
        if( stat == 0 ) {
          count[n1] = 2;
          count[n2] = 0;
          return false;
        } else if( stat < 0 ) {
          count[n1] = 1;
          count[n2] = 1;
          return false;
        } else /* stat > 0 */ {
          count[n1] = 1;
          count[n2] = 1;
          base[0] = n2;   /* temp[0] = n2; */
          base[1] = n1;   /* temp[1] = n1; */
          return false;   /* return True;  */
        }
      }

      n1 = nel/2;    n2 = nel-n1;
      b1 = base;     t1 = temp;
      b2 = base+n1;  t2 = temp+n1;

      if( hanoi(b1,n1,t1,count,compar) ) {
        if( hanoi(b2,n2,t2,count,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = false;
        ptr = base;
        s1 = t1;
      } else {
        if( hanoi(b2,n2,t2,count,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = true;
        ptr = temp;
        s1 = b1;
      }

      while( true ) {
        int stat = (*compar)(*s1,*s2);
        //std::cerr<<"    while: "<<*s1<<"-"<<*s2<<" "<<stat<<" ("<<s1<<","<<s2<<")"<<std::endl;
        if( stat == 0 ) {
          assert(count[*s1]>0);
          assert(count[*s2]>0);
          int len1 = count[*s1];
          int len2 = count[*s2];
          count[*s1] = len1+len2;
          count[*s2] = 0;

          memmove(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memmove(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;
          memmove(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memmove(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len2;
        } else if( stat < 0 ) {
          assert(count[*s1]>0);
          int len = count[*s1];
          memmove(ptr,s1,len*sizeof(int));
          ptr += len;  n1 -= len;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memmove(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len;
        } else /* stat > 0 */ {
          assert(count[*s2]>0);
          int len = count[*s2];
          memmove(ptr,s2,len*sizeof(int));
          ptr += len; n2 -= len;
          //std::cerr<<"       len: "<<len<<" "<<n2<<" "<<n1<<std::endl;
          if( n2 == 0 ) {
            memmove(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len;
        }
      }
    }
#else
    template <typename CompareFunc>
    bool hanoi( int* base, int nel, int *temp, int *count, CompareFunc compar ) {
      //std::cerr<<"  hanoi: "<<nel<<std::endl;
      register int *b1,*b2;
      register int *t1,*t2;
      register int *s1,*s2;
      register int n1,n2;
      register int result;
      register int *ptr;



      if( nel == 1 ) {
        count[base[0]] = 1;
        return false;
      } else if( nel == 2 ) {
        n1 = base[0];
        n2 = base[1];
        int stat = compar(n1,n2);
        if( stat == 0 ) {
          count[n1] = 2;
          count[n2] = 0;
          return false;
        } else if( stat < 0 ) {
          count[n1] = 1;
          count[n2] = 1;
          return false;
        } else /* stat > 0 */ {
          count[n1] = 1;
          count[n2] = 1;
          base[0] = n2;   /* temp[0] = n2; */
          base[1] = n1;   /* temp[1] = n1; */
          return false;   /* return True;  */
        }
      }

      n1 = nel/2;    n2 = nel-n1;
      b1 = base;     t1 = temp;
      b2 = base+n1;  t2 = temp+n1;

      if( hanoi(b1,n1,t1,count,compar) ) {
        if( hanoi(b2,n2,t2,count,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = false;
        ptr = base;
        s1 = t1;
      } else {
        if( hanoi(b2,n2,t2,count,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = true;
        ptr = temp;
        s1 = b1;
      }

#if 1
      while( true ) {
        assert(*s1!=*s2);
        int stat = compar(*s1,*s2);
        int len1 = count[*s1];
        int len2 = count[*s2];
        assert(len1>0);
        assert(len2>0);
        if( stat == 0 ) {
          count[*s1] = len1+len2;
          count[*s2] = 0;
          memmove(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memmove(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;

          //std::cerr<<"  cpy: "<<*s1<<" "<<*s2<<" "<<len2<<std::endl;
          memmove(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memmove(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len2;
        } else if( stat < 0 && len1>0 ) {
          memmove(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memmove(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;
        } else if (stat > 0  && len2>0) /* stat > 0 */ {
          memmove(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memmove(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len2;
        } else {
          assert(0);
        }
      }
#elif 1
      while( true ) {
        if( compar(*s1,*s2) <= 0 ) {
          *ptr++ = *s1++;
          n1--;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memmove(ptr,s2,n2*sizeof(int));
            return result;
          }
        } else {
          *ptr++ = *s2++;
          n2--;
          if( n2 == 0 ) {
            memmove(ptr,s1,n1*sizeof(int));
            return result;
          }
        }
      }
#else
      do {
        if( compar(*s1,*s2) <= 0 ) {
          *ptr++ = *s1++;
          n1--;
        } else {
          *ptr++ = *s2++;
          n2--;
        }
      } while( (n1>0) && (n2>0) );

      if( n1 > 0 ) {
        memmove(ptr,s1,n1*sizeof(int));
      } else if( ptr != s2 )
        memmove(ptr,s2,n2*sizeof(int));
      return result;
#endif    
    }
#endif

    template <typename CompareFunc>
    void hanoisort( int* base, int nel, int *count, CompareFunc compar )
    {
      register int *temp;

      temp = (int*)malloc(nel*sizeof(int));
      if( hanoi(base,nel,temp,count,compar) )
        memmove(base,temp,nel*sizeof(int));
      free(temp);
    }

    void CreateSinglePartition(unsigned int nAtoms,
                               int *order,
                               int *count,
                               canon_atom *atoms);

    void ActivatePartitions(unsigned int nAtoms,
                            int *order,
                            int *count,
                            int &activeset,int *next);    


    struct bondholder {
      Bond::BondType bondType;
      unsigned int nbrIdx;
      bondholder() : bondType(Bond::UNSPECIFIED), nbrIdx(0) {};
      bondholder(Bond::BondType bt,unsigned int ni) : bondType(bt),nbrIdx(ni) {};
      bool operator<(const bondholder &o) const {
        if(bondType!=o.bondType)
          return bondType<o.bondType;
        return nbrIdx<o.nbrIdx;
      }
      static int compare(const bondholder &x,const bondholder &y){
        if(x.bondType<y.bondType)
          return -1;
        else if(x.bondType>y.bondType)
          return 1;
        return x.nbrIdx-y.nbrIdx;
      }
    };
    class AtomCompareFunctor {
      Canon::canon_atom *dp_atoms;
      const ROMol *dp_mol;
      void getAtomNeighborhood(unsigned int i,std::vector<bondholder> &nbrs) const{
        unsigned int res=0;
        const Atom *at=dp_mol->getAtomWithIdx(i);
        nbrs.resize(0);
        nbrs.reserve(at->getDegree());
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = dp_mol->getAtomBonds(at);
        while(beg!=end){
          const BOND_SPTR bond=(*dp_mol)[*beg];
          nbrs.push_back(bondholder(bond->getBondType(),dp_atoms[bond->getOtherAtomIdx(i)].index));
          ++beg;
        }
        std::sort(nbrs.begin(),nbrs.end());
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
      AtomCompareFunctor() : dp_atoms(NULL), dp_mol(NULL), df_useNbrs(false) {};
      AtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m) : dp_atoms(atoms), dp_mol(&m), df_useNbrs(false) {};
      int operator()(int i,int j) const {
        PRECONDITION(dp_atoms,"no atoms");
        PRECONDITION(dp_mol,"no molecule");
        PRECONDITION(i!=j,"bad call");
        int v=basecomp(i,j);
        if(v) return v;
        TEST_ASSERT(i!=j);

        if(df_useNbrs){
          std::vector<bondholder> nbrsi,nbrsj;
          getAtomNeighborhood(i,nbrsi);
          getAtomNeighborhood(j,nbrsj);
          TEST_ASSERT(i!=j);
          for(unsigned int ii=0;ii<nbrsi.size();++ii){
            int cmp=bondholder::compare(nbrsi[ii],nbrsj[ii]);
            if(cmp) return cmp;
          }
        }
        return 0;
      }
    };


    template <typename CompareFunc>
    void RefinePartitions(const ROMol &mol,
                          canon_atom *atoms,
                          CompareFunc compar, int mode,
                          int *order,
                          int *count,
                          int &activeset, int *next){
      unsigned int nAtoms=mol.getNumAtoms();
      register int partition;
      register int symclass;
      register int *start;
      register int offset;
      register int index;
      register int len;
      register int i;

      while( activeset != -1 ) {
        //std::cerr<<"ITER: "<<activeset<<" next: "<<next[activeset]<<std::endl;
        // std::cerr<<" next: ";
        // for(unsigned int ii=0;ii<nAtoms;++ii){
        //   std::cerr<<ii<<":"<<next[ii]<<" ";
        // }
        // std::cerr<<std::endl; 
        // for(unsigned int ii=0;ii<nAtoms;++ii){
        //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index: "<<atoms[order[ii]].index<<std::endl;
        // }

        partition = activeset;
        activeset = next[partition];
        next[partition] = -2;

        len = count[partition]; 
        offset = atoms[partition].index;
        start = order+offset;
        //std::cerr<<"  sort: "<<atoms[partition].index<<" "<<len<<std::endl;
        hanoisort(start,len,count,compar);
        // for(unsigned int ii=0;ii<nAtoms;++ii){
        //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index: "<<atoms[order[ii]].index<<std::endl;
        // }

        index = start[0];
        for( i=count[index]; i<len; i++ ) {
          index = start[i];
          if( count[index] )
            symclass = offset+i;
          atoms[index].index = symclass;

          if( mode ) {
            ROMol::ADJ_ITER nbrIdx,endNbrs;
            boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(index));
            while(nbrIdx!=endNbrs){
              int nbor=mol[*nbrIdx]->getIdx();
              ++nbrIdx;
              int nbroffset = atoms[nbor].index;
              partition = order[nbroffset];
              if( (count[partition]>1) &&
                  (next[partition]==-2) ) {
                next[partition] = activeset;
                activeset = partition;
              }
            }
          }
        }
      }
    } // end of RefinePartitions()

    template <typename CompareFunc>
    void BreakTies(const ROMol &mol,
                   canon_atom *atoms,
                   CompareFunc compar, int mode,
                   int *order,
                   int *count,
                   int &activeset, int *next){
      unsigned int nAtoms=mol.getNumAtoms();
      register int partition;
      register int *start;
      register int offset;
      register int index;
      register int len;
      register int i;

      for( i=0; i<nAtoms; i++ ) {
        partition = order[i];
        while( count[partition] > 1 ) {
          len = count[partition];
          offset = i+len-1;
          index = order[offset];
          atoms[index].index = offset;
          count[partition] = len-1;
          count[index] = 1;

          ROMol::ADJ_ITER nbrIdx,endNbrs;
          boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(index));
          while(nbrIdx!=endNbrs){
            int nbor=mol[*nbrIdx]->getIdx();
            ++nbrIdx;

            offset = atoms[nbor].index;
            int npart = order[offset];
            if( (count[npart]>1) &&
                (next[npart]==-2) ) {
              next[npart] = activeset;
              activeset = npart;
            }
          }
          RefinePartitions(mol,atoms,compar,mode,order,count,activeset,next);
        } 
      }
    } // end of BreakTies()

    void RankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                      bool includeChirality=true,
                      bool includeIsotopes=true);

  } // end of Canon namespace
} // end of RDKit namespace
