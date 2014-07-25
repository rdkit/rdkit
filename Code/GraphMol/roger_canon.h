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
      unsigned int invar;
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

          memcpy(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memcpy(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;
          memcpy(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memcpy(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len2;
        } else if( stat < 0 ) {
          assert(count[*s1]>0);
          int len = count[*s1];
          memcpy(ptr,s1,len*sizeof(int));
          ptr += len;  n1 -= len;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memcpy(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len;
        } else /* stat > 0 */ {
          assert(count[*s2]>0);
          int len = count[*s2];
          memcpy(ptr,s2,len*sizeof(int));
          ptr += len; n2 -= len;
          //std::cerr<<"       len: "<<len<<" "<<n2<<" "<<n1<<std::endl;
          if( n2 == 0 ) {
            memcpy(ptr,s1,n1*sizeof(int));
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
        int stat = compar(*s1,*s2);
        int len1 = count[*s1];
        int len2 = count[*s2];
        assert(len1>0);
        assert(len2>0);
        if( stat == 0 ) {
          count[*s1] = len1+len2;
          count[*s2] = 0;
          memcpy(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memcpy(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;
          memcpy(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memcpy(ptr,s1,n1*sizeof(int));
            return result;
          }
          s2 += len2;
        } else if( stat < 0 && len1>0 ) {
          memcpy(ptr,s1,len1*sizeof(int));
          ptr += len1;  n1 -= len1;
          if( n1 == 0 ) {
            if( ptr != s2 )
              memcpy(ptr,s2,n2*sizeof(int));
            return result;
          }
          s1 += len1;
        } else if (stat > 0  && len2>0) /* stat > 0 */ {
          memcpy(ptr,s2,len2*sizeof(int));
          ptr += len2; n2 -= len2;
          if( n2 == 0 ) {
            memcpy(ptr,s1,n1*sizeof(int));
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
              memcpy(ptr,s2,n2*sizeof(int));
            return result;
          }
        } else {
          *ptr++ = *s2++;
          n2--;
          if( n2 == 0 ) {
            memcpy(ptr,s1,n1*sizeof(int));
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
        memcpy(ptr,s1,n1*sizeof(int));
      } else if( ptr != s2 )
        memcpy(ptr,s2,n2*sizeof(int));
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
        memcpy(base,temp,nel*sizeof(int));
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


    class BaseAtomCompareFunctor {
      canon_atom *d_atoms;
    public:
      BaseAtomCompareFunctor() : d_atoms(NULL) {};
      BaseAtomCompareFunctor(canon_atom *atoms) : d_atoms(atoms) {};
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

    class NbrAtomCompareFunctor {
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
      NbrAtomCompareFunctor() : dp_atoms(NULL), dp_mol(NULL) {};
      NbrAtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m) : dp_atoms(atoms), dp_mol(&m) {};
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
        //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" "<<atoms[order[ii]].invar<<" index: "<<atoms[order[ii]].index<<std::endl;
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
        //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" "<<atoms[order[ii]].invar<<" index: "<<atoms[order[ii]].index<<std::endl;
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
