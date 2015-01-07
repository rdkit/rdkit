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
#include <GraphMol/RingInfo.h>
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cstring>
#include <iostream>
#include <cassert>
#include <vector>

// #define VERBOSE_CANON 1

namespace RDKit {
  namespace Canon{
    struct canon_atom{
      const Atom *atom;
      int index;
      unsigned int degree;
      bool hasRingNbr;
      const std::string *p_symbol; // if provided, this is used to order atoms
    };

    template <typename CompareFunc>
    bool hanoi( int* base, int nel, int *temp, int *count, int *changed, CompareFunc compar ) {
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
        int stat = (!changed || changed[n1]||changed[n2]) ? compar(n1,n2) : 0;
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

      if( hanoi(b1,n1,t1,count,changed,compar) ) {
        if( hanoi(b2,n2,t2,count,changed,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = false;
        ptr = base;
        s1 = t1;
      } else {
        if( hanoi(b2,n2,t2,count,changed,compar) ) {
          s2 = t2;
        } else s2 = b2;
        result = true;
        ptr = temp;
        s1 = b1;
      }

      while( true ) {
        assert(*s1!=*s2);
        int stat = (!changed || changed[*s1]||changed[*s2]) ? compar(*s1,*s2) : 0;
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
    }

    template <typename CompareFunc>
    void hanoisort( int* base, int nel, int *count, int *changed, CompareFunc compar )
    {
      register int *temp;

      temp = (int*)malloc(nel*sizeof(int));
      if( hanoi(base,nel,temp,count,changed,compar) )
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
                            int &activeset,int *next,
                            int *changed);    

    struct bondholder {
      Bond::BondType bondType;
      unsigned int bondStereo;
      unsigned int nbrIdx;
      bondholder() : bondType(Bond::UNSPECIFIED),
                     bondStereo(static_cast<unsigned int>(Bond::STEREONONE)),
                     nbrIdx(0) {};
      bondholder(Bond::BondType bt,Bond::BondStereo bs,unsigned int ni) : bondType(bt),
                                                                          bondStereo(static_cast<unsigned int>(bs)),
                                                                          nbrIdx(ni) {};
      bondholder(Bond::BondType bt,unsigned int bs,unsigned int ni) : bondType(bt),
                                                                   bondStereo(bs),
                                                                   nbrIdx(ni) {};
      bool operator<(const bondholder &o) const {
        if(bondType!=o.bondType)
          return bondType<o.bondType;
        if(bondStereo!=o.bondStereo)
          return bondStereo<o.bondStereo;
        return nbrIdx<o.nbrIdx;
      }

      static int compare(const bondholder &x,const bondholder &y,unsigned int div=1){
        if(x.bondType<y.bondType)
          return -1;
        else if(x.bondType>y.bondType)
          return 1;
        if(x.bondStereo<y.bondStereo)
          return -1;
        else if(x.bondStereo>y.bondStereo)
          return 1;

        return x.nbrIdx/div-y.nbrIdx/div;
      }

    };

    class AtomCompareFunctor {
      void getAtomNeighborhood(unsigned int i,std::vector<bondholder> &nbrs) const{
        nbrs.resize(0);
        if(dp_atomsInPlay && !(*dp_atomsInPlay)[i])
          return;
        const Atom *at=dp_atoms[i].atom;
        nbrs.reserve(at->getDegree());
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = dp_mol->getAtomBonds(at);
        // std::cerr<< "                   nbrs: "<<i+1;
        while(beg!=end){
          const BOND_SPTR bond=(*dp_mol)[*beg];
          ++beg;
          if(dp_bondsInPlay && !(*dp_bondsInPlay)[bond->getIdx()])
            continue;

          Bond::BondStereo stereo=Bond::STEREONONE;
          if(df_useChirality){
            stereo=bond->getStereo();
          }
          bondholder bh(bondholder(bond->getBondType(),stereo,
                                   dp_atoms[bond->getOtherAtomIdx(i)].index));
          nbrs.insert(std::lower_bound(nbrs.begin(),nbrs.end(),bh),1,bh);
          // std::cerr<<" "<<bond->getOtherAtomIdx(i)+1<<"("<<dp_atoms[bond->getOtherAtomIdx(i)].index<<")";
        }
        // std::cerr<<std::endl;
        std::reverse(nbrs.begin(),nbrs.end());
      }
      // EFF: it's stupid inefficient to be calling this frequently
      // however: since we use the current class of the nbr, there's some
      // heavy bookkeeping to really remove it
      unsigned int getAtomRingNbrCode(unsigned int i) const {
        if(!dp_atoms[i].hasRingNbr) return 0;

        const Atom *at=dp_atoms[i].atom;
        ROMol::ADJ_ITER beg,end;
        boost::tie(beg,end) = dp_mol->getAtomNeighbors(at);
        while(beg!=end){
          const ATOM_SPTR nbr=(*dp_mol)[*beg];
          ++beg;
          if((nbr->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
              nbr->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) &&
             nbr->hasProp("_ringStereoAtoms")){
            // we have a neighbor with stereochem
            unsigned int nidx=1;
            ROMol::ADJ_ITER nbeg,nend;
            boost::tie(nbeg,nend) = dp_mol->getAtomNeighbors(nbr);
            while(nbeg!=nend){
              const ATOM_SPTR nnbr=(*dp_mol)[*nbeg];
              if(nnbr->getIdx()==at->getIdx()){
                //std::cerr<<"  "<<at->getIdx()<<"-"<<nbr->getIdx()<<" "<<nidx<<std::endl;
                return dp_atoms[nbr->getIdx()].index*1000000+nidx;
              }
              ++nidx;
              ++nbeg;
            }
          }
        }
        return 0;
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
        ivi= dp_atoms[i].degree;
        ivj= dp_atoms[j].degree;
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;

        if(dp_atoms[i].p_symbol &&
           dp_atoms[j].p_symbol ) {
          if( *(dp_atoms[i].p_symbol) < *(dp_atoms[j].p_symbol) )
            return -1;
          else if( *(dp_atoms[i].p_symbol) > *(dp_atoms[j].p_symbol) )
            return 1;
          else
            return 0;
        }
        // move onto atomic number
        ivi= dp_atoms[i].atom->getAtomicNum();
        ivj= dp_atoms[j].atom->getAtomicNum();
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;

        // isotopes if we're using them
        if(df_useIsotopes){
          ivi=dp_atoms[i].atom->getIsotope();
          ivj=dp_atoms[j].atom->getIsotope();
          if(ivi<ivj)
            return -1;
          else if(ivi>ivj)
            return 1;
        }

        // nHs
        ivi=dp_atoms[i].atom->getTotalNumHs();
        ivj=dp_atoms[j].atom->getTotalNumHs();
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;
        
        // charge
        ivi=dp_atoms[i].atom->getFormalCharge();
        ivj=dp_atoms[j].atom->getFormalCharge();
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;
        
        // ring membership
        // initial passes at this were just checking "isInRing" to allow
        // a possible future more efficient check. These break on this
        // lovely double-diamond pathological case:
        //   *12*3*1*3*4*5*4*52
        //
        // probably going to regret allowing this to be skipped some day
        if(dp_mol->getRingInfo()->isInitialized()){ 
          ivi=dp_mol->getRingInfo()->numAtomRings(dp_atoms[i].atom->getIdx());
          ivj=dp_mol->getRingInfo()->numAtomRings(dp_atoms[j].atom->getIdx());
          if(ivi<ivj)
            return -1;
          else if(ivi>ivj)
            return 1;
        }
        
        // chirality if we're using it
        if(df_useChirality){
          // first atom stereochem:
          ivi=0;
          ivj=0;
          if(dp_atoms[i].atom->hasProp("_CIPCode")){
            std::string cipCode;
            dp_atoms[i].atom->getProp("_CIPCode",cipCode);
            ivi=cipCode=="R"?2:1;
          }
          if(dp_atoms[j].atom->hasProp("_CIPCode")){
            std::string cipCode;
            dp_atoms[j].atom->getProp("_CIPCode",cipCode);
            ivj=cipCode=="R"?2:1;
          }
          if(ivi<ivj)
            return -1;
          else if(ivi>ivj)
            return 1;

          // can't actually use values here, because they are arbitrary
          ivi=dp_atoms[i].atom->getChiralTag()!=0;
          ivj=dp_atoms[j].atom->getChiralTag()!=0;
          if(ivi<ivj)
            return -1;
          else if(ivi>ivj)
            return 1;
          
          // ring stereochemistry
          if(dp_mol->getRingInfo()->isInitialized() &&
             (dp_mol->getRingInfo()->numAtomRings(dp_atoms[i].atom->getIdx()) &&
              dp_mol->getRingInfo()->numAtomRings(dp_atoms[j].atom->getIdx()) ) ){
            ivi=getAtomRingNbrCode(i);
            ivj=getAtomRingNbrCode(j);
            if(ivi<ivj)
              return -1;
            else if(ivi>ivj)
              return 1;
          }
          // bond stereo is taken care of in the neighborhood comparison
        }
        return 0;
      }
    public:
      Canon::canon_atom *dp_atoms;
      const ROMol *dp_mol;
      const boost::dynamic_bitset<> *dp_atomsInPlay,*dp_bondsInPlay;
      bool df_useNbrs;
      bool df_useIsotopes;
      bool df_useChirality;
      AtomCompareFunctor() : dp_atoms(NULL), dp_mol(NULL),
                             dp_atomsInPlay(NULL), dp_bondsInPlay(NULL),
                             df_useNbrs(false),
                             df_useIsotopes(true), df_useChirality(true) {
      };
      AtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m,
                         const boost::dynamic_bitset<> *atomsInPlay=NULL,
                         const boost::dynamic_bitset<> *bondsInPlay=NULL ) :
        dp_atoms(atoms), dp_mol(&m),
        dp_atomsInPlay(atomsInPlay),dp_bondsInPlay(bondsInPlay),
        df_useNbrs(false),df_useIsotopes(true), df_useChirality(true) {
      };
      int operator()(int i,int j) const {
        PRECONDITION(dp_atoms,"no atoms");
        PRECONDITION(dp_mol,"no molecule");
        PRECONDITION(i!=j,"bad call");
        if(dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])){
          return 0;
        }
        int v=basecomp(i,j);
        // std::cerr<<"           bc: "<<i+1<<"-"<<j+1<<": "<<v<<std::endl;
        if(v) return v;

        if(df_useNbrs){
          std::vector<bondholder> nbrsi,nbrsj;
          getAtomNeighborhood(i,nbrsi);
          getAtomNeighborhood(j,nbrsj);
          for(unsigned int ii=0;ii<nbrsi.size() && ii<nbrsj.size();++ii){
            int cmp=bondholder::compare(nbrsi[ii],nbrsj[ii]);
            // std::cerr<<"              : "<<nbrsi[ii].nbrIdx<<"-"<<nbrsj[ii].nbrIdx<<": "<<cmp<<std::endl;
            if(cmp) return cmp;
          }
          if(nbrsi.size()<nbrsj.size()){
            return -1;
          } else if(nbrsi.size()>nbrsj.size()) {
            return 1;
          }
        }
        return 0;
      }
    };

    const unsigned int ATNUM_CLASS_OFFSET=10000;
    class ChiralAtomCompareFunctor {
      void getAtomNeighborhood(unsigned int i,std::vector<bondholder> &nbrs) const{
        unsigned int res=0;
        const Atom *at=dp_atoms[i].atom;
        nbrs.resize(0);
        nbrs.reserve(8);
        //std::cerr<<"gac: "<<i<<" "<<at->getTotalNumHs()<<std::endl;
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = dp_mol->getAtomBonds(at);
        while(beg!=end){
          const BOND_SPTR bond=(*dp_mol)[*beg++];
          unsigned int nbrIdx=bond->getOtherAtomIdx(i);
          const Atom *nbr=dp_mol->getAtomWithIdx(nbrIdx);
          //std::cerr<<"    "<<bond->getOtherAtom(at)->getIdx()<<std::endl;
          unsigned int stereo=0;
          switch(bond->getStereo()){
          case Bond::STEREOZ:
            stereo=1;
            break;
          case Bond::STEREOE:
            stereo=2;
            break;
          default:
            stereo=0;
          }
          unsigned int nReps;

          if(bond->getBondType()==Bond::DOUBLE &&
             nbr->getAtomicNum()==15 &&
             (nbr->getDegree()==4 || nbr->getDegree()==3) ) {
            // a special case for chiral phophorous compounds
            // (this was leading to incorrect assignment of
            // R/S labels ):
            nReps=1;

            // general justification of this is:
            // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
            // Multiple-Bond Unsaturation and Aromaticity". It contains several
            // conventions of which convention (b) is the one applying here:
            // "(b) Contibutions by d orbitals to bonds of quadriligant atoms are
            // neglected."
            // FIX: this applies to more than just P
          } else {
            nReps = static_cast<unsigned int>(floor(2.*bond->getBondTypeAsDouble()));
          }
          //std::cerr<<"          "<<nReps<<std::endl;
          // we push on both the atomic number and the index using the ATNUM_CLASS_OFFSET
          // pre-multiplier. This allows us to first compare the full set of neighbors based
          // on atomic number and then based on the ranks that have already been computed.
          bondholder bh(Bond::SINGLE,stereo,
                         nbr->getAtomicNum()*ATNUM_CLASS_OFFSET+
                         dp_atoms[nbrIdx].index+1);
          std::vector<bondholder>::iterator iPos=std::lower_bound(nbrs.begin(),nbrs.end(),bh);
          nbrs.insert(iPos,nReps,bh);
        }

        // FIX: don't want to be doing this long-term
        std::reverse(nbrs.begin(),nbrs.end());

        // add Hs that aren't in the graph, we know they go at the end:
        if(!at->needsUpdatePropertyCache()){
          for(unsigned int ii=0;ii<at->getTotalNumHs();++ii){
            nbrs.push_back(bondholder(Bond::SINGLE,Bond::STEREONONE,ATNUM_CLASS_OFFSET));
            nbrs.push_back(bondholder(Bond::SINGLE,Bond::STEREONONE,ATNUM_CLASS_OFFSET));
          }
        }
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

        // move onto atomic number
        ivi= dp_atoms[i].atom->getAtomicNum();
        ivj= dp_atoms[j].atom->getAtomicNum();
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;

        // isotopes:
        ivi=dp_atoms[i].atom->getIsotope();
        ivj=dp_atoms[j].atom->getIsotope();
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;

        // atom stereochem:
        ivi=0;
        ivj=0;
        if(dp_atoms[i].atom->hasProp("_CIPCode")){
          std::string cipCode;
          dp_atoms[i].atom->getProp("_CIPCode",cipCode);
          ivi=cipCode=="R"?2:1;
        }
        if(dp_atoms[j].atom->hasProp("_CIPCode")){
          std::string cipCode;
          dp_atoms[j].atom->getProp("_CIPCode",cipCode);
          ivj=cipCode=="R"?2:1;
        }
        if(ivi<ivj)
          return -1;
        else if(ivi>ivj)
          return 1;

        // bond stereo is taken care of in the neighborhood comparison
        return 0;
      }
    public:
      Canon::canon_atom *dp_atoms;
      const ROMol *dp_mol;
      bool df_useNbrs;
      ChiralAtomCompareFunctor() : dp_atoms(NULL), dp_mol(NULL), df_useNbrs(false) {
      };
      ChiralAtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m) : dp_atoms(atoms), dp_mol(&m), df_useNbrs(false) {
      };
      int operator()(int i,int j) const {
        PRECONDITION(dp_atoms,"no atoms");
        PRECONDITION(dp_mol,"no molecule");
        PRECONDITION(i!=j,"bad call");
        int v=basecomp(i,j);
        if(v) return v;

        if(df_useNbrs){
          std::vector<bondholder> nbrsi,nbrsj;
          getAtomNeighborhood(i,nbrsi);
          getAtomNeighborhood(j,nbrsj);

          // we do two passes through the neighbor lists. The first just uses the
          // atomic numbers (by passing the optional 10000 to bondholder::compare),
          // the second takes the already-computed index into account
          for(unsigned int ii=0;ii<nbrsi.size() && ii<nbrsj.size();++ii){
            int cmp=bondholder::compare(nbrsi[ii],nbrsj[ii],ATNUM_CLASS_OFFSET);
            if(cmp) return cmp;
          }
          for(unsigned int ii=0;ii<nbrsi.size() && ii<nbrsj.size();++ii){
            int cmp=bondholder::compare(nbrsi[ii],nbrsj[ii]);
            if(cmp) return cmp;
          }
          if(nbrsi.size()<nbrsj.size()){
            return -1;
          } else if(nbrsi.size()>nbrsj.size()) {
            return 1;
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
                          int &activeset, int *next,
                          int *changed){
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
        // std::cerr<<"\n\n**************************************************************"<<std::endl;
        // std::cerr<<"  sort: "<<atoms[partition].index<<" "<<len<<std::endl;
        //  for(unsigned int ii=0;ii<nAtoms;++ii){
        //    std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index: "<<atoms[order[ii]].index<<std::endl;
        //  }
        hanoisort(start,len,count,changed, compar);
        // std::cerr<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<std::endl;
        // for(unsigned int ii=0;ii<nAtoms;++ii){
        //    std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index: "<<atoms[order[ii]].index<<std::endl;
        //  }
        index = start[0];
        // std::cerr<<"  len:"<<len<<" index:"<<index<<" count:"<<count[index]<<std::endl;
        for( i=count[index]; i<len; i++ ) {
          index = start[i];
          if( count[index] )
            symclass = offset+i;
          atoms[index].index = symclass;
        }
        if( mode ) {
          index=start[0];
          for( i=count[index]; i<len; i++ ) {
            index=start[i];
            ROMol::ADJ_ITER nbrIdx,endNbrs;
            boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(atoms[index].atom);
            while(nbrIdx!=endNbrs){
              int nbor=mol[*nbrIdx]->getIdx();
              ++nbrIdx;
              int nbroffset = atoms[nbor].index;
              changed[nbor]=1;
              partition = order[nbroffset];
              // std::cerr<<"            changed: "<<index<<" "<<nbor<<" "<<nbroffset<<" count["<<partition<<"]:"<<count[partition]<<" "<<next[partition]<<std::endl;
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
                   int &activeset, int *next,
                   int *changed){
      unsigned int nAtoms=mol.getNumAtoms();
      register int partition;
      register int *start;
      register int offset;
      register int index;
      register int len;

      for(unsigned int i=0; i<nAtoms; i++ ) {
        partition = order[i];
        while( count[partition] > 1 ) {
          len = count[partition];
          offset = i+len-1;
          index = order[offset];
          atoms[index].index = offset;
          count[partition] = len-1;
          count[index] = 1;

          ROMol::ADJ_ITER nbrIdx,endNbrs;
          boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(atoms[index].atom);
          while(nbrIdx!=endNbrs){
            int nbor=mol[*nbrIdx]->getIdx();
            ++nbrIdx;

            offset = atoms[nbor].index;
            changed[nbor]=1;
            int npart = order[offset];
            if( (count[npart]>1) &&
                (next[npart]==-2) ) {
              next[npart] = activeset;
              activeset = npart;
            }
          }
          RefinePartitions(mol,atoms,compar,mode,order,count,activeset,next,changed);
        } 
      }
    } // end of BreakTies()

    void rankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res,
                      bool breakTies=true,
                      bool includeChirality=true,
                      bool includeIsotopes=true);
    void rankFragmentAtoms(const ROMol &mol,std::vector<unsigned int> &res, 
                           const boost::dynamic_bitset<> &atomsInPlay,
                           const boost::dynamic_bitset<> &bondsInPlay,
                           const std::vector<std::string> *atomSymbols=NULL,
                           bool breakTies=true,
                           bool includeChirality=true,
                           bool includeIsotopes=true);
    void chiralRankMolAtoms(const ROMol &mol,std::vector<unsigned int> &res);

  } // end of Canon namespace
} // end of RDKit namespace
