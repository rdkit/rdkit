// $Id$
//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "SubstructUtils.h"
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>

#include <argedit.h>


namespace RDKit{

bool atomCompat(Atom const *a1,Atom const *a2){
  bool res = a1->Match(a2);
  //std::cout << "\t\tatomCompat: "<< a1->getIdx() << "-" << a2->getIdx() << ": " << res << std::endl;
  return res;
}

bool bondCompat(Bond const *b1,Bond const *b2){
  bool res = b1->Match(b2);
  //std::cout << "\t\tbondCompat: "<< b1->getIdx() << "-" << b2->getIdx() << ": " << res << std::endl;
  return res;
}

void MolToVFGraph(const ROMol &mol,ARGEdit *vgEd){
  ROMol::ConstAtomIterator atomIt;
  for(atomIt=mol.beginAtoms();atomIt!=mol.endAtoms();atomIt++){
    vgEd->InsertNode((void *)*atomIt);
  }
  ROMol::ConstBondIterator bondIt;
  for(bondIt=mol.beginBonds();bondIt!=mol.endBonds();bondIt++){
    Bond const *bond = *bondIt;
    int idx1=bond->getBeginAtomIdx(),idx2=bond->getEndAtomIdx();
    vgEd->InsertEdge(idx1,idx2,(void *)bond);
    // FIX: this maybe ought to be changed to include other dative bond types?
    if(bond->getBondType() != Bond::DATIVE){
      vgEd->InsertEdge(idx2,idx1,(void *)bond);
    }
  }
}

bool substructVisitor(int n, node_id ni1[],node_id ni2[],void *mvp)
{
  std::vector< MatchVectType > *matchV = (std::vector< MatchVectType > *)mvp;
  MatchVectType locV;
  locV.resize(n);
  for(int i=0;i<n;i++){
    locV[i] = std::pair<int,int>(ni1[i],ni2[i]);
  }
  matchV->push_back(locV);
  return false;
}

bool substructHeadVisitor(int n, node_id ni1[],node_id ni2[],void *mvp)
{
  std::vector< int > *matchV = (std::vector< int > *)mvp;
  for(int i=0;i<n;i++){
    if(ni1[i]==0){
      matchV->push_back(ni2[i]);
      break;
    }
  }
  return false;
}

double toPrime(const MatchVectType &v){
  double res = 1.0;
  MatchVectType::const_iterator ci;
  for(ci=v.begin();ci!=v.end();ci++){
    int idx=ci->second;;
    //std::cout << " " << idx;
    CHECK_INVARIANT(idx<NUM_PRIMES_AVAIL,"number too large");
    res *= firstThousandPrimes[idx];
  }
  return res;
}

void removeDuplicates(std::vector<MatchVectType> &v){
  //
  //  This works by calculating a product of primes based on the
  //  indices of the atoms in each match vector.  This can lead to
  //  unexpected behavior when looking at rings and queries that don't
  //  specify bond orders.  For example querying this molecule:
  //    C1CCC=1
  //  with the pattern constructed from SMARTS C~C~C~C will return a
  //  single match, despite the fact that there are 4 different paths
  //  when valence is considered.  The defense of this behavior is
  //  that the 4 paths are equivalent in the semantics of the query.
  //  Also, OELib returns the same results
  //
  DOUBLE_VECT seen;
  std::vector<MatchVectType> res;
  std::vector<MatchVectType>::iterator i;
  for(i=v.begin();i!=v.end();i++){
    //std::cout << "Path: ";
    double val = toPrime(*i);
    //std::cout << " " << val << std::endl;
    if(std::find(seen.begin(),seen.end(),val) == seen.end()){
      // it's something new
      //std::cout << "KEEP" << std::endl;
      res.push_back(*i);
      seen.push_back(val);
    }
  }
  v = res;
}

}
