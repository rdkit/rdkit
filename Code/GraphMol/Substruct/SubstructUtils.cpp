// $Id$
//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstructUtils.h"
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <boost/dynamic_bitset.hpp>

namespace RDKit{

  bool atomCompat(const ATOM_SPTR &a1,const ATOM_SPTR &a2,bool useQueryQueryMatches){
    PRECONDITION(a1,"bad atom");
    PRECONDITION(a2,"bad atom");
    //std::cerr << "\t\tatomCompat: "<< a1 << " " << a1->getIdx() << "-" << a2 << " " << a2->getIdx() << std::endl;
    bool res;
    if(useQueryQueryMatches && a1->hasQuery() && a2->hasQuery()){
      res = static_cast<QueryAtom *>(a1.get())->QueryMatch(static_cast<QueryAtom *>(a2.get()));
    } else {
      res = a1->Match(a2);
    }
    return res;
  }

  bool chiralAtomCompat(const ATOM_SPTR &a1,const ATOM_SPTR &a2){
    PRECONDITION(a1,"bad atom");
    PRECONDITION(a2,"bad atom");
    //std::cerr << "\t\tatomCompat: "<< a1 << " " << a1->getIdx() << "-" << a2 << " " << a2->getIdx() << std::endl;
    bool res = a1->Match(a2);
    if(res){
      if(a1->hasProp("_CIPCode") || a2->hasProp("_CIPCode")){
        // if either atom has a CIPCode, they need to both have it and match:
        if(a1->hasProp("_CIPCode") && a2->hasProp("_CIPCode")){
          std::string s1,s2;
          a1->getProp("_CIPCode",s1);
          a2->getProp("_CIPCode",s2);
          if(s1!=s2) res=false;
        } else {
          res=false;
        }
      }
    }
    return res;
  }

  bool bondCompat(const BOND_SPTR &b1,const BOND_SPTR &b2,bool useQueryQueryMatches){
    PRECONDITION(b1,"bad bond");
    PRECONDITION(b2,"bad bond");
    bool res;
    if(useQueryQueryMatches && b1->hasQuery() && b2->hasQuery()){
      res = static_cast<QueryBond *>(b1.get())->QueryMatch(static_cast<QueryBond *>(b2.get()));
    } else {
      res = b1->Match(b2);
    }
    //std::cout << "\t\tbondCompat: "<< b1->getIdx() << "-" << b2->getIdx() << ": " << res << std::endl;
    return res;
  }

  void removeDuplicates(std::vector<MatchVectType> &v,unsigned int nAtoms){
    //
    //  This works by tracking the indices of the atoms in each match vector.  
    //  This can lead to unexpected behavior when looking at rings and queries 
    //  that don't specify bond orders.  For example querying this molecule:
    //    C1CCC=1
    //  with the pattern constructed from SMARTS C~C~C~C will return a
    //  single match, despite the fact that there are 4 different paths
    //  when valence is considered.  The defense of this behavior is
    //  that the 4 paths are equivalent in the semantics of the query.
    //  Also, OELib returns the same results
    //
    std::vector< boost::dynamic_bitset<> > seen;
    std::vector<MatchVectType> res;
    for(std::vector<MatchVectType>::const_iterator i=v.begin();i!=v.end();++i){
      boost::dynamic_bitset<> val(nAtoms);
      for(MatchVectType::const_iterator ci=i->begin();ci!=i->end();++ci){
        val.set(ci->second);
      }
      if(std::find(seen.begin(),seen.end(),val)==seen.end()){
        // it's something new
        res.push_back(*i);
        seen.push_back(val);
      }
    }
    v = res;
  }
}
