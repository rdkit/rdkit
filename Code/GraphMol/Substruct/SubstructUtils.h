//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBSTRUCT_UTILS_H_
#define _RD_SUBSTRUCT_UTILS_H_

#include "SubstructMatch.h"
#include <boost/smart_ptr.hpp>

#ifdef USE_VFLIB
class ARGEdit;
#endif

namespace RDKit{
  class ROMol;
  class Atom;
  class Bond;
  typedef boost::shared_ptr<Atom>    ATOM_SPTR;
  typedef boost::shared_ptr<Bond>    BOND_SPTR;
  
  double toPrime(const MatchVectType &v);
  void removeDuplicates(std::vector<MatchVectType> &v);
#ifdef USE_VFLIB
  bool atomCompat(const Atom *a1,const Atom *a2);
  bool chiralAtomCompat(const Atom *a1,const Atom *a2);
  bool bondCompat(const Bond *b1,const Bond *b2);
  void MolToVFGraph(const ROMol &mol,ARGEdit *vgEd);
  bool substructVisitor(int n, node_id ni1[],node_id ni2[],void *mvp);
  bool substructHeadVisitor(int n, node_id ni1[],node_id ni2[],void *mvp);
#else
  bool atomCompat(const ATOM_SPTR a1,const ATOM_SPTR a2);
  bool chiralAtomCompat(const ATOM_SPTR a1,const ATOM_SPTR a2);
  bool bondCompat(const BOND_SPTR b1,const BOND_SPTR b2);
#endif

}


#endif
