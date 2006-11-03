//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBSTRUCT_UTILS_H_
#define _RD_SUBSTRUCT_UTILS_H_

#include "SubstructMatch.h"

class ARGEdit;
namespace RDKit{
  class Atom;
  class Bond;
  class ROMol;
  
  bool atomCompat(Atom const *a1,Atom const *a2);
  bool chiralAtomCompat(Atom const *a1,Atom const *a2);
  bool bondCompat(Bond const *b1,Bond const *b2);
  void MolToVFGraph(const ROMol &mol,ARGEdit *vgEd);
  bool substructVisitor(int n, node_id ni1[],node_id ni2[],void *mvp);
  bool substructHeadVisitor(int n, node_id ni1[],node_id ni2[],void *mvp);
  double toPrime(const MatchVectType &v);
  void removeDuplicates(std::vector<MatchVectType> &v);

}


#endif
