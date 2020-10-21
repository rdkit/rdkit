//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_RDPARENTCHILD_H
#define RDKIT_RDPARENTCHILD_H

#include <unordered_set>

namespace RDKit {

class RDChild;

class RDKIT_RDGENERAL_EXPORT RDParent {
  friend class RDChild;

 public:
  RDParent() {
  }
  ~RDParent();

 private:
  bool addChild(RDChild *child);
  bool removeChild(RDChild *child);
  std::unordered_set<RDChild *> d_children;
};

class RDKIT_RDGENERAL_EXPORT RDChild {
  friend class RDParent;

 public:
  RDChild(RDParent *parent);
  virtual ~RDChild();

 protected:
  bool hasParent() const { return (d_parent != nullptr); }

 private:
  RDParent *parent() const { return d_parent; }
  void unsetParent() { d_parent = nullptr; }
  RDParent *d_parent;
};

}  // namespace RDKit
#endif
