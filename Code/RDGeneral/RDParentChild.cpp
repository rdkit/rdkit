//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDParentChild.h>

#ifdef RDK_THREADSAFE_SSS
#include <mutex>
#endif

#ifdef RDK_THREADSAFE_SSS
namespace {
std::mutex &propmutex_get() {
  // create on demand
  static std::mutex _mutex;
  return _mutex;
}

void propmutex_create() {
  std::mutex &mutex = propmutex_get();
  std::lock_guard<std::mutex> test_lock(mutex);
}

std::mutex &GetPropMutex() {
  static std::once_flag flag;
  std::call_once(flag, propmutex_create);
  return propmutex_get();
}
}  // namespace
#endif

namespace RDKit {

RDParent::~RDParent() {
#ifdef RDK_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  for (auto child : d_children) {
    child->unsetParent();
  }
}

bool RDParent::addChild(RDChild *child) {
  return d_children.insert(child).second;
}

bool RDParent::removeChild(RDChild *child) {
  auto it = d_children.find(child);
  bool res = (it != d_children.end());
  if (res) {
    d_children.erase(it);
  }
  return res;
}

RDChild::RDChild(RDParent *parent) : d_parent(parent) {
#ifdef RDK_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  if (d_parent) {
    d_parent->addChild(this);
  }
}

RDChild::~RDChild() {
#ifdef RDK_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(GetPropMutex());
#endif
  if (d_parent) {
    d_parent->removeChild(this);
  }
}

}  // namespace RDKit
