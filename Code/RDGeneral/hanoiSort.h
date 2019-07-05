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

#include <RDGeneral/export.h>
#ifndef HANOISORT_H_
#define HANOISORT_H_

#include <cstring>
#include <iostream>
#include <cassert>
#include <cstdlib>

#if defined(_MSC_VER)
#pragma warning(push, 1)
#pragma warning(disable : 4800)
#endif
namespace RDKit {
template <typename CompareFunc>
bool hanoi(int *base, int nel, int *temp, int *count, int *changed,
           CompareFunc compar) {
  assert(base);
  assert(temp);
  assert(count);
  assert(changed);
  // std::cerr<<"  hanoi: "<<nel<< " start " << (*base)+1 << std::endl;
  int *b1, *b2;
  int *t1, *t2;
  int *s1, *s2;
  int n1, n2;
  int result;
  int *ptr;

  if (nel == 1) {
    count[base[0]] = 1;
    return false;
  } else if (nel == 2) {
    n1 = base[0];
    n2 = base[1];
    int stat =
        (/*!changed || */ changed[n1] || changed[n2]) ? compar(n1, n2) : 0;
    if (stat == 0) {
      count[n1] = 2;
      count[n2] = 0;
      return false;
    } else if (stat < 0) {
      count[n1] = 1;
      count[n2] = 1;
      return false;
    } else /* stat > 0 */ {
      count[n1] = 1;
      count[n2] = 1;
      base[0] = n2; /* temp[0] = n2; */
      base[1] = n1; /* temp[1] = n1; */
      return false; /* return True;  */
    }
  }

  n1 = nel / 2;
  n2 = nel - n1;
  b1 = base;
  t1 = temp;
  b2 = base + n1;
  t2 = temp + n1;

  if (hanoi(b1, n1, t1, count, changed, compar)) {
    if (hanoi(b2, n2, t2, count, changed, compar)) {
      s2 = t2;
    } else
      s2 = b2;
    result = false;
    ptr = base;
    s1 = t1;
  } else {
    if (hanoi(b2, n2, t2, count, changed, compar)) {
      s2 = t2;
    } else
      s2 = b2;
    result = true;
    ptr = temp;
    s1 = b1;
  }

  while (true) {
    assert(*s1 != *s2);
    int stat =
        (/*!changed || */ changed[*s1] || changed[*s2]) ? compar(*s1, *s2) : 0;
    int len1 = count[*s1];
    int len2 = count[*s2];
    assert(len1 > 0);
    assert(len2 > 0);
    if (stat == 0) {
      count[*s1] = len1 + len2;
      count[*s2] = 0;
      memmove(ptr, s1, len1 * sizeof(int));
      ptr += len1;
      n1 -= len1;
      if (n1 == 0) {
        if (ptr != s2) memmove(ptr, s2, n2 * sizeof(int));
        return result;
      }
      s1 += len1;

      // std::cerr<<"  cpy: "<<*s1<<" "<<*s2<<" "<<len2<<std::endl;
      memmove(ptr, s2, len2 * sizeof(int));
      ptr += len2;
      n2 -= len2;
      if (n2 == 0) {
        memmove(ptr, s1, n1 * sizeof(int));
        return result;
      }
      s2 += len2;
    } else if (stat < 0 && len1 > 0) {
      memmove(ptr, s1, len1 * sizeof(int));
      ptr += len1;
      n1 -= len1;
      if (n1 == 0) {
        if (ptr != s2) memmove(ptr, s2, n2 * sizeof(int));
        return result;
      }
      s1 += len1;
    } else if (stat > 0 && len2 > 0) /* stat > 0 */ {
      memmove(ptr, s2, len2 * sizeof(int));
      ptr += len2;
      n2 -= len2;
      if (n2 == 0) {
        memmove(ptr, s1, n1 * sizeof(int));
        return result;
      }
      s2 += len2;
    } else {
      assert(0);
    }
  }
}

template <typename CompareFunc>
void hanoisort(int *base, int nel, int *count, int *changed,
               CompareFunc compar) {
  assert(base);
  int *temp = (int *)malloc(nel * sizeof(int));
  assert(temp);
  if (hanoi(base, nel, temp, count, changed, compar))
    memmove(base, temp, nel * sizeof(int));
  free(temp);
}
}  // namespace RDKit

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

#endif /* HANOISORT_H_ */
