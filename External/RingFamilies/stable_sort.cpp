#include "stable_sort.h"

#include <algorithm>

void stable_sort(void *arr, size_t n, size_t,
                 int (*comp)(const void *, const void *)) {
  std::stable_sort(
      (void **)arr, (void **)arr + n,
      [comp](const void *a, const void *b) { return comp(&a, &b) < 0; });
}
