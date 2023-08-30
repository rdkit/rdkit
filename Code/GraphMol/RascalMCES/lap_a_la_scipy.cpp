// This is a mildly modified version of the code in SciPy's
// scipy.optimize.linear_sum_assignment, extracted from
// rectangular_lsap.cpp.
// https://github.com/scipy/scipy/blob/main/scipy/optimize/rectangular_lsap/rectangular_lsap.cpp
// As such it is subject to the following notice:
/*
Copyright (c) 2001-2002 Enthought, Inc. 2003-2023, SciPy Developers.
All rights reserved.

 Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


This code implements the shortest augmenting path algorithm for the
rectangular assignment problem.  This implementation is based on the
pseudocode described in pages 1685-1686 of:

    DF Crouse. On implementing 2D rectangular assignment algorithms.
    IEEE Transactions on Aerospace and Electronic Systems
    52(4):1679-1696, August 2016
    doi: 10.1109/TAES.2016.140952

Author: PM Larsen
*/

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace RDKit {
namespace RascalMCES {
template <typename T>
std::vector<size_t> argsortIter(const std::vector<T> &v) {
  std::vector<size_t> index(v.size());
  std::iota(index.begin(), index.end(), 0);
  std::sort(index.begin(), index.end(),
            [&v](size_t i, size_t j) { return v[i] < v[j]; });
  return index;
}

static int augmentingPath(size_t nc, std::vector<int> &cost,
                          std::vector<double> &u, std::vector<double> &v,
                          std::vector<size_t> &path,
                          std::vector<size_t> &row4col,
                          std::vector<double> &shortestPathCosts, size_t i,
                          std::vector<bool> &SR, std::vector<bool> &SC,
                          std::vector<size_t> &remaining, double *p_minVal) {
  double minVal = 0;

  // Crouse's pseudocode uses set complements to keep track of remaining
  // nodes.  Here we use a vector, as it is more efficient in C++.
  size_t numRemaining = nc;
  for (size_t it = 0; it < nc; it++) {
    // Filling this up in reverse order ensures that the solution of a
    // constant cost matrix is the identity matrix (c.f. #11602).
    remaining[it] = nc - it - 1;
  }

  std::fill(SR.begin(), SR.end(), false);
  std::fill(SC.begin(), SC.end(), false);
  std::fill(shortestPathCosts.begin(), shortestPathCosts.end(),
            std::numeric_limits<double>::max());

  // find shortest augmenting path
  int sink = -1;
  while (sink == -1) {
    // Clearly this will produce an overflow and set index to a large integer.
    // It is how the original code did it, and I assume whoever wrote it knew
    // what they were doing.
    size_t index = -1;
    double lowest = std::numeric_limits<double>::max();
    SR[i] = true;

    for (size_t it = 0; it < numRemaining; it++) {
      size_t j = remaining[it];

      double r = minVal + cost[i * nc + j] - u[i] - v[j];
      if (r < shortestPathCosts[j]) {
        path[j] = i;
        shortestPathCosts[j] = r;
      }

      // When multiple nodes have the minimum cost, we select one which
      // gives us a new sink node. This is particularly important for
      // integer cost matrices with small co-efficients.
      if (shortestPathCosts[j] < lowest ||
          (shortestPathCosts[j] == lowest &&
           row4col[j] == static_cast<size_t>(-1))) {
        lowest = shortestPathCosts[j];
        index = it;
      }
    }

    minVal = lowest;
    if (minVal ==
        std::numeric_limits<double>::max()) {  // infeasible cost matrix
      return -1;
    }

    size_t j = remaining[index];
    if (row4col[j] == static_cast<size_t>(-1)) {
      sink = j;
    } else {
      i = row4col[j];
    }

    SC[j] = true;
    remaining[index] = remaining[--numRemaining];
  }

  *p_minVal = minVal;
  return sink;
}

int lapMaximize(const std::vector<std::vector<int>> &costsMat,
                std::vector<size_t> &a, std::vector<size_t> &b) {
  if (costsMat.empty() || costsMat.front().empty()) {
    return 0;
  }
  size_t nr = costsMat.size();
  size_t nc = costsMat.front().size();
  bool transpose = nc < nr;
  std::vector<int> cost(nc * nr);
  // for maximization, take -ve of costs.
  for (size_t i = 0; i < nr; ++i) {
    for (size_t j = 0; j < nc; ++j) {
      if (transpose) {
        cost[j * nr + i] = -costsMat[i][j];
      } else {
        cost[i * nc + j] = -costsMat[i][j];
      }
    }
  }
  if (transpose) {
    std::swap(nc, nr);
  }
  // initialize variables
  std::vector<double> u(nr, 0);
  std::vector<double> v(nc, 0);
  std::vector<double> shortestPathCosts(nc);
  std::vector<size_t> path(nc, -1);
  std::vector<size_t> col4row(nr, -1);
  std::vector<size_t> row4col(nc, -1);
  std::vector<bool> SR(nr);
  std::vector<bool> SC(nc);
  std::vector<size_t> remaining(nc);

  // iteratively build the solution
  for (size_t curRow = 0; curRow < nr; curRow++) {
    double minVal;
    int sink = augmentingPath(nc, cost, u, v, path, row4col, shortestPathCosts,
                              curRow, SR, SC, remaining, &minVal);
    if (sink < 0) {
      return -1;
    }

    // update dual variables
    u[curRow] += minVal;
    for (size_t i = 0; i < nr; i++) {
      if (SR[i] && i != curRow) {
        u[i] += minVal - shortestPathCosts[col4row[i]];
      }
    }

    for (size_t j = 0; j < nc; j++) {
      if (SC[j]) {
        v[j] -= minVal - shortestPathCosts[j];
      }
    }

    // augment previous solution
    size_t j = sink;
    while (1) {
      size_t i = path[j];
      row4col[j] = i;
      std::swap(col4row[i], j);
      if (i == curRow) {
        break;
      }
    }
  }

  if (transpose) {
    size_t i = 0;
    for (auto v : argsortIter(col4row)) {
      a[i] = col4row[v];
      b[i] = v;
      i++;
    }
  } else {
    for (size_t i = 0; i < nr; i++) {
      a[i] = i;
      b[i] = col4row[i];
    }
  }

  return 0;
}
}  // namespace RascalMCES
}  // namespace RDKit
