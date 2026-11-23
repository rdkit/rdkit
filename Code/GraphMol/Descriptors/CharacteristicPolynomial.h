//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Internal to the Descriptors library: the exact characteristic polynomial
// behind calcIpc().

#ifndef RD_CHARACTERISTICPOLYNOMIAL_H
#define RD_CHARACTERISTICPOLYNOMIAL_H

#include <optional>
#include <utility>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace detail {

//! |c_0| .. |c_n| of det(xI - A), each rounded to the nearest double once.
/*!
  A is the symmetric 0/1 matrix with ones at \c edges (pairs of distinct
  vertices, each listed once) and on the diagonal where \c loops is set. Its
  characteristic polynomial has integer coefficients, which are computed
  exactly: det(xI - A) is evaluated modulo primes below 2^62 at n + 1 points
  by an elimination ordered along the graph, interpolated, and the residues
  combined by the Chinese remainder theorem. A bound on the coefficients
  fixes how many primes that takes.

  \return n + 1 magnitudes, leading coefficient first; std::nullopt if the
          bound needs more primes than the built-in table holds (about 3,000
          atoms of a molecular graph, far beyond the double range)
*/
std::optional<std::vector<double>> characteristicPolynomialMagnitudes(
    unsigned int numVertices,
    const std::vector<std::pair<unsigned int, unsigned int>> &edges,
    const std::vector<bool> &loops);

}  // namespace detail
}  // namespace Descriptors
}  // namespace RDKit

#endif
