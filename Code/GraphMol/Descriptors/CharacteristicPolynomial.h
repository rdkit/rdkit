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

//! |c_0| .. |c_n| of a characteristic polynomial, as scaled[k] * 2^exponent.
/*!
  The exponent is 0 unless the sum of the magnitudes would overflow a double.
  Then every magnitude is divided by the same power of two, which leaves their
  ratios, and so the entropy behind AvgIpc, unchanged.
*/
struct PolynomialMagnitudes {
  std::vector<double> scaled;  //!< leading coefficient first
  int exponent = 0;
};

//! |c_0| .. |c_n| of det(xI - A), each rounded to the nearest double once.
/*!
  A is the symmetric 0/1 matrix with ones at \c edges (pairs of distinct
  vertices, each listed once) and on the diagonal where \c loops is set. Its
  characteristic polynomial has integer coefficients, which are computed
  exactly: det(xI - A) is evaluated modulo primes below 2^62 at n + 1 points
  by an elimination ordered along the graph, interpolated, and the residues
  combined by the Chinese remainder theorem. A bound on the coefficients
  fixes how many primes that takes.

  \return n + 1 magnitudes; std::nullopt if the
          bound needs more primes than the built-in table holds (about 3,000
          atoms of a molecular graph, far beyond the double range)
*/
std::optional<PolynomialMagnitudes> characteristicPolynomialMagnitudes(
    unsigned int numVertices,
    const std::vector<std::pair<unsigned int, unsigned int>> &edges,
    const std::vector<bool> &loops);

}  // namespace detail
}  // namespace Descriptors
}  // namespace RDKit

#endif
