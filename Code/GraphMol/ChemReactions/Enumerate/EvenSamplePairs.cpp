//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "EvenSamplePairs.h"
#include <boost/format.hpp>
#include <cstdint>

namespace RDKit {

using namespace EnumerationTypes;
// Based on an implementation from a correspondance with Bernd Rohde.
void EvenSamplePairsStrategy::initializeStrategy(const ChemicalReaction &,
                                                 const BBS &bbs) {
  // If we fail here, someone has a ridiculous amount of memory
  PRECONDITION(
      m_numPermutations != EnumerationStrategyBase::EnumerationOverflow,
      "Cannot represent all permutations for the even sampler");

  boost::uint64_t npos = bbs.size();
  used_count.resize(npos);
  std::fill(used_count.begin(), used_count.end(), 0);

  var_used.resize(npos);
  for (boost::uint64_t i = 0; i < npos; ++i) {
    var_used[i].resize(m_permutationSizes[i]);
    std::fill(var_used[i].begin(), var_used[i].end(), 0);
  }

  boost::uint64_t nmonomers = 0;
  for (boost::uint64_t i = 0; i < bbs.size(); ++i) {
    nmonomers += m_permutationSizes[i];
  }

  pair_used.resize(nmonomers);
  for (boost::uint64_t i = 0; i < nmonomers; ++i) {
    pair_used[i].resize(nmonomers);
    std::fill(pair_used[i].begin(), pair_used[i].end(), 0);
  }

  pair_counts.resize(npos);
  for (boost::uint64_t i = 0; i < npos; i++) {
    pair_counts[i].resize(npos);
    std::fill(pair_counts[i].begin(), pair_counts[i].end(), 0);
  }

  /* Initialize random number generator */
  /* Find modulus */
  for (M = 1; M < rdcast<boost::uint64_t>(m_numPermutations); M = 2 * M) {
    ;
  }
  /* Set factor */
  a = 5;
  b = 7;

  // control of random number and heuristics
  seed = 0;
  m_numPermutationsProcessed = 0;
  nslack = 0;  // increase this to break evenness criteria
  rejected_period = 0;
  rejected_unique = 0;
  rejected_slack_condition = 0;
  rejected_bb_sampling_condition = 0;

  selected.clear();  // clear the selected (unique) set
}

// Try to add the given encoded seed position into
//  the current set of return groups.  This checks to
//  see if the BBS are evenly sampled as pairs.  If
//  they currently are not, reject the selection.
//  This is fairly suboptimal for large collections
//  of building blocks and may take a while to
//  terminate...
bool EvenSamplePairsStrategy::try_add(boost::uint64_t iseed) {
  const RGROUPS &digits = decode(iseed);
  const RGROUPS &rgroups = m_permutationSizes;
  boost::uint64_t islack = 0;
  boost::uint64_t num_rgroups = m_permutationSizes.size();

  for (boost::uint64_t i = 0; i < num_rgroups; ++i) {
    if (var_used[i][digits[i]]) {
      islack += var_used[i][digits[i]];
    }
    if (islack > nslack) {
      // add better heuristic here??
      rejected_slack_condition += 1;
      return false;
    }
  }

  islack = 0;
  boost::uint64_t ioffset = 0;
  // check that building block pairs get evenly sampled
  for (boost::uint64_t i = 0; i < num_rgroups; ++i) {
    boost::uint64_t joffset = 0;
    for (boost::uint64_t j = 0; j < num_rgroups; ++j) {
      if (j == i) {
        continue;
      }
      boost::uint64_t ii = digits[i] + ioffset;
      boost::uint64_t jj = digits[j] + joffset;
      if (pair_used[ii][jj] > 0) {
        auto numer = (double)pair_used[ii][jj];
        double denom = sqrt((double)(rgroups[i]) * (double)(rgroups[j]));
        islack = (int)(numer / denom);
      }
      joffset += rgroups[j];
    }
    ioffset += rgroups[i];
  }

  if (islack > nslack) {
    rejected_bb_sampling_condition += 1;
    return false;
  }

  // keep track of bb usage
  for (boost::uint64_t i = 0; i < num_rgroups; ++i) {
    if (var_used[i][digits[i]] == 0) {
      used_count[i]++;
    }
    var_used[i][digits[i]] += 1;
    if (used_count[i] == rdcast<boost::int64_t>(rgroups[i])) {
      // complete variable scan => initialize
      if (nslack > min_nslack && rgroups[i] > 1) {  // cleared slack on i
        nslack = min_nslack;
      }

      used_count[i] = 0;
      for (boost::uint64_t j = 0; j < rgroups[i]; ++j) {
        var_used[i][j]--;
        if (var_used[i][j] > 0) {
          used_count[i]++;
        }
      }
    }  // end scan
  }

  // keep track of BB Pair usage
  ioffset = 0;
  for (boost::uint64_t i = 0; i < num_rgroups; ioffset += rgroups[i], ++i) {
    boost::uint64_t joffset = 0;
    for (boost::uint64_t j = 0; j < num_rgroups; joffset += rgroups[j], ++j) {
      if (j == i) {
        continue;
      }
      boost::uint64_t ii = digits[i] + ioffset;
      boost::uint64_t jj = digits[j] + joffset;
      if (pair_used[ii][jj] == 0) {
        pair_counts[i][j]++;
      }
      pair_used[ii][jj]++;
      if (pair_counts[i][j] >= rgroups[i] * rgroups[j]) {  // all pairs visited
        if (nslack > min_nslack && (rgroups[i] > 1 || rgroups[j] > 1)) {
          nslack = min_nslack;
        }
        pair_counts[i][j] = 0;
        for (boost::uint64_t ii = 0; ii < rgroups[i]; ++ii) {
          for (boost::uint64_t jj = 0; jj < rgroups[j]; ++jj) {
            pair_used[ioffset + ii][joffset + jj]--;
            if (pair_used[ioffset + ii][joffset + jj] > 0) {
              pair_counts[i][j]++;
            }
          }
        }
      }
    }
  }

  selected.insert(iseed);
  return true;
}

const RGROUPS &EvenSamplePairsStrategy::next() {
  nslack = 0;
  while (m_numPermutationsProcessed <
         rdcast<boost::uint64_t>(m_numPermutations)) {
    for (boost::uint64_t l = 0; l < M; ++l) {
      seed = ((seed * a + b) % M);
      if (seed > rdcast<boost::uint64_t>(m_numPermutations)) {
        rejected_period += 1;
        continue;
      } else if (selected.find(seed) != selected.end()) {
        rejected_unique += 1;
        continue;
      } else if (try_add(seed)) {
        m_numPermutationsProcessed++;
        return decode(seed);
      }
    }

    // loosen heuristic
    nslack += 1;
    min_nslack += 1;
  }

  throw EnumerationStrategyException("Ran out of molecules");
}

std::string EvenSamplePairsStrategy::stats() const {
  std::ostringstream ss;

  boost::uint64_t npos = m_permutationSizes.size();
  const RGROUPS &nvars = m_permutationSizes;
  boost::uint64_t i, l, j, ii, jj, ioffset, joffset;
  ss << "#BEGIN# BBSTAT\n";
  for (i = 0; i < npos; i++) {
    boost::uint64_t maxcount = 0;
    if (nvars[i] == 1) {
      continue;
    }
    for (j = 0; j < nvars[i]; j++) {
      if (maxcount < var_used[i][j]) {
        maxcount = var_used[i][j];
      }
    }
    ss << boost::format("%lu\t%lu\t%6.2f") % (i + 1) % nvars[i] %
              ((double)m_numPermutationsProcessed / nvars[i]);

    for (l = 0; l <= maxcount; l++) {
      boost::uint64_t n = 0;
      for (j = 0; j < nvars[i]; j++) {
        if (var_used[i][j] == l) {
          n++;
        }
        if (n > 0) {
          ss << boost::format("\t%lu|%lu") % l % n;
        }
      }
      ss << std::endl;
    }
  }
  ss << "#END# BBSTAT\n";

  ss << "#BEGIN# PAIRSTAT\n";
  for (i = 0, ioffset = 0; i < npos; ioffset += nvars[i], i++) {
    if (nvars[i] == 1) {
      continue;
    }
    for (j = 0, joffset = 0; j < npos; joffset += nvars[j], j++) {
      boost::uint64_t maxcount = 0;
      if (nvars[j] == 1) {
        continue;
      }
      if (j <= i) {
        continue;
      }
      for (ii = 0; ii < nvars[i]; ii++) {
        for (jj = 0; jj < nvars[j]; jj++) {
          if (maxcount < pair_used[ii + ioffset][jj + joffset]) {
            maxcount = pair_used[ii + ioffset][jj + joffset];
          }
        }
      }
      ss << boost::format("%lu\t%lu\t%lu\t%lu\t%6.2f") % (i + 1) % (j + 1) %
                nvars[i] % nvars[j] %
                ((double)m_numPermutationsProcessed / (nvars[i] * nvars[j]));
      for (l = 0; l <= maxcount; l++) {
        int n = 0;
        for (ii = 0; ii < nvars[i]; ii++) {
          for (jj = 0; jj < nvars[j]; jj++) {
            if (l == pair_used[ii + ioffset][jj + joffset]) {
              n++;
            }
          }
        }
        if (n > 0) {
          ss << boost::format("\t%ld|%d") % l % n;
        }
      }
      ss << boost::format("\n");
    }
  }
  ss << "#END# PAIRSTAT\n";

  ss << "Rejected Period: " << rejected_period << std::endl;
  ss << "Rejected (dupes): " << rejected_unique << std::endl;
  ss << "Rejected Slack Conditions: " << rejected_slack_condition << std::endl;
  ss << "Rejected Pair Sampling: " << rejected_bb_sampling_condition
     << std::endl;
  return ss.str();
}
}  // namespace RDKit
