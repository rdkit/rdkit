//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file Charges.cpp
/// \brief charge normalisation and the threaded batch entry points

#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>
#include <numeric>

#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>

#include "DASHTreeImpl.h"

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <exception>
#include <mutex>
#include <thread>
#endif

namespace RDKit {
namespace DASH {

using namespace detail;

namespace {

//! \brief Neumaier (compensated) summation
/*!
  The normalisation divides by these sums, so the last bits of the result depend
  on how they were accumulated. Compensated summation tracks the rounding error
  a naive loop throws away, which both tightens the result and makes it agree
  exactly with the reference python implementation -- CPython's \c sum() has
  used the same algorithm for floats since 3.12.
*/
double compensatedSum(const std::vector<double> &values) {
  double sum = 0.0;
  double compensation = 0.0;
  for (const double value : values) {
    const double running = sum + value;
    if (std::abs(sum) >= std::abs(value)) {
      compensation += (sum - running) + value;
    } else {
      compensation += (value - running) + sum;
    }
    sum = running;
  }
  return sum + compensation;
}

//! \brief spreads the difference to the molecule's formal charge over the atoms
/*!
  The tree stores each atom's charge independently, so the raw values do not
  quite add up to the molecule's charge. SYMMETRIC hands every atom the same
  share of the difference; STD_WEIGHTED gives more of it to the atoms the tree
  is least certain about, which is the published default.
*/
void normalizeCharges(const ROMol &mol, const std::vector<double> &rawValues,
                      const std::vector<double> &stds,
                      ChargeNormalization normalization,
                      std::vector<double> &res) {
  const std::size_t nAtoms = rawValues.size();
  if (normalization == ChargeNormalization::NONE || nAtoms == 0) {
    res = rawValues;
    return;
  }
  res.resize(nAtoms);

  const double treeTotal = compensatedSum(rawValues);
  int molTotal = 0;
  for (const auto atom : mol.atoms()) {
    molTotal += atom->getFormalCharge();
  }
  const double deficit = molTotal - treeTotal;

  if (normalization == ChargeNormalization::SYMMETRIC) {
    const double share = deficit / static_cast<double>(nAtoms);
    for (std::size_t i = 0; i < nAtoms; ++i) {
      res[i] = rawValues[i] + share;
    }
    return;
  }

  const double stdTotal = compensatedSum(stds);
  for (std::size_t i = 0; i < nAtoms; ++i) {
    res[i] = rawValues[i] + deficit * (stds[i] / stdTotal);
  }
}

//! the per-molecule charge assignment, shared by every entry point
void assignCharges(MolMatcher &matcher, const ROMol &mol,
                   const PropertyColumn &valueColumn,
                   const PropertyColumn &stdColumn,
                   const ChargeOptions &options, std::vector<double> &res,
                   std::vector<double> &rawValues, std::vector<double> &stds,
                   std::vector<unsigned int> &matchDepths) {
  const unsigned int nAtoms = mol.getNumAtoms();
  matcher.setMolecule(mol);
  rawValues.resize(nAtoms);
  stds.resize(nAtoms);
  matchDepths.resize(nAtoms);

  for (unsigned int i = 0; i < nAtoms; ++i) {
    matcher.match(i, options.params);
    rawValues[i] = matcher.pathValue(valueColumn);
    if (std::isnan(rawValues[i])) {
      // No node on this atom's path carries a value, so there is nothing to
      // normalise and a silent NaN would spread to every other atom. This is
      // the case where the reference python implementation raises too, though
      // it does so from inside pandas.
      const Atom *atom = mol.getAtomWithIdx(i);
      throw ValueErrorException(
          "the DASH tree holds no '" + options.valueProperty + "' for atom " +
          std::to_string(i) + " (" + atom->getSymbol() + ", degree " +
          std::to_string(atom->getDegree()) + ", charge " +
          std::to_string(atom->getFormalCharge()) +
          "): its atom type is in the tree but carries no data");
    }
    const double deviation = matcher.pathValue(stdColumn);
    // a zero deviation would give the atom no share of the normalisation, and
    // a missing one (NaN) would poison the whole molecule
    stds[i] = (deviation > 0.0) ? deviation : options.defaultStdValue;
    matchDepths[i] = static_cast<unsigned int>(matcher.path().size());
  }
  normalizeCharges(mol, rawValues, stds, options.normalization, res);
}

// ---------------------------------------------------------------------------
//  Work is handed out one molecule at a time through an atomic counter rather
//  than sliced up in advance: real batches vary in molecule size by an order of
//  magnitude, and a fixed slicing leaves threads idle at the end of the run.
//  Each thread keeps one matcher for the whole batch, so the per-molecule
//  tables are the only thing reallocated, and writes straight into its own slot
//  of the result, so nothing has to be merged afterwards.
// ---------------------------------------------------------------------------

//! \brief runs \c work(threadIdx, i) over [0, \p n) on \p nThreads threads
/*!
  Stops early on a Ctrl-C while a ControlCHandler is installed, as the python
  wrappers do around a batch, raising KeyboardInterrupt once they have the GIL
  back.
*/
template <typename Work>
void forEachIndex(std::size_t n, unsigned int nThreads, Work work) {
  if (nThreads <= 1 || n <= 1) {
    for (std::size_t i = 0; i < n; ++i) {
      if (ControlCHandler::getGotSignal()) {
        return;
      }
      work(0u, i);
    }
    return;
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::atomic<std::size_t> next(0);
  std::vector<std::thread> threads;
  threads.reserve(nThreads);
  std::exception_ptr failure;
  std::mutex failureMutex;
  for (unsigned int t = 0; t < nThreads; ++t) {
    threads.emplace_back([&, t]() {
      try {
        for (;;) {
          const std::size_t i = next.fetch_add(1);
          if (i >= n || ControlCHandler::getGotSignal()) {
            break;
          }
          work(t, i);
        }
      } catch (...) {
        const std::lock_guard<std::mutex> lock(failureMutex);
        if (!failure) {
          failure = std::current_exception();
        }
      }
    });
  }
  for (auto &thread : threads) {
    thread.join();
  }
  if (failure) {
    std::rethrow_exception(failure);
  }
#else
  for (std::size_t i = 0; i < n; ++i) {
    work(0u, i);
  }
#endif
}

//! how many threads to actually spawn: never more than there is work for
unsigned int resolveThreadCount(int numThreads, std::size_t nWork) {
  const unsigned int wanted = getNumThreadsToUse(numThreads);
  return static_cast<unsigned int>(
      std::min<std::size_t>(wanted, std::max<std::size_t>(nWork, 1)));
}

//! one reusable matcher per thread
std::vector<std::unique_ptr<MolMatcher>> makeMatchers(
    const DASHTree::Impl &tree, unsigned int nThreads) {
  std::vector<std::unique_ptr<MolMatcher>> matchers(nThreads);
  for (auto &matcher : matchers) {
    matcher = std::make_unique<MolMatcher>(tree);
  }
  return matchers;
}

}  // namespace

void DASHTree::getPartialCharges(const ROMol &mol, std::vector<double> &res,
                                 std::vector<double> &rawValues,
                                 std::vector<double> &stds,
                                 std::vector<unsigned int> &matchDepths,
                                 const ChargeOptions &options) const {
  const PropertyColumn &valueColumn = d_impl->property(options.valueProperty);
  const PropertyColumn &stdColumn = d_impl->property(options.stdProperty);
  MolMatcher matcher(*d_impl);
  assignCharges(matcher, mol, valueColumn, stdColumn, options, res, rawValues,
                stds, matchDepths);
}

void DASHTree::getPartialCharges(const ROMol &mol, std::vector<double> &res,
                                 const ChargeOptions &options) const {
  std::vector<double> rawValues, stds;
  std::vector<unsigned int> matchDepths;
  getPartialCharges(mol, res, rawValues, stds, matchDepths, options);
}

void DASHTree::getMolPropertyBatch(const std::vector<const ROMol *> &mols,
                                   const std::string &property,
                                   std::vector<std::vector<double>> &res,
                                   const DASHParams &params,
                                   int numThreads) const {
  const PropertyColumn &column = d_impl->property(property);
  res.clear();
  res.resize(mols.size());

  const unsigned int nThreads = resolveThreadCount(numThreads, mols.size());
  auto matchers = makeMatchers(*d_impl, nThreads);
  forEachIndex(mols.size(), nThreads, [&](unsigned int t, std::size_t i) {
    if (!mols[i]) {
      return;
    }
    MolMatcher &matcher = *matchers[t];
    const ROMol &mol = *mols[i];
    matcher.setMolecule(mol);
    res[i].resize(mol.getNumAtoms());
    for (unsigned int a = 0; a < mol.getNumAtoms(); ++a) {
      matcher.match(a, params);
      res[i][a] = matcher.pathValue(column);
    }
  });
}

void DASHTree::getPartialChargesBatch(const std::vector<const ROMol *> &mols,
                                      std::vector<std::vector<double>> &res,
                                      const ChargeOptions &options,
                                      int numThreads) const {
  const PropertyColumn &valueColumn = d_impl->property(options.valueProperty);
  const PropertyColumn &stdColumn = d_impl->property(options.stdProperty);
  res.clear();
  res.resize(mols.size());

  const unsigned int nThreads = resolveThreadCount(numThreads, mols.size());
  auto matchers = makeMatchers(*d_impl, nThreads);
  std::vector<std::vector<double>> rawScratch(nThreads), stdScratch(nThreads);
  std::vector<std::vector<unsigned int>> depthScratch(nThreads);
  forEachIndex(mols.size(), nThreads, [&](unsigned int t, std::size_t i) {
    if (!mols[i]) {
      return;
    }
    assignCharges(*matchers[t], *mols[i], valueColumn, stdColumn, options,
                  res[i], rawScratch[t], stdScratch[t], depthScratch[t]);
  });
}

}  // namespace DASH
}  // namespace RDKit
