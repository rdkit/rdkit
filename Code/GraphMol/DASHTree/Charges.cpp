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

#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>

#include "DASHTreeImpl.h"

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
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

//! \brief spreads the difference to a molecule-wide target over the atoms
/*!
  The tree stores each atom's value independently, so the raw values do not
  quite add up to what the molecule as a whole should carry -- its formal
  charge, for partial charges. SYMMETRIC hands every atom the same share of the
  difference; STD_WEIGHTED gives more of it to the atoms the tree is least
  certain about, which is the published default for charges.
*/
void normalizeValues(double target, const std::vector<double> &rawValues,
                     const std::vector<double> &stds,
                     Normalization normalization, std::vector<double> &res) {
  const std::size_t nAtoms = rawValues.size();
  if (normalization == Normalization::NONE || nAtoms == 0) {
    res = rawValues;
    return;
  }
  res.resize(nAtoms);

  const double deficit = target - compensatedSum(rawValues);

  if (normalization == Normalization::SYMMETRIC) {
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

void checkOptions(const NormalizationOptions &options) {
  if (options.normalization == Normalization::STD_WEIGHTED &&
      !(options.defaultStdValue > 0.0)) {
    throw ValueErrorException(
        "defaultStdValue must be positive: it stands in for a missing "
        "deviation, and STD_WEIGHTED divides by the sum of the deviations");
  }
}

//! \brief the deviation column \p options asks for
/*!
  STD_WEIGHTED cannot do without it; the other normalisations report it when
  it was resolved and fall back to \c defaultStdValue when it does not, so a
  column without a deviation column of its own can still be normalised.
*/
const PropertyColumn *stdColumnFor(const DASHTree::Impl &impl,
                                   const NormalizationOptions &options) {
  if (options.normalization == Normalization::STD_WEIGHTED ||
      impl.propertyIndex(options.stdProperty) >= 0) {
    return &impl.property(options.stdProperty);
  }
  return nullptr;
}

//! the per-molecule normalised assignment, shared by every entry point
void assignNormalized(MolMatcher &matcher, const ROMol &mol, double target,
                      const PropertyColumn &valueColumn,
                      const PropertyColumn *stdColumn,
                      const NormalizationOptions &options,
                      std::vector<double> &res, std::vector<double> &rawValues,
                      std::vector<double> &stds,
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
    // a zero deviation would give the atom no share of the normalisation, and
    // a missing one (NaN) would poison the whole molecule
    const double deviation =
        stdColumn ? matcher.pathValue(*stdColumn) : options.defaultStdValue;
    stds[i] = (deviation > 0.0) ? deviation : options.defaultStdValue;
    matchDepths[i] = static_cast<unsigned int>(matcher.path().size());
  }
  normalizeValues(target, rawValues, stds, options.normalization, res);
}

//! \brief runs \c work(matcher, mol, i) over the molecules of \p mols,
//! skipping null entries
/*!
  Molecules are handed out one at a time through an atomic counter, so a batch
  of unequal molecules still balances across the threads. Each thread keeps one
  matcher for the whole batch.
*/
template <typename Work>
void forEachMol(const DASHTree::Impl &tree,
                const std::vector<const ROMol *> &mols, int numThreads,
                Work work) {
  std::atomic<std::size_t> next(0);
  auto run = [&]() {
    MolMatcher matcher(tree);
    try {
      for (std::size_t i = next++; i < mols.size(); i = next++) {
        if (mols[i]) {
          work(matcher, *mols[i], i);
        }
      }
    } catch (...) {
      next = mols.size();  // the other threads stop at their current molecule
      throw;
    }
  };
  const auto nThreads =
      std::min<std::size_t>(mols.size(), getNumThreadsToUse(numThreads));
#ifdef RDK_BUILD_THREADSAFE_SSS
  if (nThreads > 1) {
    std::vector<std::future<void>> tasks;
    tasks.reserve(nThreads);
    for (std::size_t t = 0; t < nThreads; ++t) {
      tasks.emplace_back(std::async(std::launch::async, run));
    }
    for (auto &task : tasks) {
      task.get();
    }
    return;
  }
#else
  RDUNUSED_PARAM(nThreads);
#endif
  run();
}

}  // namespace

void DASHTree::getNormalizedMolProperty(
    const ROMol &mol, double target, std::vector<double> &res,
    std::vector<double> &rawValues, std::vector<double> &stds,
    std::vector<unsigned int> &matchDepths,
    const NormalizationOptions &options) const {
  checkOptions(options);
  const PropertyColumn &valueColumn = d_impl->property(options.valueProperty);
  const PropertyColumn *stdColumn = stdColumnFor(*d_impl, options);
  MolMatcher matcher(*d_impl);
  assignNormalized(matcher, mol, target, valueColumn, stdColumn, options, res,
                   rawValues, stds, matchDepths);
}

void DASHTree::getNormalizedMolProperty(
    const ROMol &mol, double target, std::vector<double> &res,
    const NormalizationOptions &options) const {
  std::vector<double> rawValues, stds;
  std::vector<unsigned int> matchDepths;
  getNormalizedMolProperty(mol, target, res, rawValues, stds, matchDepths,
                           options);
}

void DASHTree::getPartialCharges(const ROMol &mol, std::vector<double> &res,
                                 std::vector<double> &rawValues,
                                 std::vector<double> &stds,
                                 std::vector<unsigned int> &matchDepths,
                                 const NormalizationOptions &options) const {
  getNormalizedMolProperty(mol, MolOps::getFormalCharge(mol), res, rawValues,
                           stds, matchDepths, options);
}

void DASHTree::getPartialCharges(const ROMol &mol, std::vector<double> &res,
                                 const NormalizationOptions &options) const {
  getNormalizedMolProperty(mol, MolOps::getFormalCharge(mol), res, options);
}

void DASHTree::getMolPropertyBatch(const std::vector<const ROMol *> &mols,
                                   const std::string &property,
                                   std::vector<std::vector<double>> &res,
                                   const DASHParams &params,
                                   int numThreads) const {
  const PropertyColumn &column = d_impl->property(property);
  res.clear();
  res.resize(mols.size());
  forEachMol(*d_impl, mols, numThreads,
             [&](MolMatcher &matcher, const ROMol &mol, std::size_t i) {
               matcher.setMolecule(mol);
               res[i].resize(mol.getNumAtoms());
               for (unsigned int a = 0; a < mol.getNumAtoms(); ++a) {
                 matcher.match(a, params);
                 res[i][a] = matcher.pathValue(column);
               }
             });
}

void DASHTree::getNormalizedMolPropertyBatch(
    const std::vector<const ROMol *> &mols, const std::vector<double> &targets,
    std::vector<std::vector<double>> &res, const NormalizationOptions &options,
    int numThreads) const {
  if (targets.size() != mols.size()) {
    throw ValueErrorException(
        "one target per molecule: " + std::to_string(targets.size()) +
        " targets for " + std::to_string(mols.size()) + " molecules");
  }
  checkOptions(options);
  const PropertyColumn &valueColumn = d_impl->property(options.valueProperty);
  const PropertyColumn *stdColumn = stdColumnFor(*d_impl, options);
  res.clear();
  res.resize(mols.size());
  forEachMol(*d_impl, mols, numThreads,
             [&](MolMatcher &matcher, const ROMol &mol, std::size_t i) {
               std::vector<double> rawValues, stds;
               std::vector<unsigned int> matchDepths;
               assignNormalized(matcher, mol, targets[i], valueColumn,
                                stdColumn, options, res[i], rawValues, stds,
                                matchDepths);
             });
}

void DASHTree::getPartialChargesBatch(const std::vector<const ROMol *> &mols,
                                      std::vector<std::vector<double>> &res,
                                      const NormalizationOptions &options,
                                      int numThreads) const {
  std::vector<double> targets(mols.size(), 0.0);
  for (std::size_t i = 0; i < mols.size(); ++i) {
    if (mols[i]) {
      targets[i] = MolOps::getFormalCharge(*mols[i]);
    }
  }
  getNormalizedMolPropertyBatch(mols, targets, res, options, numThreads);
}

}  // namespace DASH
}  // namespace RDKit
