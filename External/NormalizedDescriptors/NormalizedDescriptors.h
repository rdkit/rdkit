//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Histogram normalization of molecular descriptors, ported from the
//  RDKit2DHistogramNormalized descriptors of descriptastorus
//  (https://github.com/bp-kelley/descriptastorus).
//
#include <RDGeneral/export.h>
#ifndef RDKIT_NORMALIZEDDESCRIPTORS_H
#define RDKIT_NORMALIZEDDESCRIPTORS_H

#include <GraphMol/Descriptors/Property.h>

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
namespace NormalizedDescriptors {

//! A cumulative histogram of a descriptor over a reference set of molecules.
/*!
  The table is a list of bins, each with a left edge and the fraction of the
  reference set in that bin or below it, sorted by edge. A descriptor value v
  normalizes to the fraction of the first bin whose edge is >= v, or to 1.0
  when v is greater than every edge. There is no interpolation and no clipping.
  NaN normalizes to 0.0.

  This reproduces descriptastorus' RDKit2DHistogramNormalized, which looks
  the value up with python's bisect:
    bins[bisect(bins, (v,))][1]
  (see tools/make_normalized_histograms.py for how the tables are made).
*/
class RDKIT_NORMALIZEDDESCRIPTORS_EXPORT CumulativeHistogram {
 public:
  CumulativeHistogram() = default;
  //! \c edges must be non-empty, sorted ascending, and the same size as
  //! \c fractions
  CumulativeHistogram(std::vector<double> edges, std::vector<double> fractions);

  //! returns the normalized value for \c value
  double normalize(double value) const;
  double operator()(double value) const { return normalize(value); }

  //! the left edge of each bin
  const std::vector<double> &getEdges() const { return d_edges; }
  //! the cumulative fraction of the reference set for each bin
  const std::vector<double> &getFractions() const { return d_fractions; }

 private:
  std::vector<double> d_edges;
  std::vector<double> d_fractions;
};

//! A named collection of histogram tables, one per descriptor.
/*!
  The text format read and written by this class (the one produced by
  tools/make_normalized_histograms.py) has, for each descriptor, a header line
  followed by nbins lines of bins sorted by edge:

    histogram <name> <nbins>
    <edge> <cumulative fraction>
    ...

  Fields are whitespace separated. Blank lines and lines starting with '#'
  are ignored.
*/
class RDKIT_NORMALIZEDDESCRIPTORS_EXPORT HistogramSet {
 public:
  HistogramSet() = default;

  //! adds (or replaces) the table for descriptor \c name
  /*!
    A new table gets the next index (see getTableIndex()); replacing a table
    keeps its index.
  */
  void addTable(const std::string &name, CumulativeHistogram table);
  bool hasTable(const std::string &name) const;
  //! throws a KeyErrorException if \c name has no table
  const CumulativeHistogram &getTable(const std::string &name) const;
  //! returns the table with index \c idx, throws an IndexErrorException if
  //! \c idx is out of range
  const CumulativeHistogram &getTable(size_t idx) const;
  //! returns the index of the table for \c name, or -1 if there is none
  int getTableIndex(const std::string &name) const;
  //! returns the descriptor names, in index order
  const std::vector<std::string> &getNames() const { return d_names; }
  size_t size() const { return d_tables.size(); }

  //! normalizes \c value with the table for \c name.
  /*!
    Following descriptastorus, a descriptor without a table normalizes to 0.0.
  */
  double normalize(const std::string &name, double value) const;
  //! normalizes \c value with the table with index \c idx, skipping the name
  //! lookup. Throws an IndexErrorException if \c idx is out of range.
  double normalize(size_t idx, double value) const;

  //! reads tables from \c inStream, adding them to this set
  /*! throws a ValueErrorException on malformed input */
  void loadFromStream(std::istream &inStream);
  //! reads tables from the file \c fileName, adding them to this set
  void loadFromFile(const std::string &fileName);
  //! writes the tables in this set to \c outStream, in index order
  void writeToStream(std::ostream &outStream) const;

 private:
  std::vector<CumulativeHistogram> d_tables;
  std::vector<std::string> d_names;
  std::map<std::string, size_t> d_index;
};

//! returns the path of the histogram tables built from the reference
//! distributions of descriptastorus' RDKit2DHistogramNormalized:
//!   $RDBASE/External/NormalizedDescriptors/data/normalized_descriptor_histograms.txt
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::string getDefaultTablePath();

//! returns the default histogram tables, loading them from
//! getDefaultTablePath() on first use
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT const HistogramSet &getDefaultTables();

//! Computes descriptors with the RDKit property registry
//! (Descriptors::Properties) and normalizes them with histogram tables.
/*!
  By default it computes the descriptors of getNormalizedDescriptorNames(),
  the descriptors in rdkit.Chem.Descriptors._descList in the same order, and
  normalizes them with getDefaultTables().

  Note that the names here in getPropertyNames are mapped to the ones
  in _descList which are NOT the same as C++ registeres the properties.
  
  Descriptors are named as in rdkit.Chem.Descriptors._descList. Where the
  property registry already had the descriptor under another name (e.g.
  "exactmw" for "ExactMolWt") that property is used, and its table is looked
  up by the _descList name.

  computeProperties() returns normalized values. A descriptor that fails to
  compute, or that has no table, gets the failure value, which is 0.0 by
  default as in descriptastorus (see setFailureValue()).
*/
class RDKIT_NORMALIZEDDESCRIPTORS_EXPORT NormalizedProperties
    : public Descriptors::Properties {
 public:
  //! the descriptors of getNormalizedDescriptorNames(), with the default
  //! tables
  NormalizedProperties();
  //! the descriptors of getNormalizedDescriptorNames(), with \c tables
  explicit NormalizedProperties(const HistogramSet &tables);
  //! the descriptors \c names, with \c tables. Throws a KeyErrorException
  //! for a name that is neither a _descList name nor a registered property.
  NormalizedProperties(const std::vector<std::string> &names,
                       const HistogramSet &tables);
  //! as above, sharing \c tables instead of copying them
  NormalizedProperties(const std::vector<std::string> &names,
                       std::shared_ptr<const HistogramSet> tables);

  //! returns the normalized descriptor values for \c mol; with
  //! \c annotate, also sets them as properties on \c mol, named by
  //! getDescriptorNames()
  std::vector<double> computeProperties(const ROMol &mol,
                                        bool annotate = false) const override;
  //! returns the raw (unnormalized) descriptor values for \c mol; a
  //! descriptor that fails to compute is NaN
  std::vector<double> computeRawProperties(const ROMol &mol) const;

  //! the descriptor names (_descList names where there is one), in the order
  //! values are returned
  std::vector<std::string> getPropertyNames() const override {
    return d_descriptorNames;
  }
  const HistogramSet &getTables() const { return *d_tables; }

 private:
  std::vector<std::string> d_descriptorNames;
  std::shared_ptr<const HistogramSet> d_tables;
  // table index for each descriptor, -1 if it has none
  std::vector<int> d_tableIndex;
};

//! returns the names of the descriptors calculated by
//! calcNormalizedDescriptors(), in the order they are returned.
/*!
  These are the descriptors in rdkit.Chem.Descriptors._descList, in the same
  order.
*/
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT const std::vector<std::string> &
getNormalizedDescriptorNames();

//! returns the name of the registered property (Descriptors::Properties)
//! that computes the descriptor \c name, e.g. "exactmw" for "ExactMolWt";
//! names without an alias are returned unchanged
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::string getPropertyName(
    const std::string &name);

//! calculates the raw (unnormalized) values of the descriptors named by
//! getNormalizedDescriptorNames(); a descriptor that cannot be calculated
//! for \c mol is NaN
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::vector<double> calcDescriptorValues(
    const ROMol &mol);

//! calculates the descriptors named by getNormalizedDescriptorNames() and
//! normalizes them with \c tables
/*!
  Following descriptastorus, a descriptor that cannot be
  calculated, or that has no table, is 0.0.
*/
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::vector<double>
calcNormalizedDescriptors(const ROMol &mol, const HistogramSet &tables);

//! calculates the normalized descriptors using getDefaultTables()
RDKIT_NORMALIZEDDESCRIPTORS_EXPORT std::vector<double>
calcNormalizedDescriptors(const ROMol &mol);

}  // namespace NormalizedDescriptors
}  // namespace RDKit

#endif
