// $Id$
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PeriodicTable.h"
#include <string>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <sstream>
#include <locale>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace RDKit {

class std::unique_ptr<PeriodicTable> PeriodicTable::ds_instance = nullptr;

PeriodicTable::PeriodicTable() {
  // it is assumed that the atomic atomData string constains atoms
  // in sequence and no atoms are missing in between
  byanum.clear();
  byname.clear();

  boost::char_separator<char> eolSep("\n");
  tokenizer tokens(periodicTableAtomData, eolSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    if (*token != " ") {
      atomicData adata(*token);
      std::string enam = adata.Symbol();
      byname[enam] = adata.AtomicNum();
      // there are, for backwards compatibility reasons, some duplicate rows for
      // atomic numbers in the atomic_data data structure. It's ok to have
      // multiple symbols map to the same atomic number (above), but we need to
      // be sure that we only store one entry per atomic number.
      // Note that this only works because the first atom in the adata list is
      // the dummy atom (atomic number 0). This was #2784
      if (rdcast<size_t>(adata.AtomicNum()) == byanum.size()) {
        byanum.push_back(adata);
      }
    }
  }

  unsigned int lidx = 0;
  std::istringstream istr;
  istr.imbue(std::locale("C"));
  while (isotopesAtomData[lidx] != "" && isotopesAtomData[lidx] != "EOS") {
    tokenizer lines(isotopesAtomData[lidx++], eolSep);
    boost::char_separator<char> spaceSep(" \t");
    for (tokenizer::iterator line = lines.begin(); line != lines.end();
         ++line) {
      if (*line != " ") {
        tokenizer tokens(*line, spaceSep);
        tokenizer::iterator token = tokens.begin();
        int anum;
        istr.clear();
        istr.str(*token);
        istr >> anum;
        atomicData &adata = byanum[anum];
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        unsigned int isotope;
        istr.clear();
        istr.str(*token);
        istr >> isotope;
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        double mass;
        istr.clear();
        istr.str(*token);
        istr >> mass;
        ++token;
        if (token == tokens.end()) {
          continue;
        }
        double abundance;
        istr.clear();
        istr.str(*token);
        istr >> abundance;
        adata.d_isotopeInfoMap[isotope] = std::make_pair(mass, abundance);
      }
    }
  }
}

void PeriodicTable::initInstance() {
  ds_instance = std::unique_ptr<PeriodicTable>(new PeriodicTable());
}

PeriodicTable *PeriodicTable::getTable() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static std::once_flag pt_init_once;
  std::call_once(pt_init_once, initInstance);
#else
  if (!ds_instance) {
    initInstance();
  }
#endif
  return ds_instance.get();
}

}  // namespace RDKit
