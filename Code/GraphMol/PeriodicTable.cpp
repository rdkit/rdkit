//
//  Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PeriodicTable.h"

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace RDKit {

class std::unique_ptr<PeriodicTable> PeriodicTable::ds_instance = nullptr;

PeriodicTable::PeriodicTable() {
  // it is assumed that the atomicDataVals vector contains atoms
  // in sequence and no atoms are missing in between
  byanum.clear();
  byname.clear();
  for (const auto &adata : atomicDataVals) {
    byname[adata.Symbol()] = adata.AtomicNum();
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
  for (const auto &idata : isotopeDataVals) {
    atomicData &adata = byanum[idata.atomicNumber];
    adata.d_isotopeInfoMap[idata.isotope] =
        std::make_pair(idata.mass, idata.abundance);
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
