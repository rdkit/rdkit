//
// Copyright (c) 2016 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/tuple/tuple_comparison.hpp>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

#include "MultiFPBReader.h"
#include <algorithm>

namespace RDKit {

namespace detail {
std::uint8_t *bitsetToBytes(const boost::dynamic_bitset<> &bitset);
}

namespace {
auto tplSorter = [](const MultiFPBReader::ResultTuple &v1,
                    const MultiFPBReader::ResultTuple &v2) {
  if (v1.get<0>() == v2.get<0>()) {
    if (v1.get<2>() == v2.get<2>()) {
      return v1.get<1>() < v2.get<1>();
    } else {
      return v1.get<2>() < v2.get<2>();
    }
  } else {
    return v1.get<0>() > v2.get<0>();
  }
};

auto pairSorter = [](const auto &v1, const auto &v2) {
  if (v1.first == v2.first) {
    return v1.second < v2.second;
  } else {
    return v1.first < v2.first;
  }
};

struct sim_args {
  const std::uint8_t *bv;
  double ca, cb;
  double threshold;
  const std::vector<FPBReader *> &readers;
  std::vector<std::vector<MultiFPBReader::ResultTuple>> *res;
  bool initOnSearch;
};

void tversky_helper(unsigned int threadId, unsigned int numThreads,
                    const sim_args *args) {
  for (unsigned int i = threadId; i < args->readers.size(); i += numThreads) {
    if (args->initOnSearch) {
      args->readers[i]->init();
    }
    std::vector<std::pair<double, unsigned int>> r_res =
        args->readers[i]->getTverskyNeighbors(args->bv, args->ca, args->cb,
                                              args->threshold);
    (*args->res)[i].clear();
    (*args->res)[i].reserve(r_res.size());
    for (std::vector<std::pair<double, unsigned int>>::const_iterator rit =
             r_res.begin();
         rit != r_res.end(); ++rit) {
      (*args->res)[i].push_back(
          MultiFPBReader::ResultTuple(rit->first, rit->second, i));
    }
  }
}
void tani_helper(unsigned int threadId, unsigned int numThreads,
                 const sim_args *args) {
  for (unsigned int i = threadId; i < args->readers.size(); i += numThreads) {
    if (args->initOnSearch) {
      args->readers[i]->init();
    }
    std::vector<std::pair<double, unsigned int>> r_res =
        args->readers[i]->getTanimotoNeighbors(args->bv, args->threshold);
    (*args->res)[i].clear();
    (*args->res)[i].reserve(r_res.size());
    for (std::vector<std::pair<double, unsigned int>>::const_iterator rit =
             r_res.begin();
         rit != r_res.end(); ++rit) {
      (*args->res)[i].push_back(
          MultiFPBReader::ResultTuple(rit->first, rit->second, i));
    }
  }
}

template <typename T>
void generic_nbr_helper(std::vector<MultiFPBReader::ResultTuple> &res, T func,
                        const sim_args &args, unsigned int numThreads) {
  res.clear();
  res.resize(0);
  numThreads = getNumThreadsToUse(numThreads);
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::vector<std::future<void>> tg;
#endif
  if (numThreads == 1) {
    func(0, 1, &args);
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    for (unsigned int tid = 0; tid < numThreads && tid < args.readers.size();
         ++tid) {
      tg.emplace_back(
          std::async(std::launch::async, func, tid, numThreads, &args));
    }
    for (auto &fut : tg) {
      fut.get();
    }
  }
#endif

  for (unsigned int i = 0; i < args.readers.size(); ++i) {
    res.reserve(res.size() + (*args.res).size());
    res.insert(res.end(), (*args.res)[i].begin(), (*args.res)[i].end());
  }
  std::sort(res.begin(), res.end(), tplSorter);
}
void get_tani_nbrs(const std::vector<FPBReader *> &d_readers,
                   const std::uint8_t *bv, double threshold,
                   std::vector<MultiFPBReader::ResultTuple> &res,
                   int numThreads, bool initOnSearch) {
  std::vector<std::vector<MultiFPBReader::ResultTuple>> accum(d_readers.size());
  sim_args args = {bv, 0., 0., threshold, d_readers, &accum, initOnSearch};
  generic_nbr_helper(res, tani_helper, args, numThreads);
}

void get_tversky_nbrs(const std::vector<FPBReader *> &d_readers,
                      const std::uint8_t *bv, double a, double b,
                      double threshold,
                      std::vector<MultiFPBReader::ResultTuple> &res,
                      int numThreads, bool initOnSearch) {
  std::vector<std::vector<MultiFPBReader::ResultTuple>> accum(d_readers.size());
  sim_args args = {bv, a, b, threshold, d_readers, &accum, initOnSearch};
  generic_nbr_helper(res, tversky_helper, args, numThreads);
}

void contain_helper(unsigned int threadId, unsigned int numThreads,
                    const std::uint8_t *bv,
                    const std::vector<FPBReader *> *readers,
                    std::vector<std::vector<unsigned int>> *accum,
                    bool initOnSearch) {
  for (unsigned int i = threadId; i < readers->size(); i += numThreads) {
    if (initOnSearch) {
      (*readers)[i]->init();
    }
    (*accum)[i] = (*readers)[i]->getContainingNeighbors(bv);
  }
}

void get_containing_nbrs(
    const std::vector<FPBReader *> &d_readers, const std::uint8_t *bv,
    std::vector<std::pair<unsigned int, unsigned int>> &res,
    unsigned int numThreads, bool initOnSearch) {
  numThreads = getNumThreadsToUse(numThreads);
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::vector<std::future<void>> tg;
#endif

  std::vector<std::vector<unsigned int>> accum(d_readers.size());
  if (numThreads == 1) {
    contain_helper(0, 1, bv, &d_readers, &accum, initOnSearch);
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    for (unsigned int tid = 0; tid < numThreads && tid < d_readers.size();
         ++tid) {
      tg.emplace_back(std::async(std::launch::async, contain_helper, tid,
                                 numThreads, bv, &d_readers, &accum,
                                 initOnSearch));
    }
    for (auto &fut : tg) {
      fut.get();
    }
  }
#endif

  res.clear();
  for (unsigned int i = 0; i < d_readers.size(); ++i) {
    std::vector<unsigned int> &r_res = accum[i];
    for (auto ri : r_res) {
      res.emplace_back(ri, i);
    }
  }

  std::sort(res.begin(), res.end(), pairSorter);
}

}  // end of anonymous namespace

void MultiFPBReader::init() {
  unsigned int nBits = 0;
  for (auto *rdr : d_readers) {
    rdr->init();
    if (!nBits) {
      nBits = rdr->nBits();
    } else {
      if (rdr->nBits() != nBits) {
        throw ValueErrorException("bit lengths of child readers don't match");
      }
    }
  }
  df_init = true;
};

MultiFPBReader::MultiFPBReader(std::vector<FPBReader *> &readers,
                               bool takeOwnership, bool initOnSearch) {
  df_init = false;
  df_takeOwnership = takeOwnership;
  df_initOnSearch = initOnSearch;
  for (auto *rdr : readers) {
    PRECONDITION(rdr != nullptr, "bad reader");
  }
  d_readers = readers;
}

FPBReader *MultiFPBReader::getReader(unsigned int which) {
  URANGE_CHECK(which, d_readers.size());
  return d_readers[which];
}
unsigned int MultiFPBReader::nBits() const {
  PRECONDITION(d_readers.size(), "no readers");
  PRECONDITION(df_init, "not initialized");
  return d_readers[0]->nBits();
}

std::vector<MultiFPBReader::ResultTuple> MultiFPBReader::getTanimotoNeighbors(
    const std::uint8_t *bv, double threshold, int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");

  std::vector<MultiFPBReader::ResultTuple> res;
  get_tani_nbrs(d_readers, bv, threshold, res, numThreads, df_initOnSearch);
  return res;
}

std::vector<MultiFPBReader::ResultTuple> MultiFPBReader::getTanimotoNeighbors(
    const ExplicitBitVect &ebv, double threshold, int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");
  std::vector<MultiFPBReader::ResultTuple> res;
  std::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  get_tani_nbrs(d_readers, bv, threshold, res, numThreads, df_initOnSearch);
  delete[] bv;
  return res;
}

std::vector<MultiFPBReader::ResultTuple> MultiFPBReader::getTverskyNeighbors(
    const std::uint8_t *bv, double ca, double cb, double threshold,
    int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");
  std::vector<MultiFPBReader::ResultTuple> res;
  get_tversky_nbrs(d_readers, bv, ca, cb, threshold, res, numThreads,
                   df_initOnSearch);
  return res;
}

std::vector<MultiFPBReader::ResultTuple> MultiFPBReader::getTverskyNeighbors(
    const ExplicitBitVect &ebv, double ca, double cb, double threshold,
    int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");
  std::vector<MultiFPBReader::ResultTuple> res;
  std::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  get_tversky_nbrs(d_readers, bv, ca, cb, threshold, res, numThreads,
                   df_initOnSearch);
  delete[] bv;
  return res;
}

std::vector<std::pair<unsigned int, unsigned int>>
MultiFPBReader::getContainingNeighbors(const std::uint8_t *bv,
                                       int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");
  std::vector<std::pair<unsigned int, unsigned int>> res;
  get_containing_nbrs(d_readers, bv, res, numThreads, df_initOnSearch);
  return res;
}

std::vector<std::pair<unsigned int, unsigned int>>
MultiFPBReader::getContainingNeighbors(const ExplicitBitVect &ebv,
                                       int numThreads) const {
  PRECONDITION(df_init || df_initOnSearch, "not initialized");
  std::vector<std::pair<unsigned int, unsigned int>> res;
  std::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  get_containing_nbrs(d_readers, bv, res, numThreads, df_initOnSearch);
  delete[] bv;
  return res;
}

}  // namespace RDKit
