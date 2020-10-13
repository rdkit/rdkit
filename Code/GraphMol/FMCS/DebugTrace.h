//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <ctime>
#include <iostream>
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#include <Windows.h>  // for Winmm.lib timeGetTime()
#ifdef _DEBUG         // check memory leaks
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#ifndef new
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif
#endif
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif
#endif

// SELECT ALGORITHM OPTIONS by comment some lines to exclude additional or
// experimental optimisations:

#define SEED_GROW_DEEP  // fast and works much times faster (but it can depend
                        // on molecules)
//#define EXCLUDE_WRONG_COMPOSITION   // fast but with a little effect, because
// amount of external bonds usually is very small.
// Exclude mismatched bonds combinations during seed growing (2^N-1 stage)

#define FAST_SUBSTRUCT_CACHE  // based on a hash of Morgan code
#define DUP_SUBSTRUCT_CACHE   // based on list of query atoms and bonds. For
                              // rings where seeds growing in both directions
                              // throw the same ring.

#define FAST_INCREMENTAL_MATCH  // fast and some time very useful. request
                                // PRECOMPUTED_TABLES_MATCH
// previous match result based match checking without finding new matched
// substructure location in the target

#define VERBOSE_STATISTICS_ON

#ifdef _MSC_VER
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL

struct timezone {
  int tz_minuteswest;  // minutes W of Greenwich
  int tz_dsttime;      // type of dst correction
};

static inline int gettimeofday(struct timeval *tv, struct timezone *tz) {
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;

  if (nullptr != tv) {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    // converting file time to unix epoch
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
    tmpres /= 10;  // convert into microseconds
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (nullptr != tz) {
    if (!tzflag) {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
  return 0;
}
#endif

static inline unsigned long long nanoClock(
    void) {  // actually returns microseconds
  struct timeval t;
  gettimeofday(&t, (struct timezone *)nullptr);
  return t.tv_usec + t.tv_sec * 1000000ULL;
}

namespace RDKit {
namespace FMCS {

#ifdef VERBOSE_STATISTICS_ON

// compute statistics of really very very fast calls.
// It a bit decrease overal performance, but might be interested for
// investigation purpose (only)
//#define VERBOSE_STATISTICS_FASTCALLS_ON

struct ExecStatistics {
  unsigned TotalSteps{0}, MCSFoundStep{0};
  unsigned long long MCSFoundTime;
  unsigned InitialSeed{0}, MismatchedInitialSeed{0};
  unsigned Seed{0}, RemainingSizeRejected{0};
  unsigned SeedCheck{0}, SingleBondExcluded{0};
  unsigned MatchCall{0}, MatchCallTrue{0};
  unsigned FastMatchCall{0}, FastMatchCallTrue{0}, SlowMatchCallTrue{0};
  unsigned ExactMatchCall{0}, ExactMatchCallTrue{0};  // hash cache
  unsigned FindHashInCache{0}, HashKeyFoundInCache{0};
  unsigned AtomCompareCalls{0}, BondCompareCalls{0};  // long long
  unsigned AtomFunctorCalls{0}, BondFunctorCalls{0};  // long long
  unsigned WrongCompositionRejected{0}, WrongCompositionDetected{0};
  unsigned DupCacheFound{0}, DupCacheFoundMatch{0};

  ExecStatistics() : MCSFoundTime(nanoClock()) {}
};
#endif
}  // namespace FMCS
}  // namespace RDKit
