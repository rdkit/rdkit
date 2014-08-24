//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <iostream>
#ifdef WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <Windows.h> // for Winmm.lib timeGetTime()
#ifdef _DEBUG   // check memory leaks
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#ifndef new
#define new new( _NORMAL_BLOCK, __FILE__, __LINE__)
#endif
#endif
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

// SELECT ALGORITHM OPTIONS by comment some lines to exclude additional or experimental optimisations:

#define SEED_GROW_DEEP              // fast and works much times faster (but it can depend on molecules)
//#define EXCLUDE_WRONG_COMPOSITION   // fast but with a little effect, because amount of external bonds usually is very small.
// Exclude mismatched bonds combinations during seed growing (2^N-1 stage)

#define FAST_SUBSTRUCT_CACHE        // based on a hash of Morgan code
#define DUP_SUBSTRUCT_CACHE         // based on list of query atoms and bonds. For rings where seeds growing in both directions throw the same ring.

#define FAST_INCREMENTAL_MATCH      // fast and some time very usefull. request PRECOMPUTED_TABLES_MATCH
// previous match result based match checking without finding new matched substructure location in the target

#define VERBOSE_STATISTICS_ON


#ifdef WIN32
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL

struct timezone {
    int  tz_minuteswest; // minutes W of Greenwich
    int  tz_dsttime;     // type of dst correction
};

static inline int gettimeofday(struct timeval *tv, struct timezone *tz) {
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv) {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        //converting file time to unix epoch
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tmpres /= 10;  //convert into microseconds
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz) {
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

static inline unsigned long long nanoClock (void) { // actually returns microseconds
    struct timeval t;
    gettimeofday(&t, (struct timezone*)0);
    return t.tv_usec + t.tv_sec * 1000000ULL;
}

namespace RDKit {
    namespace FMCS {

#ifdef VERBOSE_STATISTICS_ON

// compute statistics of really very very fast calls.
// It a bit decrease overal performance, but might be interested for investigation purpose (only)
//#define VERBOSE_STATISTICS_FASTCALLS_ON

        struct ExecStatistics {
            unsigned TotalSteps, MCSFoundStep;
            unsigned long long   MCSFoundTime;
            unsigned InitialSeed, MismatchedInitialSeed;
            unsigned Seed, RemainingSizeRejected;
            unsigned SeedCheck, SingleBondExcluded;
            unsigned MatchCall, MatchCallTrue;
            unsigned FastMatchCall, FastMatchCallTrue, SlowMatchCallTrue;
            unsigned ExactMatchCall, ExactMatchCallTrue;   // hash cache
            unsigned FindHashInCache, HashKeyFoundInCache;
            unsigned AtomCompareCalls, BondCompareCalls;    // long long
            unsigned AtomFunctorCalls, BondFunctorCalls;    // long long
            unsigned WrongCompositionRejected, WrongCompositionDetected;
            unsigned DupCacheFound, DupCacheFoundMatch;

            ExecStatistics() : TotalSteps(0), MCSFoundStep(0)
                , MCSFoundTime(nanoClock())
                , InitialSeed(0), MismatchedInitialSeed(0)
                , Seed(0), RemainingSizeRejected(0)
                , SeedCheck(0), SingleBondExcluded(0), MatchCall(0), MatchCallTrue(0)
                , FastMatchCall(0), FastMatchCallTrue(0), SlowMatchCallTrue(0)
                , ExactMatchCall(0), ExactMatchCallTrue(0)
                , FindHashInCache(0), HashKeyFoundInCache(0)
                , AtomCompareCalls(0), BondCompareCalls(0)
                , AtomFunctorCalls(0), BondFunctorCalls(0)
                , WrongCompositionRejected(0), WrongCompositionDetected(0)
                , DupCacheFound(0), DupCacheFoundMatch(0)
            {}
        };
#endif

    }
}
