#pragma once
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <iostream>
#ifdef WIN32
    #include <Windows.h> // for Winmm.lib timeGetTime()
#else
    #include <unistd.h>
    #include <fcntl.h>
#endif

// SELECT ALGORITHM OPTIONS by comment some lines to exclude additional or experimental optimisations:

#define SEED_GROW_DEEP              // fast and works much times faster (but it can depend on molecules)
//#define EXCLUDE_WRONG_COMPOSITION   // fast but with a little effect, because amount of external bonds usually is very small.
                                    // Exclude mismatched bonds combinations during seed growing (2^N stage)

#define FAST_SUBSTRUCT_CACHE        // based on Morgan code hash
#define DUP_SUBSTRUCT_CACHE         // based on list of query atoms and bonds. For rings where seeds growing in both directions throw the same ring.

#define PRECOMPUTED_TABLES_MATCH    // Improves overal performance about 20%, especially in hard cases.
                                    // Takes some extra memory (Vt*vq+Et*eq)/8 bytes for each target G(Vt,Et) matched with query G(vq,eq).

//#define FAST_INCREMENTAL_MATCH      // NOT IMPLEMENTED YET fast and should be usefull
                                    // history based check match without finding new matched substructure location in the target


#define VERBOSE_STATISTICS_ON

// Enable / Disable DEBUG TRACE output
#ifdef WIN32__xx__TRACE_ON
    #define TRACE_ON
#else
#endif

#ifdef VERBOSE_STATISTICS_ON
// compute statistics of really very very fast calls. 
// It a bit decrease overal performance, but might be interested for investigation purpose (only)
//#define VERBOSE_STATISTICS_FASTCALLS_ON

    struct ExecStatistics
    {
        unsigned TotalSteps, MCSFoundStep;
        unsigned InitialSeed, MismatchedInitialSeed;
        unsigned Seed, RemainingSizeRejected;
        unsigned SeedCheck;
        unsigned MatchCall, MatchCallTrue;
        unsigned FastMatchCall, FastMatchCallTrue;
        unsigned ExactMatchCall, ExactMatchCallTrue;   // hash cache
        unsigned FindHashInCache, HashKeyFoundInCache;
#ifdef SMILES_CACHE
        unsigned FindInMatchedCache, FoundInMatchedCache;
        unsigned FindInDoNotMatchedCache, FoundInDoNotMatchedCache;
#endif
        unsigned AtomCompareCalls, BondCompareCalls;    // long long
        unsigned AtomFunctorCalls, BondFunctorCalls;    // long long
        unsigned WrongCompositionRejected, WrongCompositionDetected;
        unsigned DupCacheFound, DupCacheFoundMatch;

        ExecStatistics() : TotalSteps(0), MCSFoundStep(0)
                         , InitialSeed(0), MismatchedInitialSeed(0)
                         , Seed(0), RemainingSizeRejected(0)
                         , SeedCheck(0), MatchCall(0), MatchCallTrue(0)
                         , FastMatchCall(0), FastMatchCallTrue(0)
                         , ExactMatchCall(0), ExactMatchCallTrue(0)
                         , FindHashInCache(0), HashKeyFoundInCache(0)
#ifdef SMILES_CACHE
                         , FindInMatchedCache(0), FoundInMatchedCache(0)
                         , FindInDoNotMatchedCache(0), FoundInDoNotMatchedCache(0)
#endif
                         , AtomCompareCalls(0), BondCompareCalls(0)
                         , AtomFunctorCalls(0), BondFunctorCalls(0)
                         , WrongCompositionRejected(0), WrongCompositionDetected(0)
                         , DupCacheFound(0), DupCacheFoundMatch(0)
        {}
    };
    extern ExecStatistics extstat;
#endif

namespace RDKit
{
 namespace FMCS
 {
     extern bool ConsoleOutputEnabled;
#ifdef TRACE_ON
     std::ostream& TRACE(unsigned info=1)
     {
         if(info > 0)   // print out time
         {
            struct tm *now;
            char strt[32], strts[12];
            #ifdef WIN32
                static unsigned long t0 = timeGetTime();                        see gettimeofday for Win in test.cpp
                typedef struct timespec
                {
	                time_t tv_sec;  // Seconds since 00:00:00 GMT = 1 January 1970
	                long   tv_nsec; // Additional nanoseconds since Windows started - assume since tv_sec
                } timespec;

                struct timespec ts;
	            //time( &(ts.tv_sec) );
                unsigned long t = timeGetTime() - t0; // t from app started
                ts.tv_sec  =  t / 1000LU;
	            ts.tv_nsec = (t % 1000LU)*1000000LU;     // nanoseconds.
            #else
                struct timespec ts;
                //clock_gettime(CLOCK_REALTIME, &ts);
                clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
            #endif

            now = gmtime(&ts.tv_sec);
            if(now != NULL)
                strftime(strt, sizeof(strt), "%H:%M:%S", now); //strftime(strt, sizeof(strt), "%Y-%m-%d %H:%M:%S", now);

            #ifdef WIN32
                sprintf(strts, ".%03lu ", ts.tv_nsec / 1000000LU);
            #else
                sprintf(strts, ".%05lu ", ts.tv_nsec / 10000LU);
            #endif

            std::cout << strt << strts;
        }
        return std::cout;
     }
#else
     class NullStream
     {
     public:
        template<class T>
        NullStream& operator << (const T&) {return *this;}  // DO NOTHING
     };
     static inline NullStream& TRACE(unsigned info=1) { static NullStream z; return z; }
#endif
}}
