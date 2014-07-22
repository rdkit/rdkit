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
#include <vector>
#include <stdexcept>
#include "FMCS.h"
#include "Seed.h"
#include "Thread.h"


namespace RDKit {
    namespace FMCS {
        class MaximumCommonSubgraph;

        class SeedGrowThread : public Thread {
            bool            StopSignal;
            SyncEvent*      NewJobEvent;
            SyncEvent*      JobFinishedPoolEvent;
            const Seed*     SeedToProceed;
            MaximumCommonSubgraph* Mcs;
            virtual void run(void);

        public:
            SeedGrowThread () : StopSignal(false), NewJobEvent(0), JobFinishedPoolEvent(0), SeedToProceed(0), Mcs(0) {}
            ~SeedGrowThread() {
                stop();
                join();
                delete NewJobEvent;
            }
            void setParameters(MaximumCommonSubgraph* mcs, SyncEvent* jobFinishedPoolEvent) {
                Mcs = mcs;
                NewJobEvent = new SyncEvent(); // avoid making a copy of pthread condition object
                JobFinishedPoolEvent = jobFinishedPoolEvent;
            }
            bool isBusy()const              {
                return SeedToProceed != 0;
            }
            void setJob(const Seed* seed)   {
                SeedToProceed = seed;
                NewJobEvent->setSignal();
            }
            void stop() {
                StopSignal = true;
                if(NewJobEvent) NewJobEvent->setSignal();
            }
        };

        class ThreadPool {
            bool                        StopSignal;
            std::vector<SeedGrowThread> ThreadSet;
            SyncEvent                   JobFinishedEvent;
        public:
            ThreadPool(size_t size, MaximumCommonSubgraph* mcs) : StopSignal(false), ThreadSet(size) {
                for(std::vector<SeedGrowThread>::iterator s = ThreadSet.begin(); s != ThreadSet.end(); s++) {
                    s->setParameters(mcs, &JobFinishedEvent);
                    s->startThread();
                    s->setThreadPriority(Thread::HIGH);
                }
                JobFinishedEvent.setSignal();   // ready to work
            }
            ~ThreadPool() {
                StopSignal = true;
                for(std::vector<SeedGrowThread>::iterator s = ThreadSet.begin(); s != ThreadSet.end(); s++)
                    s->stop();
                ThreadSet.clear();  // join() before JobFinishedEvent will be destroyed
            }
            bool addJob(const Seed* seed) {
                while(!StopSignal) {
                    for(std::vector<SeedGrowThread>::iterator s = ThreadSet.begin(); s != ThreadSet.end(); s++)
                        if(!s->isBusy()) {
                            s->setJob(seed);
                            return true;
                        }
                    JobFinishedEvent.wait(9000);//Sleep0(); // actually never
                }
                return false; // never
            }

            bool waitReady(unsigned milliseconds) {
                return JobFinishedEvent.wait(milliseconds);
            }

            bool isBusy()const {
                for(std::vector<SeedGrowThread>::const_iterator s = ThreadSet.begin(); s != ThreadSet.end(); s++)
                    if(!s->isBusy())
                        return false;
                return true;
            }

        };
    }
} // namespace RDKit
