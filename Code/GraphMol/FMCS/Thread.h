// $Id$
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
#ifndef __THREAD_H__INCLUDED
#define __THREAD_H__INCLUDED
/**
@file   Thread.h
@brief  Cross-platform (POSIX/Windows) thread.
*/
#ifdef WIN32
    #include <windows.h>
    #include <process.h>
#else
    #include <unistd.h>
    #include <pthread.h>
    #include <sys/time.h>
#endif
#include <time.h>
#include <errno.h>
#include <string>
#include <cstring>
#include <stdio.h>

namespace RDKit
{
 namespace FMCS
 {
#ifndef WIN32 // POSIX.
    static inline void Sleep(int ms)
    {
        struct timespec t;
        t.tv_sec = ms/1000;
        t.tv_nsec=(ms%1000)*1000000L;
        while(-1==nanosleep(&t,&t) && errno==EINTR)
            sleep(0);
        sleep(0);
    }

    static inline void Sleep0()
    {
        usleep(0);
    }
#else
    static inline void Sleep0()
    {
        Sleep(0);
    }
#endif


    static inline long getThreadId()    // move it to class Thread
    {
        #if defined(WIN32)
            return (long)GetCurrentThreadId();
        #elif defined(_POSIX_THREADS)
            return (long)pthread_self();
        #else
            #error "You must define thread operations appropriate for your platform"
        #endif
    }


class SyncEvent;

class Mutex
{
private:
    #ifdef WIN32
	    CRITICAL_SECTION    mutex;
    #else // POSIX compatible OS
    	pthread_mutex_t     mutex;
        friend class        SyncEvent;  // pthread_cond_timedwait() need access to mutex handle
    #endif

public:
    inline Mutex()
    {
    #ifdef WIN32
    	InitializeCriticalSection(&mutex);
    #else // POSIX compatible OS

        // MAC/SUN OS - default type is PTHREAD_MUTEX_NORMAL:
        pthread_mutexattr_t attr;
        pthread_mutexattr_init   (&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);

        pthread_mutex_init(&mutex,&attr);

        pthread_mutexattr_destroy(&attr);
        //throw exception if error occured during initializing
    #endif
    }

    inline ~Mutex()
    {
    #ifdef WIN32
		DeleteCriticalSection(&mutex);
    #else // POSIX compatible OS
        pthread_mutex_destroy(&mutex);
    #endif
    }

	inline bool lock()
    {
    #ifdef WIN32
        EnterCriticalSection(&mutex);
        return true;
    #else // POSIX compatible OS
        return 0==pthread_mutex_lock(&mutex);
    #endif
    }

	inline bool trylock()
    {
    #ifdef WIN32
        return TRUE==::TryEnterCriticalSection(&mutex); // NT only
    #else // POSIX compatible OS
        return 0==pthread_mutex_trylock(&mutex);
    #endif
    }

    inline void unlock()
    {
    #ifdef WIN32
        LeaveCriticalSection(&mutex);
    #else // POSIX compatible OS
        pthread_mutex_unlock(&mutex);
    #endif
    }
};
//---------------------------------------------------------------------

    template<class SynchT=Mutex>
    class Lock
    {
        SynchT&  Synch;
    public:
        inline explicit Lock (SynchT& mutex) : Synch(mutex)
        {
            Synch.lock(); 
        }

        inline ~Lock () 
        {
            Synch.unlock();
        }
    };

    typedef Lock<Mutex> Guard;
//---------------------------------------------------------------------

    class Thread
    {
        #ifdef WIN32
            HANDLE      handle;
        #else // POSIX compatible OS
            pthread_t   handle;
        #endif
    public:
        inline Thread() : handle(0)
        {
        }

        virtual ~Thread()
        {
            join();
            handle = 0;
        }

        inline bool startThread(unsigned stackSize=0)
        {
        #ifdef WIN32
            #ifdef _WIN32_WCE
                DWORD id;
	            handle = CreateThread(NULL, stackSize, &Thread::_run, this, 0, &id);
                if(NULL==handle)
                    return false;
            #else
            {
                unsigned threadID;
                handle =  (HANDLE) _beginthreadex(NULL, stackSize, (unsigned int (__stdcall *)(void *)) &Thread::_run,this, 0, &threadID);
            }
            #endif
        #else // POSIX compatible OS
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            if(0!=stackSize)
                pthread_attr_setstacksize(&attr, stackSize);
            int r = pthread_create(&handle, &attr, _run, this);
            pthread_attr_destroy(&attr);
            if(0 != r)
                return false;
        #endif
            return true;
        }

        enum Priority
        {
            IDLE,
            LOW,
            NORMAL,
            HIGH,
            REALTIME
        };

        inline void setThreadPriority(Priority p)
        {
        #ifdef WIN32
            int val[5] = 
            {
                THREAD_PRIORITY_IDLE,
                THREAD_PRIORITY_LOWEST,
                THREAD_PRIORITY_NORMAL,
                THREAD_PRIORITY_HIGHEST,
                THREAD_PRIORITY_TIME_CRITICAL
            };
	            SetThreadPriority(handle, val[p]);
        #else // POSIX compatible OS
            int policy, min_priority, max_priority;
            struct sched_param param;
            pthread_getschedparam(pthread_self(), &policy, &param);
            min_priority = sched_get_priority_min(policy);
            max_priority = sched_get_priority_max(policy);
            switch(p)
            {
                default:
                case IDLE:
                    param.sched_priority = min_priority; 
                    break;
                case LOW:
                    param.sched_priority = min_priority + (max_priority - min_priority) * 3/10; 
                    break;
                case NORMAL:
                    param.sched_priority = min_priority + (max_priority - min_priority) / 2; 
                    break;
                case HIGH:
                    param.sched_priority = min_priority + (max_priority - min_priority) * 9/10; 
                    break;
                case REALTIME:
                    param.sched_priority = max_priority; 
                    break;
            };
            pthread_setschedparam(pthread_self(), policy, &param);
        #endif
        }

    protected:
        virtual void run(void)=0;   // must be overriden
    private:
        static inline 
        #ifdef WIN32
            #ifdef _WIN32_WCE
                DWORD WINAPI _run(void *param)
            #else
                void _run(void *param)
            #endif
        #else // POSIX compatible OS
                void* _run(void *param)
        #endif
        {
            ((Thread*)param)->run();
            #if !defined(WIN32) || defined(_WIN32_WCE)
                return 0;
            #endif
        }
    public:
        inline void join()
        {
        #ifdef WIN32
            if(-1!=(int)handle && NULL!=handle && (handle))
            {
                WaitForSingleObject(handle,10*1000/*INFINITE*/);    //!=WAIT_FAILED);
                CloseHandle(handle);
                handle = 0;
            }
        #else // POSIX compatible OS
            void* exitcode;
            if(0!=handle)
                pthread_join(handle, &exitcode);
        #endif
        }
    };
//---------------------------------------------------------------------
    class SyncEvent
    {
    #ifdef WIN32
        HANDLE handle;
    #else // POSIX compatible OS
        Mutex           mutex;
        pthread_cond_t  condition;
        bool            signalled;
    #endif
    public:
        inline SyncEvent()
        {
        #ifdef WIN32
            handle = CreateEvent(0, FALSE, FALSE, 0);
        #else // POSIX compatible OS
            signalled = false;
            pthread_cond_init(&condition, NULL);
        #endif
        }

        inline ~SyncEvent() 
        {
        #ifdef WIN32
            CloseHandle(handle);
        #else // POSIX compatible OS
            pthread_cond_destroy(&condition);
        #endif
        }

        inline bool wait(unsigned microseconds) // microseconds
        {
        #ifdef WIN32
                return WAIT_OBJECT_0 == WaitForSingleObject(handle, microseconds/1000+(microseconds<1000 ? 1:0));
        #else // POSIX compatible OS
            struct timeval tv;
            gettimeofday(&tv, (struct timezone*)0);
            unsigned long long t = tv.tv_sec * 1000000ULL + tv.tv_usec + microseconds;
            timespec timeout;
            timeout.tv_sec  = time_t(t/1000);
            timeout.tv_nsec = long  (t%1000 * 1000000UL);
            mutex.lock();
            while (!signalled)
                if(ETIMEDOUT == pthread_cond_timedwait(&condition, &mutex.mutex, &timeout))
                {
                    mutex.unlock();
                    return false;
                }
            signalled = false;
            mutex.unlock();
            return true;
        #endif
        }

        inline void setSignal () 
        {
        #ifdef WIN32
            SetEvent(handle);
        #else // POSIX compatible OS
            mutex.lock();
            signalled = true;
            mutex.unlock();
            pthread_cond_signal(&condition);
        #endif
        }
    };
//---------------------------------------------------------------------

}}
#endif//__THREAD_H__INCLUDED
