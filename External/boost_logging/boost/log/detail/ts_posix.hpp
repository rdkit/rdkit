// defs.hpp

// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 John Torjo (john@torjo.com)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.


// Copyright (C) 2001-2003
// William E. Kempf
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  William E. Kempf makes no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.


#ifndef JT_BOOST_LOG_TS_posix_HPP
#define JT_BOOST_LOG_TS_posix_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <errno.h>
#include <pthread.h>
#include <stdexcept>
#include <cassert>

namespace boost { namespace logging {

namespace threading {

class mutex {
    mutex & operator = ( const mutex & Not_Implemented);
    mutex( const mutex & From);
public:
    mutex() : m_mutex() {
        pthread_mutexattr_t attr;
        int res = pthread_mutexattr_init(&attr);
        assert(res == 0);

        res = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        assert(res == 0);

        res = pthread_mutex_init(&m_mutex, &attr);
        {
            int res = pthread_mutexattr_destroy(&attr);
            assert(res == 0);
        }
        if (res != 0)
            throw std::runtime_error("could not create mutex");
    }
    ~mutex() {
        int res = 0;
        res = pthread_mutex_destroy(&m_mutex);
        assert(res == 0);
    }

    void Lock() {
        int res = 0;
        res = pthread_mutex_lock(&m_mutex);
        assert(res == 0);
        if (++m_count > 1)
        {
            res = pthread_mutex_unlock(&m_mutex);
            assert(res == 0);
        }
    }
    void Unlock() {
        if (--m_count == 0)
        {
            int res = 0;
            res = pthread_mutex_unlock(&m_mutex);
            assert(res == 0);
        }
    }
private:
    pthread_mutex_t m_mutex;
    unsigned m_count;
};

class scoped_lock {
    scoped_lock operator=( scoped_lock & Not_Implemented);
    scoped_lock( const scoped_lock & Not_Implemented);
public:
    scoped_lock( mutex & cs) : m_cs( cs)                { m_cs.Lock(); }
    ~scoped_lock()                                      { m_cs.Unlock(); }
private:
    mutex & m_cs;
};


}

}}


#endif

