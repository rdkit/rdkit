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

#ifndef JT_BOOST_LOG_TS_WIN32_HPP
#define JT_BOOST_LOG_TS_WIN32_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace boost { namespace logging { namespace threading {

class mutex {
    mutex & operator = ( const mutex & Not_Implemented);
    mutex( const mutex & From);
public:
    mutex() {
        InitializeCriticalSection( GetCriticalSectionPtr() );
    }
    void Lock() {
        EnterCriticalSection( GetCriticalSectionPtr());
    }
    void Unlock() {
        LeaveCriticalSection( GetCriticalSectionPtr());
    }
private:
    LPCRITICAL_SECTION GetCriticalSectionPtr() const { return &m_cs; }
    mutable CRITICAL_SECTION m_cs;
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

}}}

#endif

