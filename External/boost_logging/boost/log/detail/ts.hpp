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

#ifndef JT_BOOST_LOG_TS_HPP
#define JT_BOOST_LOG_TS_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif


#ifdef BOOST_HAS_THREADS
    
    #ifdef BOOST_LOG_USE_BOOST_THREADS
        #include <boost/log/detail/ts_boost.hpp>
    #else
        #ifdef BOOST_LOG_WIN32
        #include <boost/log/detail/ts_win32.hpp>
        #else
        #include <boost/log/detail/ts_posix.hpp>
        #endif
    #endif

#else
    // no threads
    #include <boost/log/detail/ts_none.hpp>
#endif



#endif

