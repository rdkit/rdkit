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

#ifndef JT_BOOST_LOG_TS_boost_HPP
#define JT_BOOST_LOG_TS_boost_HPP


#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/thread/mutex.hpp>


namespace boost { namespace logging {

namespace threading {

    typedef boost::mutex mutex;
    typedef mutex::scoped_lock scoped_lock;
}

}}


#endif

