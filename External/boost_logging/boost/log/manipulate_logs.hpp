// manipulate_logs.hpp

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

/*
    Note: the library is thread-safe.
*/

#ifndef JT_BOOST_LOG_manipulate_logs_hpp
#define JT_BOOST_LOG_manipulate_logs_hpp

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/log/detail/defs.hpp>
#include <boost/config.hpp>
#include <boost/log/detail/ts.hpp>



namespace boost { namespace logging {


/** 
    Helper class, allow easier manipulation of logs:
    - enable/disabling logs
    - add/del appenders
    - add/del modifiers

    You can chain multiple manipulation commands, like this:

    @code
    manipulate_logs("*")
        .add_modifier(prepend_time("$hh:$mm:$ss "), DEFAULT_INDEX + 1 )
        .add_modifier(&prepend_prefix)
        .add_appender(&write_to_cout, "write_to_cout");
    @endcode

*/
struct manipulate_logs {
    typedef manipulate_logs & self;
    manipulate_logs(const logging_types::log_name_string_type & spec) : spec(spec) {}

    self enable(level_type lvl = level::default_) {
        enable_logs(spec, lvl);
        return *this;
    }
    self disable() {
        disable_logs(spec);
        return *this;
    }

    self add_appender(appender_func f,  const logging_types::log_name_string_type & name = logging_types::log_name_string_type(), int idx = DEFAULT_INDEX) {
        ::boost::logging::add_appender(spec, f, name, idx);
        return *this;
    }
    self add_appender(appender_func f, int idx ) {
        ::boost::logging::add_appender(spec, f, logging_types::log_name_string_type(), idx);
        return *this;
    }

    self add_modifier(modifier_func f,  const logging_types::log_name_string_type & name = logging_types::log_name_string_type(), int idx = DEFAULT_INDEX) {
        ::boost::logging::add_modifier(spec, f, name, idx);
        return *this;
    }
    self add_modifier(modifier_func f, int idx ) {
        ::boost::logging::add_modifier(spec, f, logging_types::log_name_string_type(), idx);
        return *this;
    }

    self del_appender(const logging_types::log_name_string_type & name) {
        ::boost::logging::del_appender(spec, name);
        return *this;
    }
    self del_modifier(const logging_types::log_name_string_type & name) {
        ::boost::logging::del_modifier(spec, name);
        return *this;
    }

private:
    logging_types::log_name_string_type spec;
};

// in case you prefer to type less...
typedef manipulate_logs log_sink;
}}


#endif
