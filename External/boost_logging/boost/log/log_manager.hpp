// log_manager.hpp

// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 John Torjo (john@torjo.com)
// Copyright (C) 2004-2005 Darryl Green (darryl.green@aanet.com.au)
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


#ifndef JT_BOOST_LOG_LOG_MANAGER_HPP
#define JT_BOOST_LOG_LOG_MANAGER_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#ifndef JT_BOOST_LOG_LOG_EXT_HPP
#error Donot include this directly. Include <boost/log/log_impl.hpp> instead.
#endif


namespace boost { namespace logging {

struct logger_impl;

// may be used from functions:
// - boost::logging::add_appender
// - boost::logging::add_modifier
// - boost::logging::enable_logs

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Required functions

BOOST_LOG_DECL void add_appender(default_log_manager & manager, const logging_types::log_name_string_type & str, appender_func f, const logging_types::log_name_string_type & name, int idx );
// adds a function that modifies messages (like, adds a prefix and/or suffix)
BOOST_LOG_DECL void add_modifier(default_log_manager & manager, const logging_types::log_name_string_type & str, modifier_func f, const logging_types::log_name_string_type & name , int idx );

BOOST_LOG_DECL void del_modifier(default_log_manager & manager, const logging_types::string & spec, const logging_types::log_name_string_type & name) ;
BOOST_LOG_DECL void del_appender(default_log_manager & manager, const logging_types::string & spec, const logging_types::log_name_string_type & name) ;


BOOST_LOG_DECL void enable_logs(default_log_manager & manager, const logging_types::log_name_string_type & logs_spec, level_type lvl);

BOOST_LOG_DECL void write_msg(default_log_manager & manager, logger & l, const logging_types::string & msg, level_type lvl);

BOOST_LOG_DECL boost::shared_ptr<logger_impl> find_log_by_name( default_log_manager & manager, const logging_types::log_name_string_type & log_name);

typedef boost::function2<void, const logging_types::string &, level_type> on_exception_writer_func;
BOOST_LOG_DECL void on_exception_writer(default_log_manager & , on_exception_writer_func func) ;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Default-manager functions

const int DEFAULT_CACHE_LIMIT = 256;

BOOST_LOG_DECL void set_log_caching(bool do_cache, int cache_limit = DEFAULT_CACHE_LIMIT );
BOOST_LOG_DECL void flush_log_cache();


}}


#endif
