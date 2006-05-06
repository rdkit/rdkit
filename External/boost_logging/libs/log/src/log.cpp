// log.cpp 

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
    Note: I supplied this as a static library, because otherwise 
    compilation times would be too big.
*/


#include <boost/log/log.hpp>
#include <boost/log/log_impl.hpp>
#include <sstream>
#include <algorithm>

namespace boost { namespace logging {

BOOST_LOG_DECL logger_impl & logger_to_logger_impl(logger& l) {
    // FIXME - in case the logger is destroyed, make sure we return something else!
    return *l.m_impl.get();
}


logging_types::string logger::name() const {
    return m_impl->name;
}

bool logger::is_enabled(level_type lvl) const {
    return !m_destroyed && m_impl->is_enabled(lvl);
}

void logger::enable( level_type lvl) {
    if ( !m_destroyed)
        if ( m_impl)
            m_impl->enable(lvl);
}

bool logger::still_exists() const {
    return !m_destroyed;
}


logger::~logger() {
    m_destroyed = true;
    // out-of-line, so that the code for m_impl's destruction goes here 
    // had this been inline, the logger would need logger_impl's definition...
}

namespace {
    // sort by index
    struct by_idx {
        template<class type> bool operator()(const type & first, const type & second) {
            return first.idx < second.idx;
        }
    };
}

logger_impl::logger_impl(const logging_types::log_name_string_type & name) : name(name), m_is_compile_time(false), m_level(level::default_) {
}
logger_impl::~logger_impl() {
    m_level = level::disable_all;
}

bool logger_impl::is_enabled(level_type lvl) const {
    // note: intentionally not locked (for speed)
    return lvl >= m_level;
}


void logger_impl::add_appender(appender_func a, const logging_types::log_name_string_type & name, int idx) {
    logging_types::lock lk(cs);

    appender_info info(name, idx);
    appender_array::iterator at = std::lower_bound( m_appenders.begin(), m_appenders.end(), info, by_idx() );
    info.appender = a;
    m_appenders.insert(at, info);
}

void logger_impl::add_modifier(modifier_func m, const logging_types::log_name_string_type & name, int idx) {
    logging_types::lock lk(cs);

    modifier_info info(name, idx);
    modifier_array::iterator at = std::lower_bound( m_modifiers.begin(), m_modifiers.end(), info, by_idx() );
    info.modifier = m;
    m_modifiers.insert(at, info);
}



namespace {
    struct has_name {
        has_name(const logging_types::log_name_string_type & name) : name(name) {}
        template<class modifier_or_appender> bool operator()(const modifier_or_appender & val) {
            return val.name == name;
        }
        logging_types::log_name_string_type name;
    };
}

// deletes the appender(s) that have this name
void logger_impl::del_appender(const logging_types::log_name_string_type & name) {
    logging_types::lock lk(cs);
    m_appenders.erase( 
        std::remove_if( m_appenders.begin(), m_appenders.end(), has_name(name)), 
        m_appenders.end() );
}

void logger_impl::del_modifier(const logging_types::log_name_string_type & name) {
    logging_types::lock lk(cs);
    m_modifiers.erase( 
        std::remove_if( m_modifiers.begin(), m_modifiers.end(), has_name(name)), 
        m_modifiers.end() );
}


void logger_impl::enable(level_type lvl) {
    logging_types::lock lk(cs);
    if ( m_is_compile_time)
        // compile-time logs: we always consider them enabled;
        // if a compile-time log is disabled, the whole "BOOST_LOG(..) << ..." line 
        // is optimized away while pre-processing
        m_level = level::enable_all;
    else
        m_level = lvl;
}

void logger_impl::set_is_compile_time(bool val) {
    logging_types::lock lk(cs);
    m_is_compile_time = val;
    // compile-time logs: we always consider them enabled
    // if a compile-time log is disabled, the whole "BOOST_LOG(..) << ..." line 
    // is optimized away while pre-processing
    if ( m_is_compile_time)
        m_level = level::enable_all;
}


void logger_impl::write_msg(const logging_types::string & msg, level_type lvl) {
    logging_types::string copy = msg; // it could be modified...

    logging_types::lock lk(cs);
    if ( !is_enabled(lvl) )
        if ( !m_is_compile_time)
            // now, we're locked - thus, really make sure we're enabled... (in other words - the is_enabled() test
            // is not locked - for speed, but this one is)
            return; 
        else
            ; // this is a compile-time enabled log
              // if we ended up here, it means we're enabled at compiled time

    // run modifiers
    { typedef modifier_array array;
      for (array::const_iterator b = modifiers().begin(), e = modifiers().end(); b != e; ++b)
        (b->modifier)( name, copy); }

    // run appenders
    { typedef appender_array array;
      for (array::const_iterator b = appenders().begin(), e = appenders().end(); b != e; ++b)
        (b->appender)( name, copy); }

}


}}

