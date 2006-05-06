// log_impl.hpp

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

#ifndef JT_BOOST_LOG_LOG_EXT_HPP
#define JT_BOOST_LOG_LOG_EXT_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif


#include <boost/log/detail/defs.hpp>
#include <boost/log/log.hpp>
#include <sstream>
#include <vector>

/*
    Note: you might wonder what's the purpose of this header file:
    why not include all this in log.hpp?

    The reason is that internally I'm using Boost.Function - for Appenders/Modifier Functions.
    #including <boost/function.hpp> is time-costly. 
    
    So better include this header only when you're manipulating the logs.
    When you're writing to the logs, you'll include only <boost/log/log.hpp>, which compiles VERY FAST.
*/

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/log/log.hpp>

namespace boost { namespace logging {

// logs a message
typedef boost::function<void (const logging_types::string&,const logging_types::string&)> appender_func;
// modifies a message
typedef boost::function<void (const logging_types::string&,logging_types::string&)> modifier_func;

}}

#include <boost/log/log_manager.hpp>

namespace boost { namespace logging {

struct BOOST_LOG_DECL logger_impl {
    friend struct logger;

    logger_impl(const logging_types::log_name_string_type & name);
    ~logger_impl();

    void add_appender(appender_func a, const logging_types::log_name_string_type & name, int idx);
    void add_modifier(modifier_func m, const logging_types::log_name_string_type & name, int idx);

    void del_appender(const logging_types::log_name_string_type & name);
    void del_modifier(const logging_types::log_name_string_type & name);

    void enable(level_type lvl);
    bool is_enabled(level_type lvl) const;

    void write_msg(const logging_types::string & msg, level_type lvl);

    void set_is_compile_time(bool val);

private:
    struct appender_info {
        appender_info( const logging_types::log_name_string_type & name, int idx) : name(name), idx(idx) {}
        appender_func appender;
        logging_types::log_name_string_type name;
        int idx;
    };
public:
    typedef std::vector<appender_info> appender_array;
private:
    appender_array m_appenders;

    struct modifier_info {
        modifier_info(const logging_types::log_name_string_type & name, int idx) : name(name), idx(idx) {}
        modifier_func modifier;
        logging_types::log_name_string_type name;
        int idx;
    };
public:
    typedef std::vector<modifier_info> modifier_array;
private:
    modifier_array m_modifiers;

public:
    const appender_array & appenders() const { return m_appenders; }
    const modifier_array & modifiers() const { return m_modifiers; }

    // the name of this log; example: "app.gui.charts"
    //
    // note: the name is ALWAYS a logging_types::log_name_string_type, no matter what strings you're writing!
    const logging_types::log_name_string_type name;

private:
    // if writing a message with this level or a higher level, it's logged
    // otherwise, not
    level_type m_level;

    // if true, its enabled state can only be set at compile-time.
    bool m_is_compile_time;

    logging_types::mutex cs;
};




const int DEFAULT_INDEX = 100;
// adds a function that logs messages
inline void add_appender(const logging_types::string & str, appender_func f,  const logging_types::log_name_string_type & name = logging_types::log_name_string_type(), int idx = DEFAULT_INDEX) {
    add_appender( logger::manager(), str, f, name, idx);
}

// adds a function that modifies messages (like, adds a prefix and/or suffix)
//
// name - the name of the modifier. So that it could eventually be later on removed
inline void add_modifier(const logging_types::string & spec, modifier_func f, const logging_types::log_name_string_type & name = logging_types::log_name_string_type(), int idx = DEFAULT_INDEX) {
    add_modifier( logger::manager(), spec, f, name, idx);
}

inline void del_modifier(const logging_types::string & spec, const logging_types::log_name_string_type & name) {
    del_modifier( logger::manager(), spec, name);
}

inline void del_appender(const logging_types::string & spec, const logging_types::log_name_string_type & name) {
    del_appender( logger::manager(), spec, name);
}


inline void enable_logs(const logging_types::string & logs_spec, level_type lvl = level::default_) {
    enable_logs( logger::manager(), logs_spec, lvl);
}

// helper
inline void disable_logs(const logging_types::string & logs_spec) {
    enable_logs( logger::manager(), logs_spec, level::disable_all);
}

inline void on_exception_writer(on_exception_writer_func func) {
    on_exception_writer( logger::manager(), func);
}


}}

#include <boost/log/manipulate_logs.hpp>


#endif
