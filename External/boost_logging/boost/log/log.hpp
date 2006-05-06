// log.hpp

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

// FIXME make it possible to write to rotating logs (make logger)

#ifndef JT_BOOST_LOG_LOG_HPP
#define JT_BOOST_LOG_LOG_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <sstream>
#include <boost/log/detail/defs.hpp>
#include <boost/log/log_fwd.hpp>



///////////////////////////////////////////////////////////////////////////////////////////////////
// Levels

#define BOOST_LOG_DEFINE_LEVEL(lvl, value) namespace boost { namespace logging { namespace level { const ::boost::logging::level_type lvl = (::boost::logging::level_type)(value); }}}

namespace boost { namespace logging {

// what we use to hold a level
typedef unsigned int level_type;
}}

BOOST_LOG_DEFINE_LEVEL(disable_all, -1)
BOOST_LOG_DEFINE_LEVEL(default_, 1000)
BOOST_LOG_DEFINE_LEVEL(enable_all, 0)

BOOST_LOG_DEFINE_LEVEL(fatal, 2000)
BOOST_LOG_DEFINE_LEVEL(err, 1600)
BOOST_LOG_DEFINE_LEVEL(dbg, 1400)
BOOST_LOG_DEFINE_LEVEL(warn, 1200)
BOOST_LOG_DEFINE_LEVEL(info, 1000)



namespace boost { namespace logging {

// forward declare
struct logger;
struct logger_impl;
struct enabled_logger;
struct logger_stream;

// types that are used for logging. You can override them.
namespace logging_types {

    typedef log_manager<0> log_manager_class;
    typedef log_manager_class::type log_manager_type;
    typedef log_manager_type::char_t char_t;
    typedef log_manager_type::string string;
    typedef log_manager_type::stream stream;

    // used to hold the name of a log. Note that you can choose a better 'string' type
    // as your string; however, log names are always held like this.
    //
    // having this typedef allows to move the implementation of logger_impl in a source file
    // (otherwise, it wouldn't be possible)
    typedef std::basic_string<char_t> log_name_string_type;

    typedef boost::logging::threading::mutex mutex;
    typedef boost::logging::threading::scoped_lock lock;

};


// forward declare - in case you provide your own log manager, you have to define it yourself...
BOOST_LOG_DECL void write_msg(logging_types::log_manager_type & manager, logger & l, const logging_types::string & msg, level_type lvl);
BOOST_LOG_DECL boost::shared_ptr<logger_impl> find_log_by_name( logging_types::log_manager_type & manager, const logging_types::log_name_string_type & log_name);

BOOST_LOG_DECL void mark_log_as_compile_time(const logging_types::log_name_string_type & log_name);


BOOST_LOG_DECL logger_impl & logger_to_logger_impl(logger&);

/* 
    represents a logger. It has a name, it can be enabled/disabled, 
    it can have appenders and/or modifiers

    multiple logs with the same name share the same impl.
*/
struct BOOST_LOG_DECL logger {
    enum { is_compile_time = 0 };

    friend struct enabled_logger;
    friend BOOST_LOG_DECL logger_impl & logger_to_logger_impl(logger&);

    logger(const logging_types::log_name_string_type & log_name) : m_destroyed(false) {
        m_impl = find_log_by_name( manager(), log_name);
    }

    logging_types::log_name_string_type name() const;
    ~logger();

    bool is_enabled(level_type lvl) const;
    void enable( level_type lvl);

    bool still_exists() const;

    bool operator()(level_type lvl) const {
        return is_enabled(lvl);
    }

    operator bool() const {
        return is_enabled(level::default_);
    }

    inline logger_stream stream(level_type lvl = level::default_);
private:
    // this creates a manager, shared by all log instances
    typedef logging_types::log_manager_type manager_type;
public:
    // FIXME - in the future, make private
    //
    // note: - this is just for overload resolution - in case you have your own manager (log_manager.hpp)
    //       - don't keep any persistent data in the manager class.
    static manager_type & manager() {
        static manager_type m;
        return m;
    }

private:
    boost::shared_ptr<logger_impl> m_impl;
    // simple way to detect if we've been destroyed...
    // (it's possible that while we're destroyed - other threads still try to write...)
protected:
    bool m_destroyed;
};


/*
    Helper - for keeping both a logger and the level of your message.
    Useful when defining the macros
*/
struct logger_and_level {
    logger_and_level(logger & l, level_type lvl) : l(l), lvl(lvl) {}
    mutable logger &l;
    level_type lvl;
};


/*
    allows writing to a log, only if it's enabled

    You should prefer using BOOST_LOG, but you can use the following as well:
    if ( enabled_logger l = some_log)
        l.stream() << whatever << you << wish << to << write;

    You should note that  '<< whatever << you << wish << to << write' will be executed
    *only if the logger is enabled!*

*/
struct enabled_logger {
    // when the logger is enabled at compile-time 
    struct force_enabled {};

    enabled_logger(logger & l, const force_enabled &) : m_stream(0), m_logger(l), m_level(0) {
        // logger is enabled at compile-time
        m_stream = new logging_types::stream ;
    }

    enabled_logger(const logger_and_level & l) : m_stream(0), m_logger(l.l), m_level(l.lvl) {
        if (  m_logger(m_level) ) 
            m_stream = new logging_types::stream ;
        else ; // log is not enabled - ignore
    }
    enabled_logger(const enabled_logger & copy) : m_stream(copy.m_stream), m_logger(copy.m_logger), m_level(copy.m_level) {
        copy.m_stream = 0;
    }

    ~enabled_logger() {
        // the actual writing of the message
        if ( m_stream) {
            try {
                write_msg( logger::manager(), m_logger, m_stream->str(), m_level );
            } catch(...) {
                // write_msg is guaranteed not to throw
                // so the only thing that can throw is m_stream->str() (couldn't allocate memory?)
                // if so, all bets are off
            }
            delete m_stream;
        }
    }

    operator bool() const {
        return m_stream != 0;
    }

    logging_types::stream * stream() const { return m_stream; }

private:
    mutable logging_types::stream * m_stream;
    logger & m_logger;
    level_type m_level;

private:
    enabled_logger & operator=(const enabled_logger&) const;
};


/*
    allows writing to a log, only if it's enabled

    Used from BOOST_LOG
*/
struct simple_logger_keeper : enabled_logger {
    simple_logger_keeper(const logger_and_level & l) : enabled_logger(l) {}

    operator bool() const {
        return ! (enabled_logger::operator bool() );
    }
};


template<bool is_compile_time, bool is_enabled > struct logger_keeper : enabled_logger {
    logger_keeper(const logger_and_level & l) : enabled_logger(l) {}

    operator bool() const {
        return ! (enabled_logger::operator bool() );
    }
};


// compile-time logger, enabled
template<> struct logger_keeper<true,true> : enabled_logger {
    logger_keeper(const logger_and_level & l) : enabled_logger(l.l, force_enabled() ) {}

    operator bool() const {
        return false;
    }
};

// compile-time logger, disabled
template<> struct logger_keeper<true,false> {
    logger_keeper(const logger_and_level & ) {}

    operator bool() const {
        return true;
    }
    logging_types::stream * stream() const { return 0; }
};



/*
    allows writing to a log

    You should prefer using BOOST_LOG, or enabled_logger, but you can use the following as well:

    logger l = ...;
    l.stream() << whatever << you << wish << to << write;

    You should note that in this case, '<< whatever << you << wish << to << write' will be executed
    *even if the log is disabled!*

    You can test if the log is enabled first:
    logger l = ...;
    if (l(some_level)) l.stream(some_level) << whatever << you << wish << to << write;

*/
struct logger_stream {

    logger_stream(logger & l, level_type lvl) : m_stream(0), m_logger(l), m_level(lvl) {
        m_stream = new logging_types::stream ;
    }
    logger_stream(const logger_stream & copy) : m_stream(copy.m_stream), m_logger(copy.m_logger), m_level(copy.m_level) {
        copy.m_stream = 0;
    }

    ~logger_stream() {
        // the actual writing of the message
        if ( m_stream) {
            try {
                write_msg( logger::manager(), m_logger, m_stream->str(), m_level );
            } catch(...) {
                // write_msg is guaranteed not to throw
                // so the only thing that can throw is m_stream->str() (couldn't allocate memory?)
                // if so, all bets are off
            }
            delete m_stream;
        }
    }

    operator bool() const {
        return m_stream != 0;
    }

    logging_types::stream * stream() const { return m_stream; }

private:
    mutable logging_types::stream * m_stream;
    logger & m_logger;
    level_type m_level;

private:
    logger_stream & operator=(const logger_stream&) const;
};


inline logger_stream logger::stream(level_type lvl) {
    return logger_stream(*this, lvl);
}



// helper - make BOOST_IS_LOG_ENABLED efficient, for compile-time logs
template<bool is_compile_time, bool is_enabl> struct logger_keeper_is_enabled {
    bool is_enabled(logger & l, level_type lvl) const { return l.is_enabled(lvl); }
};

template<> struct logger_keeper_is_enabled<true,true> {
    bool is_enabled(logger & l, level_type ) const { return true; }
};
template<> struct logger_keeper_is_enabled<true,false> {
    bool is_enabled(logger & l, level_type ) const { return false; }
};

template<bool is_compile_time> struct mark_log_as_compile_time_t {
    mark_log_as_compile_time_t(const logging_types::log_name_string_type & ) {}
};
template<> struct mark_log_as_compile_time_t<true> {
    // this is a compile-time log
    mark_log_as_compile_time_t(const logging_types::log_name_string_type & log_name) {
        mark_log_as_compile_time(log_name);
    }
};


}}



// ... for use in BOOST_DECLARE_LOG_DEBUG/ BOOST_DEFINE_LOG_DEBUG
#ifndef NDEBUG
#define BOOST_LOG_IS_DEBUG_MODE 1
#else
#define BOOST_LOG_IS_DEBUG_MODE 0
#endif


#define BOOST_DECLARE_LOG(log_name) \
    struct log_name ## _log_class { enum { is_compile_time = 0, is_enabled = 1 }; }; \
    ::boost::logging::logger & log_name();

// note: you need to include <boost/log/log_impl.hpp> before using this
#define BOOST_DEFINE_LOG(log_name,log_str) \
    ::boost::logging::logger & log_name() { \
    static ::boost::logging::mark_log_as_compile_time_t< (log_name ## _log_class::is_compile_time != 0) > mark_logger(log_str); \
    static ::boost::logging::logger l(log_str); return l; \
    }


#define BOOST_LOGL(log_name,lvl) \
    if (::boost::logging::logger_keeper<(log_name ## _log_class::is_compile_time != 0),(log_name ## _log_class::is_enabled != 0)> keep = ::boost::logging::logger_and_level(log_name(),::boost::logging::level::lvl) ) ; else (*keep.stream())

#define BOOST_LOG(log_name) BOOST_LOGL(log_name,default_)



#define BOOST_IS_LOG_ENABLED(log_name,lvl) \
    ::boost::logging::logger_keeper_is_enabled<(log_name ## _log_class::is_compile_time != 0),(log_name ## _log_class::is_enabled != 0)>().is_enabled( log_name(), ::boost::logging::level::lvl )
//log_name##_log().is_enabled()

#define BOOST_SCOPEDLOGL(log_name,lvl) if (::boost::logging::simple_logger_keeper keep = ::boost::logging::simple_logger_keeper(log_name,::boost::logging::level::lvl) ) ; else (*keep.stream())

#define BOOST_SCOPEDLOG(log_name) BOOST_SCOPEDLOGL(log_name,default_)

/* 
    Note: at runtime, enabling/disabling of compile-time loggers is COMPLETELY ignored!
*/

#define BOOST_DECLARE_LOG_COMPILE_TIME(log_name, compile_constant) \
    struct log_name ## _log_class { enum { is_compile_time = 1, is_enabled = ((compile_constant) != 0) }; }; \
    ::boost::logging::logger & log_name();

// helpers - logging is enabled only in debug-mode
#define BOOST_DECLARE_LOG_DEBUG(log_name) BOOST_DECLARE_LOG_COMPILE_TIME(log_name, BOOST_LOG_IS_DEBUG_MODE)





/*

          log statistics.

  FIXME
    - make it possible to EXTREMELY EASY: derive one log from another. like:
        log(dbg_gui) = log(dbg) + "gui";
        If dbg is turned off, dbg_gui is turned off.

*/



/*
Note:
You can create "copy" logs. That is, create log_id (A) that is a "shadow" of another log_id (B).
Any operation (enable/disable/add log func/etc.) that happens to B will reflect into A.

Just make A's name a sub-name of B 
(example: "app.gui.charts.dbg" (A) and "app.gui" (B) )
*/

#endif
