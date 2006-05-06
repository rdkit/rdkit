// log.cpp 

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
    Note: I supplied this as a static library, because otherwise 
    compilation times would be too big.
*/

#include <boost/log/log_impl.hpp>
#include <map>
#include <vector>
#include <boost/detail/atomic_count.hpp>
#include <algorithm>
#include <iostream>

namespace boost { namespace logging { 

namespace {
    // .. used in default_log_manager_impl
    struct logger_info {
        logger_info( boost::shared_ptr<logger_impl> impl) : impl(impl) {}
        boost::shared_ptr<logger_impl> impl;

        struct unwritten_msg {
            unwritten_msg(const logging_types::log_name_string_type & str, level_type lvl) : str(str), lvl(lvl) {}
            logging_types::log_name_string_type str;
            level_type lvl;
        };


        typedef std::vector<unwritten_msg> unwritten_array;
        unwritten_array unwritten;
    };
    typedef std::map<logging_types::log_name_string_type, logger_info > logger_coll;
    
    struct appender_info {
        appender_info(const logging_types::log_name_string_type & spec, appender_func func, const logging_types::log_name_string_type & name, int idx) : spec(spec), func(func), name(name), idx(idx) {}
        logging_types::log_name_string_type spec;
        // note: we should not keep a shared pointer
        //       each logger_impl has its own copy one appender func - this is totally ok
        appender_func func;
        logging_types::log_name_string_type name;
        int idx;
    };
    typedef std::vector<appender_info> appender_array;

    struct modifier_info {
        modifier_info(const logging_types::log_name_string_type & spec, modifier_func func, const logging_types::log_name_string_type & name, int idx) : spec(spec), func(func), name(name), idx(idx) {}
        logging_types::log_name_string_type spec;
        // note: we should not keep a shared pointer
        //       each logger_impl has its own copy one modifier func - this is totally ok
        modifier_func func;
        logging_types::log_name_string_type name;
        int idx;
    };
    typedef std::vector<modifier_info> modifier_array;

    struct enabled_info {
        enabled_info(const logging_types::log_name_string_type & spec, level_type lvl) : spec(spec), lvl(lvl) {}
        logging_types::log_name_string_type spec;
        level_type lvl;
    };
    typedef std::vector<enabled_info> enabled_array;

    // by default, if trying to write a message, and an exception occurs, this is called
    void on_exception_write_to_cout(const logging_types::string & msg, level_type lvl) {
        std::cout << msg;
    }

    struct default_log_manager_impl {
        default_log_manager_impl() : m_caching(true), m_cache_limit(DEFAULT_CACHE_LIMIT) {
            on_exception_writer() = &on_exception_write_to_cout;
        }
        ~default_log_manager_impl() {
            // just in case...
            flush_log_cache();
            s_destroyed = true;
        }

        // should never be deleted (it will be deleted by the OS itself, when the app. ends).
        // we need this to outlive everything, because we need - before every write,
        // to check out if you can write to the log
        static logging_types::mutex & s_cs() {
            static logging_types::mutex* p = new logging_types::mutex;
            return *p;
        }

        // should never be deleted (it will be deleted by the OS itself, when the app. ends).
        // we need this to outlive everything, because we never know when it'll be called
        // (after manager_impl's destruction, for instance)
        static on_exception_writer_func & on_exception_writer() {
            static on_exception_writer_func * p = new on_exception_writer_func ;
            return *p;
        }

        logger_coll m_logs;
        appender_array m_appenders;
        modifier_array m_modifiers;
        enabled_array m_enabled;

        // if true, we're caching messages written to loggers
        // (in other words, some of the loggers/loggers properties have not been set yet
        //  like: appenders/ modifiers)
        bool m_caching;

        // if true, we've been destroyed...
        static bool s_destroyed;

        // in case caching, this is the max. number of messages to cache, per one log
        int m_cache_limit;
    };
    bool default_log_manager_impl::s_destroyed = false;
//    logging_types::mutex* default_log_manager_impl::s_cs = new logging_types::mutex;


    default_log_manager_impl & manager_impl() {
        static default_log_manager_impl impl;
        return impl;
    }

    // returns true if the string matches the specification
    // ... empty string and '*' match everything
    bool name_matches_spec(const logging_types::log_name_string_type & log_name, const logging_types::log_name_string_type & spec) {
        if ( (spec == "*") || spec.empty() ) 
            return true; // matches all

        std::vector<logging_types::log_name_string_type> searches;
        logging_types::log_name_string_type::size_type pos = 0, next;
        while ( pos < spec.size() ) {
            next = spec.find('*', pos);
            if ( next == logging_types::log_name_string_type::npos) next = spec.size();
            if ( next - pos > 0)
                searches.push_back( spec.substr(pos, next - pos) );
            pos = next + 1;
        }
        if ( searches.empty())
            return true; // string was like "****" (only wildcards)

        // if fixed, when search starts/ends, the log_name must start/end
        bool fixed_start = *spec.begin() != '*';
        bool fixed_end = *spec.rbegin() != '*';

        logging_types::log_name_string_type::size_type at = 0;
        for ( std::vector<logging_types::log_name_string_type>::const_iterator begin = searches.begin(), end = searches.end(); begin != end ; ++begin)  {
            logging_types::log_name_string_type::size_type next = log_name.find(*begin, at);
            if ( next == logging_types::log_name_string_type::npos)
                return false; // doesn't match
            at = next + begin->length();
        }

        if ( fixed_start)
            if ( log_name.find(*searches.begin()) != 0)
                return false; // log_name should be like "word*" (start with 'word')
        if ( fixed_end)
            if ( at != log_name.length() )
                return false; // log_name should be like "*word" (end in 'word')

        return true;
    }


    // writes this message to the log - even if caching is on.
    //
    // caller must hold the (log_manager) lock!
    void force_write_msg(logger_impl & l, const logging_types::string & msg, level_type lvl) {
        if ( manager_impl().s_destroyed) 
            return; // we've been destroyed...

        try {
            l.write_msg(msg, lvl);
        } catch(...) {
            // error while writing...
            try {
                manager_impl().on_exception_writer()(msg, lvl);
            }
            catch(...) {
                // the on_exception_writer yielded an exception as well
                // there's nothing we can do.
            }
        }
    }

    void cache_msg(const logger & l, const logging_types::log_name_string_type & msg, level_type lvl) {
        if ( !l.is_enabled(lvl) )
            return; // log is disabled

        bool needs_flush = false;
        try {
            logging_types::lock lk(manager_impl().s_cs());

            logger_coll::iterator found = manager_impl().m_logs.find( l.name() );
            if ( found != manager_impl().m_logs.end() ) {
                found->second.unwritten.push_back( logger_info::unwritten_msg(msg,lvl) );
                needs_flush = ((int)found->second.unwritten.size() > manager_impl().m_cache_limit);
            }
            else
                assert(false); // this log should already belong to manager_impl()
        } catch(...) {
            // the only exception that could occur is on push_back (allocating memory).
            // nothing to do.
        }

        // if cache limit reached, caching is turned off.
        if ( needs_flush)
            flush_log_cache();
    }

}



BOOST_LOG_DECL void add_appender(default_log_manager & , const logging_types::log_name_string_type & logs_spec, appender_func f, const logging_types::log_name_string_type & name, int idx ) {
    logging_types::lock lk(manager_impl().s_cs());
    // every operation reflects to all logs that match logs_spec (even logs that will be created in the future)!

    appender_info info( logs_spec, f, name, idx);
    manager_impl().m_appenders.push_back( info);

    for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin)
        if ( name_matches_spec( begin->second.impl->name, logs_spec) )
            begin->second.impl->add_appender(f, name, idx);
}

// adds a function that modifies messages (like, adds a prefix and/or suffix)
BOOST_LOG_DECL void add_modifier(default_log_manager & , const logging_types::log_name_string_type & logs_spec, modifier_func f, const logging_types::log_name_string_type & name, int idx ) {
    logging_types::lock lk(manager_impl().s_cs());
    // every operation reflects to all logs that match logs_spec (even logs that will be created in the future)!

    modifier_info info( logs_spec, f, name, idx);
    manager_impl().m_modifiers.push_back( info);

    for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin)
        if ( name_matches_spec( begin->second.impl->name, logs_spec) )
            begin->second.impl->add_modifier(f, name, idx);
}

namespace {
    struct has_spec_and_name {
        has_spec_and_name (const logging_types::log_name_string_type & spec, const logging_types::log_name_string_type & name) : spec(spec), name(name) {}
        template<class modifier_or_appender> bool operator()(const modifier_or_appender & val) {
            return (val.spec == spec) && (val.name == name);
        }
        logging_types::log_name_string_type spec;
        logging_types::log_name_string_type name;
    };
}

BOOST_LOG_DECL void del_appender(default_log_manager & , const logging_types::string & spec, const logging_types::log_name_string_type & name) {
    logging_types::lock lk(manager_impl().s_cs());

    for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin)
        if ( name_matches_spec( begin->second.impl->name, spec) )
            begin->second.impl->del_appender(name);

    // existing appenders
    manager_impl().m_appenders.erase(
        std::remove_if( manager_impl().m_appenders.begin(), manager_impl().m_appenders.end(), has_spec_and_name(spec,name)),
        manager_impl().m_appenders.end() );
}

BOOST_LOG_DECL void del_modifier(default_log_manager & , const logging_types::string & spec, const logging_types::log_name_string_type & name) {
    logging_types::lock lk(manager_impl().s_cs());

    for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin)
        if ( name_matches_spec( begin->second.impl->name, spec) )
            begin->second.impl->del_modifier(name);

    // existing modifiers
    manager_impl().m_modifiers.erase(
        std::remove_if( manager_impl().m_modifiers.begin(), manager_impl().m_modifiers.end(), has_spec_and_name(spec,name)),
        manager_impl().m_modifiers.end() );
}


namespace {
    // refreshes each log's "enabled" property. Called after an "enable_logs" has happened
    void refresh_log_enabled() {
        logging_types::lock lk(manager_impl().s_cs());

        for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin) {
            level_type lvl = level::default_;
            for ( enabled_array::iterator b = manager_impl().m_enabled.begin(), e = manager_impl().m_enabled.end(); b != e; ++b)
                if ( name_matches_spec( begin->second.impl->name, b->spec) )
                    lvl = b->lvl;
            begin->second.impl->enable(lvl);
        }
    }
}

BOOST_LOG_DECL void enable_logs(default_log_manager & , const logging_types::log_name_string_type & logs_spec, level_type lvl) {
    { logging_types::lock lk(manager_impl().s_cs());
    // every operation reflects to all logs that match logs_spec (even logs that will be created in the future)!
    manager_impl().m_enabled.push_back( enabled_info(logs_spec, lvl) );
    // FIXME if the same spec is found again before, remove the before, but add it at the end
    //       if the same spec with !enabled is found before, remove it from there
    }
    refresh_log_enabled();
}


BOOST_LOG_DECL void write_msg(default_log_manager & , logger & l, const logging_types::log_name_string_type & msg, level_type lvl) {
    if ( manager_impl().s_destroyed)
        return;

    // we're caching - if loggers have not been initialized yet...
    if ( manager_impl().m_caching)
        cache_msg( l, msg, lvl);
    else {
        // double-test : we want to know if the logger been destroyed
        //               that is, when the app terminates, we can't control the destroy order.
        //               thus, BOOST_LOG could be called even after the logs have been destroyed
        if ( l.still_exists() )
            force_write_msg( logger_to_logger_impl(l), msg, lvl);
    }
}

BOOST_LOG_DECL boost::shared_ptr<logger_impl> find_log_by_name( default_log_manager & , const logging_types::log_name_string_type & log_name) {
    logging_types::lock lk(manager_impl().s_cs());
    logger_coll::iterator found = manager_impl().m_logs.find(log_name);
    if ( found != manager_impl().m_logs.end() )
        return found->second.impl;

    // this log does not yet exist, let's create it now
    boost::shared_ptr<logger_impl> new_log( new logger_impl(log_name) );

    // add appenders
    for ( appender_array::iterator begin = manager_impl().m_appenders.begin(), end = manager_impl().m_appenders.end(); begin != end; ++begin)
        if ( name_matches_spec( log_name, begin->spec))
            new_log->add_appender(begin->func, begin->name, begin->idx);

    // add modifiers
    for ( modifier_array::iterator begin = manager_impl().m_modifiers.begin(), end = manager_impl().m_modifiers.end(); begin != end; ++begin)
        if ( name_matches_spec( log_name, begin->spec))
            new_log->add_modifier(begin->func, begin->name, begin->idx);

    // find out its level
    level_type lvl = level::default_;
    for ( enabled_array::iterator begin = manager_impl().m_enabled.begin(), end = manager_impl().m_enabled.end(); begin != end; ++begin)
        if ( name_matches_spec( log_name, begin->spec))
            lvl = begin->lvl;
    new_log->enable( lvl);

    manager_impl().m_logs.insert( std::make_pair(log_name, new_log) );
    return new_log;
}

BOOST_LOG_DECL void mark_log_as_compile_time(const logging_types::log_name_string_type & log_name) {
    find_log_by_name( logger::manager(), log_name)->set_is_compile_time(true);
}



/*
    sets whether we do log caching and, if so, what's the cache limit. If not caching, the cache limit is ignored.

    If we're caching: messages sent to a log are not written to the log, but are cached.
    Once you set the caching to false, those cached messages are automatically written to their corresponding logs.

    Cache limit: how many messages we cache, for a given log. If this limit is reached, 
    the caching is automatically turned to 'false', and everything is flushed. This is so that you
    won't mistakenly forget to do caching...


    Usefulness of caching: at the beginning of the program, you might have to write 
    information to certain logs, before actually setting those logs' appenders & modifiers.
*/
BOOST_LOG_DECL void set_log_caching(bool do_cache, int cache_limit) {
    logging_types::lock lk(manager_impl().s_cs());
    manager_impl().m_caching = do_cache;
    if ( do_cache) {
        manager_impl().m_cache_limit = cache_limit;
    }
    else {
        if ( manager_impl().s_destroyed)
            return;

        // flush everything.
        for ( logger_coll::iterator begin = manager_impl().m_logs.begin(), end = manager_impl().m_logs.end(); begin != end; ++begin) {
            for ( logger_info::unwritten_array::iterator b = begin->second.unwritten.begin(), e = begin->second.unwritten.end(); b != e; ++b)
                force_write_msg( *begin->second.impl, b->str, b->lvl );
            logger_info::unwritten_array empty;
            std::swap( begin->second.unwritten, empty);
        }
    }
}

// writes everything that was is in cache - to their corresponding logs. Also, sets the log' caching to false
BOOST_LOG_DECL void flush_log_cache() {
    set_log_caching(false);
}


BOOST_LOG_DECL void on_exception_writer(default_log_manager & , on_exception_writer_func func) {
    manager_impl().on_exception_writer() = func;
}


/* 
    FIXME removing logs...
    It should be like this: when a logger goes out of scope, if its use_count() is 2
    (that is, this is the was the last reference to the logger, simply remove it. Maybe an option
    should exist to dictate this - something like, an auto_remove() property)
*/

}}


