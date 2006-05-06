// functions.hpp

// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 John Torjo (john@torjo.com)
// Copyright (C) 2004-2005 Caleb Epstein (caleb.epstein@gmail.com)
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
    Helper Logging functions

    Note: including this file brings a boost_thread dependency
*/


#ifndef JT_BOOST_appender_funcTIONS_ts_HPP
#define JT_BOOST_appender_funcTIONS_ts_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/log/log_impl.hpp>

#if defined(BOOST_HAS_THREADS)
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#ifndef BOOST_LOG_WIN32
// for do_sleep
#include <boost/thread/xtime.hpp>
#endif

namespace boost { namespace logging {

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Log functions (thread-safe)



/*
    Thread-safe appender

    Caches all messages written to it, and writes them, once at X milliseconds ON A DEDICATED THREAD.
    (writes them = calls the functor you passed at constructor with the cached messages that have been concatenated)

    If you want to forward to multiple appenders you can use an appender_array as the underlying appender.


    Note: normally, the first parameter passed to the appender function is the log's name.
    In our case, when we cache multiple messages, it would be way too complex to remember the log name for each message
    (in order to forward them, to the underlying appender).
    In this implementation, we will always pass an empty string ("") as the first argument to the underlying appender.
*/    
struct BOOST_LOG_DECL ts_appender {
    ts_appender( appender_func func, int sleep_ms = 100) : m_info(new info) {
        m_info->func = func;
        m_info->sleep_ms = sleep_ms;
    }
    ~ts_appender() {
        if ( !m_info.unique() )
            return; // we'll terminate the thread only when we're the only thread left

        // note: even if we're the only object left, we still need to use locks
        // there could be another thread: the log_thread; so we still need to be thread-safe
        int sleep_ms;
        { logging_types::lock lk( m_info->cs);
        if ( m_info->h_log == 0) 
            return; /* this thread was not even started... */ 

        sleep_ms = m_info->sleep_ms;
        m_info->terminated = true; }

        // FIXME in the future, when thread library settles down, I could use conditions

        // wait for other thread to exit - note: join() might be too much. It could end up waiting forever...
        do_sleep( (int)(m_info->sleep_ms * 1.5) );
    //    m_info->h_log->join ();
    }

    void operator()( const logging_types::string &, const logging_types::string & msg) {
        bool create_now = false;
        // FIXME maybe I can create the thread in the constructor?
        { logging_types::lock lk( m_info->cs);
        m_info->str += msg;
        if ( m_info->h_log == 0) create_now = true;
        }

        if (create_now)
            create_thread();
    }

private:
    struct info {
        info() : terminated(false), sleep_ms(0), h_log(0) {
        }
        logging_types::mutex cs;

        appender_func func;
        bool terminated;
        logging_types::string str;
        int sleep_ms;

        // FIXME leak - use boost::shared_ptr!!!!
        boost::thread * h_log;
    };

    static void log_thread(info * self) {
        logging_types::string msg;
        bool terminated = false;
        int sleep_ms;

        while ( !terminated) {
            { logging_types::lock lk( self->cs);
            msg = self->str;
            self->str.erase(); 
            sleep_ms = self->sleep_ms; }
            if ( !msg.empty() )
                try {
                    self->func("", msg);
                }
                catch (...) {}

            { logging_types::lock lk( self->cs);
            terminated = self->terminated; }
            if ( !terminated)
                do_sleep( sleep_ms);
        }
    }


    static void do_sleep(int ms) {
#ifdef BOOST_LOG_WIN32
        ::Sleep(ms);
#else
        using namespace boost;
        xtime next;
        xtime_get(&next, TIME_UTC);
        next.sec += ms / 1000;
        next.nsec += (ms % 1000) * 1000;
        thread::sleep(next);
#endif    
    }

    void create_thread() {
        {
        logging_types::lock lk( m_info->cs);
        if ( m_info->h_log != 0) return; // already created

        using namespace boost;
        m_info->h_log = new thread( bind( ts_appender::log_thread,m_info.get()) );
        }
        // make sure the thread gets created and gets a chance to access the info
        do_sleep(100); 
    }


private:
    // every data should be shared. The log_thread needs all this info.
    //
    // log_thread should not get a reference to this, because objects of type 'ts_appender' could get moved around.
    //
    // For example, every time you add a Log Func (add_appender_func), the log functions might get rearranged,
    // thus existing 'ts_appender' objects could get moved around, so a 'this' pointer could later become invalid!
    //
    boost::shared_ptr<info> m_info;
};



}}


#endif // BOOST_HAS_THREADS

#endif

