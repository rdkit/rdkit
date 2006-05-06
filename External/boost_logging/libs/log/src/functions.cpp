// functions.cpp 

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
    Note: I supplied this as a static library, because otherwise 
    compilation times would be too big.
*/

//FIXME update docs about new functionality + log manager. Tell about the macros , and that you can say:
// logger l("app.gui");
// if ( enabled_logger ls = l) ls << whatever << you << want;

// I should also allow for l.stream() << whatever << you << want;

#include <boost/log/functions.hpp>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <iostream>
#include <boost/bind.hpp>

#ifdef BOOST_LOG_WIN32
#include <windows.h>
#endif

namespace boost { namespace logging {

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Log functions

// ... Windows only 
#ifdef BOOST_LOG_WIN32
    BOOST_LOG_DECL void write_to_dbg_wnd(const logging_types::string &, const logging_types::string &msg) {
        OutputDebugString( msg.c_str() );
    }
#endif


void write_to_file::write(const logging_types::string& msg) {
//    std::ios::openmode open_type = std::ios::out | ((*overwrite_file) ? 0 : std::ios::app);
    // many thanks to Wesselink!
    std::ios::openmode open_type = (*overwrite_file) ? (std::ios::out) : (std::ios::out | std::ios::app);

    std::basic_ofstream<logging_types::char_t> out( file_name.c_str(), open_type);
    out << msg;

    *overwrite_file = false;
}

    
    
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Modifier functions



// writes the prefix
BOOST_LOG_DECL void prepend_prefix(const logging_types::string &prefix, logging_types::string & msg) {
    msg = '[' + prefix + "] " + msg;
}

// appends an enter if not already found
BOOST_LOG_DECL void append_enter(const logging_types::string &, logging_types::string & msg) {
    if ( msg.empty() ) return;
    if ( *msg.rbegin() != '\n') msg += '\n';
}



BOOST_LOG_DECL void prepend_thread_id(const logging_types::string &, logging_types::string & msg) {
    std::ostringstream out;
    out << "[Thread "
#if defined (BOOST_HAS_WINTHREADS)
        << ::GetCurrentThreadId()
#elif defined (BOOST_HAS_PTHREADS)
        << pthread_self ()
#elif defined (BOOST_HAS_MPTASKS)
        << MPCurrentTaskID()
#endif
        << "] ";
    msg = out.str() + msg;
}




namespace {
    struct index_info {
        index_info(int src_idx, int *format_idx, int size = 2) : src_idx(src_idx), format_idx(format_idx), size(size) {}
        int src_idx;
        int * format_idx;
        int size;
    };
    bool by_index(const index_info & first, const index_info & second) {
        return first.src_idx < second.src_idx;
    }
}

prepend_time::prepend_time(const logging_types::string & format) : m_day(-1), m_month(-1), m_yy(-1), m_yyyy(-1), m_hour(-1), m_min(-1), m_sec(-1) {
    // format too big
    assert( format.size() < 256);

    typedef logging_types::string::size_type uint;
    int day_idx    = (int)format.find("$dd");
    int month_idx  = (int)format.find("$MM");
    int yy_idx     = (int)format.find("$yy");
    int yyyy_idx   = (int)format.find("$yyyy");
    int hour_idx   = (int)format.find("$hh");
    int min_idx    = (int)format.find("$mm");
    int sec_idx    = (int)format.find("$ss");

    typedef std::vector<index_info> array;
    array indexes;
    if ( day_idx != logging_types::string::npos)
        indexes.push_back( index_info(day_idx, &m_day) );
    if ( month_idx != logging_types::string::npos)
        indexes.push_back( index_info(month_idx, &m_month) );

    if ( yy_idx != logging_types::string::npos || yyyy_idx != logging_types::string::npos)
        if ( yyyy_idx  != logging_types::string::npos)
            indexes.push_back( index_info(yyyy_idx, &m_yyyy, 4) );
        else
            indexes.push_back( index_info(yy_idx, &m_yy) );

    if ( hour_idx != logging_types::string::npos)
        indexes.push_back( index_info(hour_idx, &m_hour ) );
    if ( min_idx != logging_types::string::npos)
        indexes.push_back( index_info(min_idx, &m_min) );
    if ( sec_idx != logging_types::string::npos)
        indexes.push_back( index_info(sec_idx, &m_sec) );
    std::sort( indexes.begin(), indexes.end(), by_index);
    
    // create the format logging_types::string, that we can actually pass to sprintf 
    int prev_idx = 0;
    int idx = 0;
    for ( array::iterator begin = indexes.begin(), end = indexes.end(); begin != end; ++begin) {
        m_format += format.substr( prev_idx, begin->src_idx - prev_idx);
        *begin->format_idx = idx;
        m_format += (begin->size == 4) ? "%04d" : "%02d";
        prev_idx = begin->src_idx + begin->size + 1;
        ++idx;
    }

    m_format += format.substr(prev_idx);
}

void prepend_time::operator()(const logging_types::string &, logging_types::string & msg) {
    time_t t = time(0); 
    tm details = *localtime( &t);

    int vals[8];
    vals[m_day + 1]      = details.tm_mday;
    vals[m_month + 1]    = details.tm_mon + 1; // many thanks to Matthew P. Cashdollar
    vals[m_yy + 1]       = details.tm_year % 105; // many thanks to Andy Schweitzer
    vals[m_yyyy + 1]     = details.tm_year + 1900;
    vals[m_hour + 1]     = details.tm_hour;
    vals[m_min + 1]      = details.tm_min;
    vals[m_sec + 1]      = details.tm_sec;
  
    // ignore value at index 0 - it's there so that I don't have to test for an index being -1
    sprintf( m_buffer, m_format.c_str(), vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7] );
    msg = m_buffer + msg;
}



prepend_time_strf::prepend_time_strf(const logging_types::string & format, bool localtime)
    : m_format (format), m_localtime (localtime)
{}

void prepend_time_strf::operator() (const logging_types::string&, logging_types::string& msg) {
    m_t = time (0); 
    m_tm = m_localtime ? *localtime( &m_t) : *gmtime( &m_t);
    if (0 != strftime (m_buffer, sizeof (m_buffer), m_format.c_str (), &m_tm))
        msg = m_buffer + msg;
}


void appender_array::operator()( const logging_types::string & log_name, const logging_types::string & msg) {
    for( std::vector<appender_func>::const_iterator begin = m_funcs.begin(), end = m_funcs.end(); begin != end; ++begin)
        (*begin)(log_name, msg);
}


void write_to_cout(const std::string &, const std::string &msg) { std::cout << msg; }

}}

