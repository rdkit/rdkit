// Boost Logging Template library
//
// Author: John Torjo
//
// Copyright (C) 2004-2005 Caleb Epstein caleb.epstein@gmail.com
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

#ifndef CE_BOOST_LOG_FILE_APPENDER_HPP
#define CE_BOOST_LOG_FILE_APPENDER_HPP

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <boost/log/log_impl.hpp>
#include <boost/shared_ptr.hpp>

namespace boost { namespace logging {

/**
Logfile class encapsulates our file stream, filename prototype,
filename factory method and the actual generated filename
*/
struct logfile {
    typedef logging_types::log_name_string_type log_name_type;
    typedef log_name_type (*filename_factory) (const log_name_type& proto);

    filename_factory filename_generator;
    log_name_type filename_proto;
    log_name_type filename;
    std::ofstream file;
    
    logfile (const log_name_type& filename_proto,
             filename_factory factory = 0)
        : filename_proto (filename_proto),
          filename_generator (factory)
        {}

    virtual ~logfile ()                 { close (); }

    void close ()                       { file.close (); }

    void open (bool append = false)
        {
            file.clear ();

            if (filename_generator)
                filename = (*filename_generator) (filename_proto);
            else
                filename = filename_proto;

            file.open (filename.c_str (), 
                       std::ios::out |
                       (append ? std::ios::app : std::ios::trunc));
        }
};

/// Basic logging class that appends to a file specified by the user.
/// The filename can be a prototype filename which is
struct file_appender {
    typedef logging_types::log_name_string_type log_name_type;

    boost::shared_ptr<logfile> log;

    file_appender (const log_name_type& filename_proto,
                   logfile::filename_factory filename_factory = 0,
                   bool append = true)
        : log (new logfile (filename_proto, filename_factory))
        {
            log->open (append);
        }

    virtual ~file_appender () {}

    void operator () (const logging_types::string&, const logging_types::string& msg)
        { write (msg); }

protected:
    virtual void write (const logging_types::string& msg)
    { log->file.write (msg.data (), (std::streamsize)msg.size ()); }
};

} // namespace log
} // namespace boost

#endif  // macro guard


