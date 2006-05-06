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

#ifndef CE_BOOST_LOG_ROLLING_FILE_APPENDER_HPP
#define CE_BOOST_LOG_ROLLING_FILE_APPENDER_HPP

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/log/log_impl.hpp>
#include <boost/log/extra/file_appender.hpp>
#include <boost/filesystem/operations.hpp>

namespace boost { namespace logging {


/**
@brief File appender that allows user to specify a maximum file size and number of files to keep.  

When logfiles meet or exceed the maximum
file size, they are renamed to numbered backups in the following
progression: file -> file.1 -> file.2 -> ... -> file.<num_files>
-> removal.  Therefore the most recent logfile is always <file>
and the nexst oldest is file.1 etc.
*/
struct rolling_file_appender : public file_appender {
    int max_size;
    int num_files;

    rolling_file_appender (const log_name_type& proto,
                           int max_size,
                           int num_files = 0,
                           logfile::filename_factory filename_factory = 0,
                           bool append = true)
        : file_appender (proto, filename_factory, append),
          max_size (max_size),
          num_files (num_files)
        {}

protected:
    void write (const logging_types::string& msg) {
        file_appender::write (msg);

        /// If there is a file size limit and we have met or exceeded
        /// it, close the log, run the backup policy and open the
        /// logfile again
        if (max_size != 0 && log->file.tellp () >= max_size)
            {
            log->close ();
            rollover (log->filename);
            log->open (false);
            }
    }

    /// roll over the logfile
    virtual void rollover (const log_name_type& file) const {
        using namespace boost::filesystem;

        for (int suffix = num_files - 1; suffix >= 0; --suffix) {
            std::ostringstream from, to;
            from << file;
            if (suffix)
                from << '.' << suffix;
            to << file << '.' << (suffix + 1);
            path from_p (from.str ()), to_p (to.str ());
            if (exists (from_p)) {
                remove (to_p);
                rename (from_p, to_p);
            }
        }
    }
};


}} 

#endif  // macro guard


