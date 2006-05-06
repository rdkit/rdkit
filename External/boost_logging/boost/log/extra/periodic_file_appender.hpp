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

#ifndef CE_BOOST_LOG_PERIODIC_FILE_APPENDER_HPP
#define CE_BOOST_LOG_PERIODIC_FILE_APPENDER_HPP

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <ctime>

#include <boost/log/log_impl.hpp>
#include <boost/log/extra/file_appender.hpp>

namespace boost { namespace logging {

/**
@brief File appender that allows the user to specify a rollover frequency, in seconds.  

Logfiles will be closed and re-opened
every <roll_period> seconds.  It is the responsibility of the user
will to provide a <filename_factory> function which generates a
unique filename based on the rollover frequency.
*/
struct periodic_file_appender : public file_appender {
    int roll_period;            /// number of seconds we wait between rollovers
    std::time_t period_begin;   /// starting time of the current period

    /// \note we could set <period_begin> %= <roll_period>,
    /// but this ties the period start times to GMT which may
    /// not be what the user desires.  Should we allow the user to
    /// specify the period start time based on their own wall clock?
    periodic_file_appender (const log_name_type& proto,
                            int roll_period,
                            logfile::filename_factory filename_factory = 0,
                            bool append = true)
        : file_appender (proto, filename_factory, append),
          roll_period (roll_period),
          period_begin (std::time (0))
        {}

protected:
    void write (const logging_types::string& msg) {
        std::time_t now = std::time (0);

        /// If roll_period seconds has elapsed, we close the file and
        /// open a new one (hopefully with a new name, but thats up to
        /// the filename factory function)
        if (now - period_begin >= roll_period)
            {
            log->close ();
            log->open (false);
            period_begin = now;
            }

        file_appender::write (msg);
    }
};

/// roll over files every hour
periodic_file_appender
hourly_file_appender (const logfile::log_name_type& proto,
                      logfile::filename_factory filename_factory,
                      bool append = true)
{
    return periodic_file_appender (proto, 60 * 60,
                                   filename_factory, append);
}

/// roll over files once a day
periodic_file_appender
daily_file_appender (const logfile::log_name_type& proto,
                     logfile::filename_factory filename_factory,
                     bool append = true)
{
    return periodic_file_appender (proto, 24 * 60 * 60,
                                   filename_factory, append);
}

/// roll over files once a week
periodic_file_appender
weekly_file_appender (const logfile::log_name_type& proto,
                      logfile::filename_factory filename_factory,
                      bool append = true)
{
    return periodic_file_appender (proto, 7 * 24 * 60 * 60,
                                   filename_factory, append);
}

/// roll over files once a month; the period length here is rather
/// problematic
periodic_file_appender
monthly_file_appender (const logfile::log_name_type& proto,
                       logfile::filename_factory filename_factory,
                       bool append = true)
{
    return periodic_file_appender (proto, 31 * 24 * 60 * 60,
                                   filename_factory, append);
}

/// roll over files once a year
periodic_file_appender
yearly_file_appender (const logfile::log_name_type& proto,
                      logfile::filename_factory filename_factory,
                      bool append = true)
{
    return periodic_file_appender (proto, 365 * 24 * 60 * 60,
                                   filename_factory, append);
}

} // namespace log
} // namespace boost

#endif  // macro guard


