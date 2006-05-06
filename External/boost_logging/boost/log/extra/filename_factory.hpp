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

#ifndef CE_BOOST_LOG_FILENAME_FACTORY_H_INCL
#define CE_BOOST_LOG_FILENAME_FACTORY_H_INCL

#include <ctime>
#include <boost/log/log_impl.hpp>

namespace boost { namespace logging {

/// filename factory method which passes the prototype through
/// "strftime"
inline logging_types::log_name_string_type
filename_factory_strftime (const logging_types::log_name_string_type& proto)
{
    char buf[FILENAME_MAX];

    time_t t = time (0);
    tm tm = *localtime (&t);
    int len = strftime (buf, sizeof (buf), proto.c_str (), &tm);

    return (len == 0 ? proto : logging_types::log_name_string_type (buf, len));
}

} // namespace log
} // namespace boost

#endif // macro guard

