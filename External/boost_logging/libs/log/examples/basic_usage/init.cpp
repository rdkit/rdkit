// init.cpp - initialize the logs
//            (where do they output information, and how)

#include <boost/log/functions.hpp>

#include <iostream>

void init_logs() {
    using namespace boost::logging;

    manipulate_logs("*")
        // all logs prefix the message by time
        .add_modifier(prepend_time("$hh:$mm:$ss "), DEFAULT_INDEX + 1 )
        // all log' messages are prefixed by the log name ("app", or "DEBUG" or "Inf")
        .add_modifier(&prepend_prefix)
        // all messages are written to cout
        .add_appender(&write_to_cout);

    // app messages are also written to file 'out.txt'
    manipulate_logs("app").add_appender(write_to_file("out.txt") );
    // dbg messages are also written to Output Debug Window
#ifdef BOOST_LOG_WIN32
    manipulate_logs("DEBUG").add_appender(write_to_dbg_wnd );
#endif
    flush_log_cache();
}
