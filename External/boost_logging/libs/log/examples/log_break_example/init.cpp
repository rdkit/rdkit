// init.cpp - initialize the logs
//            (where do they output information, and how)

#include <boost/log/functions.hpp>

#include <log_break/log_break.hpp>

#include <iostream>

void init_logs() {
    using namespace boost::logging;
    // all logs prefix the message by time
    add_modifier("*", prepend_time("$hh:$mm:$ss "), "", DEFAULT_INDEX + 1 );
    // all log' messages are prefixed by the log name 
    add_modifier("*", &prepend_prefix);

    // VERY IMPORTANT: this is how we monitor all messages. So that we can break into debug mode,
    // if necessary!
    add_appender("*", log_break("out.txt") );

    // all messages are written to a file
    add_appender("*", write_to_file("out.txt") );
    // all messages are written to cout 
    // note: this is just for brevity, so you can easier understand what's happening
    add_appender("*", write_to_cout);
    // dbg & warn messages are also written to Output Debug Window
    // note: this is just for brevity, to see how easy it is to write to multiple locations
#ifdef BOOST_LOG_WIN32
    add_appender("app.*", write_to_dbg_wnd );
#endif

    flush_log_cache();
}
