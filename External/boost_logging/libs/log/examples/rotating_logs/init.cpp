// init.cpp - initialize the logs
//            (where do they output information, and how)

#include <boost/log/functions.hpp>
#include <boost/log/extra/appenders.hpp>

#include <iostream>

void init_logs() {
    using namespace boost::logging;
    // all logs prefix the message by time
    add_modifier("*", prepend_time("$hh:$mm:$ss "), "", DEFAULT_INDEX + 1 );
    // all log' messages are prefixed by the log name ("app", or "DEBUG" or "Inf")
    add_modifier("*", &prepend_prefix);

    // all messages are written to cout
    add_appender("*", write_to_cout);

    add_appender("*", 
        rolling_file_appender(
            "out.txt",  // file prefix
            1024 * 10,  // when this size reached, go to next file
            10          // how many files?
            ) );

    flush_log_cache();
}
