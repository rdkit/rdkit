
#include "logs.h"
#include <boost/log/functions.hpp>
#include <boost/log/extra/functions_ts.hpp>
#include <iostream>
#include <string>

BOOST_DEFINE_LOG(app, "app")
BOOST_DEFINE_LOG(dbg, "app.dbg")
BOOST_DEFINE_LOG(err, "app.err")


void init_logs() {
    using namespace boost::logging;
    // All Logs:
    // thread-safe write to file "out.txt", every 100 ms - on a dedicated thread
    add_appender("*", ts_appender( write_to_file("out.txt") ) );

    // App log only:
    // thread-safe write to console, each half a second - on a dedicated thread
    add_appender("app", ts_appender(write_to_cout, 500) );

    // Dbg & Err logs: write to Output Debug Window
#ifdef BOOST_LOG_WIN32
    add_appender("app.*", write_to_dbg_wnd);
#endif

    add_modifier("*", prepend_prefix);
    add_modifier("*", prepend_time( "$hh:$mm:$ss "), "", 200);

    flush_log_cache();
}
