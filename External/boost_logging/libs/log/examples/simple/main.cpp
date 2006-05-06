// Sample for Boost Logging
//

#include <boost/log/log_impl.hpp>
#include <boost/log/functions.hpp>
#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/limits.hpp>


// in a header file
BOOST_DECLARE_LOG(app)
BOOST_DECLARE_LOG(dbg)
BOOST_DECLARE_LOG(err)
BOOST_DECLARE_LOG(warn)
BOOST_DECLARE_LOG(info)

// ... namespaces work too...
namespace charts {
    BOOST_DECLARE_LOG(dbg)
    BOOST_DECLARE_LOG(gui)
}

// in a source file
BOOST_DEFINE_LOG(app, "app")
BOOST_DEFINE_LOG(dbg, "app.dbg")
BOOST_DEFINE_LOG(err, "app.err")
BOOST_DEFINE_LOG(warn, "app.warn")
BOOST_DEFINE_LOG(info, "info")

BOOST_DEFINE_LOG(charts::dbg, "charts.dbg")
BOOST_DEFINE_LOG(charts::gui, "charts.gui")


int main() {
    using namespace boost::logging;

    BOOST_LOG(app) << "pre-start" << " - this gets logged too!";
    BOOST_LOG(dbg) << "debug pre-start" << " - this gets logged too!";

    // Modifiers for all: 
    // [type_of_message] original_message append_enter_if_needed
    add_modifier("*", prepend_time("$hh:$mm:$ss "), "", INT_MAX );
    add_modifier("*", &append_enter);
    // Modifiers for app and its descendants
    // <time> [type_of_message] original_message append_enter_if_needed
    add_modifier("app*", &prepend_prefix);
    // Modifiers for "app" only
    // <time> [Thread ID] [type_of_message] original_message append_enter_if_needed
    add_modifier("app", &prepend_thread_id, 0);

    // Appenders
    // all messages are written to cout
    add_appender("*", write_to_cout);
    // "app*" messages are written to file as well
    add_appender("app*", write_to_file("out.txt") );
    // 'app' only and dbg messages are written to Output Debug Window as well
#ifdef BOOST_LOG_WIN32
    add_appender("app", write_to_dbg_wnd);
    add_appender("*dbg*", write_to_dbg_wnd);
#endif

    // signal that we've set the log properties (like, appenders and/or modifiers),
    // and we can start logging... - thus, any messages that got logged before now,
    // will be written to the appropriate destinations...
    flush_log_cache();

    int i = 1, j = 2, k = 3;
    BOOST_LOG(app) << "testing" << i << '-' << j << '-' << k;
    // written to both cout & the file & Output Debug Window
    BOOST_LOG(dbg) << "this is a debug message, i=" << i << std::endl;

    BOOST_LOG(info) << "I just wanted to tell you something....";
    BOOST_LOG(warn) << "Logged-on Users approach max. limit";
    BOOST_LOG(err) << "Too many users!";

    BOOST_LOG(charts::gui) << "Creating main window" << std::endl;
    BOOST_LOG(charts::dbg) << "A debug msg coming from {charts} module" ;
    BOOST_LOG(charts::gui) << "Destroying main window" << std::endl;

    // disable all descendants of 'app' (not the 'app' itself)
    disable_logs("app.*");
    BOOST_LOG(dbg) << "this won't be written" << std::endl;
    BOOST_LOG(app) << "However, this msg. will" << std::endl;
    enable_logs("app.dbg"); // specifically, only dbg log is enabled back now
    BOOST_LOG(dbg) << "this will be written - this log just got enabled" << std::endl;
    BOOST_LOG(err) << "this still won't be written" << std::endl;
    enable_logs("app.*");
    disable_logs("app.dbg");
    // now, all logs are back to the 'enabled' state
    BOOST_LOG(err) << "this will be written - all logs are enabled" << std::endl;

    // disable all logs
    disable_logs("*");
    BOOST_LOG(err) << "this won't be written" << std::endl;
    BOOST_LOG(app) << "neither will this" << std::endl;
    BOOST_LOG(dbg) << "or this..." << std::endl;
    BOOST_LOG(warn) << "or this..." << std::endl;
    BOOST_LOG(info) << "or this..." << std::endl;

    // enable all dbg logs
    enable_logs("*dbg*");
    BOOST_LOG(app) << "this won't be written" << std::endl;
    BOOST_LOG(dbg) << "this will be written" << std::endl;
    BOOST_LOG(info) << "this won't be written" << std::endl;
    // enable info log
    enable_logs("*info*");
    BOOST_LOG(info) << "a simple info" << std::endl;

    /*** Console

        12:52:18 [app] [Thread 7528] pre-start - this gets logged too!
        12:52:18 [app.dbg] debug pre-start - this gets logged too!
        12:52:18 [app] [Thread 7528] testing1-2-3
        12:52:18 [app.dbg] this is a debug message, i=1
        12:52:18 I just wanted to tell you something....
        12:52:18 [app.warn] Logged-on Users approach max. limit
        12:52:18 [app.err] Too many users!
        12:52:18 Creating main window
        12:52:18 A debug msg coming from {charts} module
        12:52:18 Destroying main window
        12:52:18 [app] [Thread 7528] However, this msg. will
        12:52:18 [app.dbg] this will be written - this log just got enabled
        12:52:18 [app.err] this will be written - all logs are enabled
        12:52:18 [app.dbg] this will be written
        12:52:18 a simple info
    */

    /*** Output Debug Window

        12:52:18 [app] [Thread 7528] pre-start - this gets logged too!
        12:52:18 [app.dbg] debug pre-start - this gets logged too!
        12:52:18 [app] [Thread 7528] testing1-2-3
        12:52:18 [app.dbg] this is a debug message, i=1
        12:52:18 A debug msg coming from {charts} module
        12:52:18 [app] [Thread 7528] However, this msg. will
        12:52:18 [app.dbg] this will be written - this log just got enabled
        12:52:18 [app.dbg] this will be written
    */

    /*** out.txt file

        12:52:18 [app] [Thread 7528] pre-start - this gets logged too!
        12:52:18 [app.dbg] debug pre-start - this gets logged too!
        12:52:18 [app] [Thread 7528] testing1-2-3
        12:52:18 [app.dbg] this is a debug message, i=1
        12:52:18 [app.warn] Logged-on Users approach max. limit
        12:52:18 [app.err] Too many users!
        12:52:18 [app] [Thread 7528] However, this msg. will
        12:52:18 [app.dbg] this will be written - this log just got enabled
        12:52:18 [app.err] this will be written - all logs are enabled
        12:52:18 [app.dbg] this will be written

    */

}



