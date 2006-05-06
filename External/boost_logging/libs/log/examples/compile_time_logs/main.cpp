// Sample for Boost Logging
//

/** 

    Shows how you can define logs that can be enabled/disabled at COMPILE-TIME

    The only thing that differs, is how the logs are declared. The rest (defining
    and using them in code) is THE SAME.
*/

#include <boost/log/log_impl.hpp>
#include <boost/log/functions.hpp>
#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/limits.hpp>


// can be enabled/disabled at runtime
BOOST_DECLARE_LOG(app)

// note: this log is enabled only in debug mode
BOOST_DECLARE_LOG_DEBUG(dbg)

// note: this log is enabled only in release mode
BOOST_DECLARE_LOG_COMPILE_TIME(rel, !BOOST_LOG_IS_DEBUG_MODE)

// in a source file
BOOST_DEFINE_LOG(app, "app")
BOOST_DEFINE_LOG(dbg, "debug")
BOOST_DEFINE_LOG(rel, "release")



int main() {
    using namespace boost::logging;

    BOOST_LOG(app) << "pre-start" << " - this gets logged too!";
    BOOST_LOG(dbg) << "debug pre-start" << " - this gets logged too!";

    // Modifiers for all: 
    // [type_of_message] original_message append_enter_if_needed
    manipulate_logs("*")
        .add_modifier( prepend_time("$hh:$mm:$ss "), INT_MAX )
        .add_modifier( &append_enter)
        .add_modifier( &prepend_prefix)
        .add_appender( &write_to_cout);

    // "app" messages are written to file as well
    manipulate_logs("app").add_appender( write_to_file("out.txt") );
    // 'app' only and dbg messages are written to Output Debug Window as well
#ifdef BOOST_LOG_WIN32
    manipulate_logs("app").add_appender( write_to_dbg_wnd);
    manipulate_logs("debug").add_appender( write_to_dbg_wnd);
#endif
    // signal that we've set the log properties (like, appenders and/or modifiers),
    // and we can start logging... - thus, any messages that got logged before now,
    // will be written to the appropriate destinations...
    flush_log_cache();

    disable_logs("app");
    BOOST_LOG(app) << "this WILL NOT be written ANYWHERE";
    enable_logs("app");

    // note: enabling/disabling is completely ignored, for compile-time logs
    disable_logs("debug");
    disable_logs("release");

    int i = 1, j = 2, k = 3;
    BOOST_LOG(app) << "testing" << i << '-' << j << '-' << k;
    // written to both cout & the file & Output Debug Window
    BOOST_LOG(dbg) << "this is a debug message, i=" << i << std::endl;

    BOOST_LOG(rel) << "I just wanted to tell you something....";

    BOOST_LOG(dbg) << "Creating main window" << std::endl;
    BOOST_LOG(dbg) << "One more debug message" ;
    BOOST_LOG(dbg) << "Destroying main window" << std::endl;

    // now, all logs are back to the 'enabled' state
    BOOST_LOG(rel) << "this appears only in release mode, i =" << i << std::endl;
    if ( BOOST_IS_LOG_ENABLED(rel, default_) )
        BOOST_LOG(rel) << "this appears only in release mode, k =" << k << std::endl;

    BOOST_LOG(app) << "the end" << std::endl;
    std::cin.get();
}



