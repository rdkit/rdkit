
/*
    Sample to show that using Boost Logging library incurs in almost no
    overhead in compilation-time.

    This file takes most to compile, because it manipulates the logs'
    behavior: add_modifier_func/add_appender/enable/disabling of logs.

    Usually, in your applications, this code will confine to one or two files.
    The rest of the files will deal with writing to the logs.


    When manipulating logs, you'll include <boost/log/log_ext.hpp>
*/

#include <boost/log/functions.hpp>
#include <boost/limits.hpp>
#include <iostream>
#include <fstream>
#include <time.h>

// *** Appenders
void my_write_to_cout(const std::string &, const std::string &msg) { 
    std::cout << msg; 
}


// *** Modifiers
void prefix_time(const std::string &, std::string & msg) {
    char time_buff[ 20]; time_t t = time(0); tm details = *localtime( &t);
    sprintf( time_buff, "%02d:%02d:%02d  ", details.tm_hour, details.tm_min, details.tm_sec);
    msg = time_buff + msg;
}

void init_logs() {
    using namespace boost::logging;

    // for all logs
    manipulate_logs("*")
        // [type_of_message] original_message append_enter_if_needed
        .add_modifier( &prefix_time, INT_MAX )
        .add_modifier( &append_enter)
        // all messages are written to cout
        .add_appender( my_write_to_cout);

    // for "app" and descenants
    manipulate_logs("app*")
        // <time> [type_of_message] original_message append_enter_if_needed
        .add_modifier( &prepend_prefix)
        // "app*" messages are written to file as well
        .add_appender( write_to_file("out.txt") );

    // for "app" only
    manipulate_logs("app")
        // <time> [Thread ID] [type_of_message] original_message append_enter_if_needed
        .add_modifier( &prepend_thread_id, 0);

    // 'app' only and dbg messages are written to Output Debug Window as well
#ifdef BOOST_LOG_WIN32
    manipulate_logs("app").add_appender( write_to_dbg_wnd);
    manipulate_logs("*dbg*").add_appender( write_to_dbg_wnd);
#endif

    flush_log_cache();
}



