// main.cpp 

#include <boost/log/log.hpp>
#include <boost/log/functions.hpp>

BOOST_DECLARE_LOG(lg)
BOOST_DEFINE_LOG(lg, "log")


int main() {
    using namespace boost::logging;

    manipulate_logs("*")
        .add_modifier(prepend_time("$hh:$mm:$ss "), "time", DEFAULT_INDEX + 1 )
        .add_modifier(&prepend_prefix, "prefix")

        .add_appender(&write_to_cout, "cout")        
        .add_appender(write_to_file("out.txt"), "outfile" )
    #ifdef BOOST_LOG_WIN32
        .add_appender(write_to_dbg_wnd, "dbg" );
    #endif
        ;
    flush_log_cache();

    BOOST_LOG(lg) << "writing to: cout, debug, out.txt file" << std::endl;
    manipulate_logs("*").del_appender("dbg");

    BOOST_LOG(lg) << "writing to: cout, out.txt file"  << std::endl;

    manipulate_logs("*").del_appender("outfile");
    BOOST_LOG(lg) << "writing to: cout only"  << std::endl;

    manipulate_logs("*").add_appender(write_to_file("out.txt",false), "outfile" );
    BOOST_LOG(lg) << "writing to: cout, out.txt file"  << std::endl;

    manipulate_logs("*").del_modifier("prefix");
    BOOST_LOG(lg) << "writing to: cout, out.txt file - without [log] prefix"  << std::endl;

    manipulate_logs("*").del_modifier("time");
    BOOST_LOG(lg) << "writing to: cout, out.txt file - without [log] prefix and without time prefix"  << std::endl;

    manipulate_logs("*").add_modifier(prepend_time("$hh:$mm:$ss "), "time", DEFAULT_INDEX + 1 );
    BOOST_LOG(lg) << "writing to: cout, out.txt file - without [log] prefix but with time prefix"  << std::endl;

}
