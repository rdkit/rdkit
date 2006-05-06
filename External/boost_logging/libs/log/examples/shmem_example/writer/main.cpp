// main.cpp 

#include <boost/log/log.hpp>
#include <boost/log/functions.hpp>
#include <boost/log/extra/shared_memory.hpp>

BOOST_DECLARE_LOG(lg)
BOOST_DEFINE_LOG(lg, "log")


int main() {
    using namespace boost::logging;

    manipulate_logs("*")
        .add_modifier(prepend_time("$hh:$mm:$ss "), "time", DEFAULT_INDEX + 1 )
        .add_modifier(&prepend_prefix, "prefix")
        .add_appender( shared_memory_appender("/shared_test") );

    flush_log_cache();

    int i = 1, j = 2, k = 3;
    BOOST_LOG(lg) << "testing " << i << '-' << j << '-' << k  << std::endl;
    BOOST_LOG(lg) << "this is a debug message, i=" << i << std::endl;

    BOOST_LOG(lg) << "I just wanted to tell you something...."  << std::endl;
    BOOST_LOG(lg) << "Logged-on Users approach max. limit"  << std::endl;
    BOOST_LOG(lg) << "Too many users!"  << std::endl;

    BOOST_LOG(lg) << "an application message" << std::endl;
    BOOST_LOG(lg) << "a simple info" << std::endl;
}
