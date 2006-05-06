// use.cpp - write to the logs
#include "declare.h"
#include <boost/log/log_impl.hpp>

using namespace boost::logging;

int main() {
    init_logs();

    int i = 1, j = 2, k = 3;

    enable_logs("*", level::dbg);
    BOOST_LOGL(lg,dbg) << "testing " << i << '-' << j << '-' << k  << std::endl;
    BOOST_LOGL(lg,dbg) << "this is a debug message, i=" << i << std::endl;

    enable_logs("*", level::warn);
    BOOST_LOGL(lg,info) << "This won't be written...."  << std::endl;
    BOOST_LOGL(lg,warn) << "This is a warning (written)"  << std::endl;
    BOOST_LOGL(lg,dbg) << "Too many users!"  << std::endl;

    enable_logs("*", level::dbg);
    BOOST_LOGL(lg,info) << "This won't be written...."  << std::endl;
    BOOST_LOGL(lg,warn) << "This won't be written"  << std::endl;
    BOOST_LOGL(lg,dbg) << "Too many users!"  << std::endl;

    enable_logs("*", level::err);
    BOOST_LOGL(lg,info) << "This won't be written...."  << std::endl;
    BOOST_LOGL(lg,warn) << "This won't be written"  << std::endl;
    BOOST_LOGL(lg,dbg) << "This won't be written"  << std::endl;
    BOOST_LOGL(lg,err) << "an error" << std::endl;

    enable_logs("*", level::fatal);
    BOOST_LOGL(lg,info) << "This won't be written...."  << std::endl;
    BOOST_LOGL(lg,warn) << "This won't be written"  << std::endl;
    BOOST_LOGL(lg,dbg) << "This won't be written"  << std::endl;
    BOOST_LOGL(lg,err) << "This won't be written" << std::endl;
    BOOST_LOGL(lg,fatal) << "a fatal error" << std::endl;
}
