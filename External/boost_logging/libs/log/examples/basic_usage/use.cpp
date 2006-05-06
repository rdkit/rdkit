// use.cpp - write to the logs
#include "declare.h"

int main() {
    init_logs();

    int i = 1, j = 2, k = 3;
    BOOST_LOG(app) << "testing " << i << '-' << j << '-' << k  << std::endl;
    BOOST_LOG(dbg) << "this is a debug message, i=" << i << std::endl;

    BOOST_LOG(info) << "I just wanted to tell you something...."  << std::endl;
    BOOST_LOG(info) << "Logged-on Users approach max. limit"  << std::endl;
    BOOST_LOG(dbg) << "Too many users!"  << std::endl;

    BOOST_LOG(app) << "an application message" << std::endl;
    BOOST_LOG(info) << "a simple info" << std::endl;
}
