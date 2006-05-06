// Sample for Boost Logging
//

/*
    Sample to show that using Boost Logging library incurs in almost no
    overhead in compilation-time.

    This file compiles VERY FAST.
*/


#include <boost/log/log.hpp>
#include "declare_logs.h"

void init_logs(); // defined in manipulate_logs.cpp

int main() {
    init_logs();

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

    BOOST_LOG(app) << "some app message" << std::endl;
    BOOST_LOG(dbg) << "You are in debug mode" << std::endl;
    BOOST_LOG(err) << "some error" << std::endl;

    BOOST_LOG(dbg) << "this will be written" << std::endl;
    BOOST_LOG(info) << "a simple info" << std::endl;
}

