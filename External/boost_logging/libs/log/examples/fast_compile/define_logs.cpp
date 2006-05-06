/*
    Sample to show that using Boost Logging library incurs in almost no
    overhead in compilation-time.

    This file compiles VERY FAST.
*/

#include "declare_logs.h"

BOOST_DEFINE_LOG(app, "app")
BOOST_DEFINE_LOG(dbg, "app.dbg")
BOOST_DEFINE_LOG(err, "app.err")
BOOST_DEFINE_LOG(warn, "app.warn")
BOOST_DEFINE_LOG(info, "info")

BOOST_DEFINE_LOG(charts::dbg, "charts.dbg")
BOOST_DEFINE_LOG(charts::gui, "charts.gui")
