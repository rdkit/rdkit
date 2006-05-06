#ifndef DECLARE_LOGS_H
#define DECLARE_LOGS_H

#include <boost/log/log.hpp>

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

#endif
