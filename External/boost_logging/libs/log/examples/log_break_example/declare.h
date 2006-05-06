// declare.h - declare the logs
#ifndef SAMPLE_DECLARE_H
#define SAMPLE_DECLARE_H

#include <boost/log/log.hpp>

BOOST_DECLARE_LOG(app)
BOOST_DECLARE_LOG(dbg)
BOOST_DECLARE_LOG(warn)

void init_logs();

#endif
