#ifndef SAMPLE_LOGS_H
#define SAMPLE_LOGS_H

#include <boost/log/log.hpp>

BOOST_DECLARE_LOG(app)
BOOST_DECLARE_LOG(dbg)
BOOST_DECLARE_LOG(err)

void init_logs();


#endif
