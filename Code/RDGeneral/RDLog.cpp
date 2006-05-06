// $Id$
//
// Copyright (C)  2005-2006 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "RDLog.h"
#include <boost/log/functions.hpp>
#if defined(BOOST_HAS_THREADS2)
#include <boost/log/extra/functions_ts.hpp>
#endif
#include <iostream>
namespace logging = boost::logging;

BOOST_DEFINE_LOG(rdAppLog,"rdApp")
BOOST_DEFINE_LOG(rdDebugLog,"rdApp.debug")
BOOST_DEFINE_LOG(rdInfoLog,"rdApp.info")
BOOST_DEFINE_LOG(rdErrorLog,"rdApp.error")
BOOST_DEFINE_LOG(rdWarningLog,"rdApp.warning")
BOOST_DEFINE_LOG(rdStatusLog,"rdApp.status")

namespace RDLog {
  void write_to_cout(const std::string &, const std::string &msg) {
    std::cout << msg; std::cout.flush();
  }
  void write_to_cerr(const std::string &, const std::string &msg) {
    std::cerr << msg; std::cerr.flush();
  }

  void InitLogs(){
    static bool callOnce=true;
    if(!callOnce) return;
    callOnce=false;

    // turn off caching:
    logging::set_log_caching(false);
    logging::manipulate_logs("rdApp.*")
      .add_modifier(logging::prepend_time("[$hh:$mm:$ss] "),
		    logging::DEFAULT_INDEX-10);
    logging::manipulate_logs("rdApp.info")
      .add_appender(write_to_cout,
		    logging::DEFAULT_INDEX+1);
#if defined(BOOST_HAS_THREADS2)
    logging::manipulate_logs("rdApp.error")
      .add_appender(logging::ts_appender(write_to_cerr,100),
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.warning")
      .add_appender(logging::ts_appender(write_to_cerr,100),
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.status")
      .add_appender(logging::ts_appender(write_to_cerr,100),
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.debug")
      .add_appender(logging::ts_appender(write_to_cerr,100),
		    logging::DEFAULT_INDEX+1);
#else
    logging::manipulate_logs("rdApp.error")
      .add_appender(write_to_cerr,
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.warning")
      .add_appender(write_to_cerr,
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.status")
      .add_appender(write_to_cerr,
		    logging::DEFAULT_INDEX+1);
    logging::manipulate_logs("rdApp.debug")
      .add_appender(write_to_cerr,
		    logging::DEFAULT_INDEX+1);
#endif
    // start with the debug log disabled:
    logging::disable_logs("rdApp.debug");
  };
}
