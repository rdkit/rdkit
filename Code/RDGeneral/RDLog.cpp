// $Id$
//
// Copyright (C)  2005-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "RDLog.h"

#if 1
#include <iomanip>
#include <time.h>

std::ostream *rdAppLog=0;
std::ostream *rdDebugLog=0;
std::ostream *rdInfoLog=0;
std::ostream *rdErrorLog=0;
std::ostream *rdWarningLog=0;
std::ostream *rdStatusLog=0;
namespace boost {
  namespace logging {
    void enable_logs(const char *arg) {};
    void enable_logs(const std::string &arg) {};
    void disable_logs(const char *arg) {};
    void disable_logs(const std::string &arg) {};
  }
}


namespace RDLog {
  void InitLogs(){
    rdAppLog=&std::cout;
    rdDebugLog=&std::cerr;
    rdInfoLog=&std::cout;
    rdErrorLog=&std::cerr;
    rdWarningLog=&std::cerr;
    rdStatusLog=&std::cout;
  }
  std::ostream &toStream(std::ostream &strm) {
    time_t t = time(0); 
    tm details = *localtime( &t);
    strm << "["<<std::setw(2)<<std::setfill('0')<<details.tm_hour<<":"<<std::setw(2)<<std::setfill('0')<<details.tm_min<<":"<<std::setw(2)<<std::setfill('0')<<int(details.tm_sec)<<"] ";
    return strm;
  }

}

#else
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
#endif
