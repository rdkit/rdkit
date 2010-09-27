// $Id$
//
// Copyright (C)  2005-2010 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDLog.h"

#if 1
#include <iomanip>
#include <string>
#include <time.h>

namespace {
  // this is a "bit" of a hack to work around shared/static library problems
  // on windows
  boost::logging::rdLogger cerrLogger(&std::cerr);  
  boost::logging::rdLogger coutLogger(&std::cout);  
}

boost::logging::rdLogger *rdAppLog=0;
boost::logging::rdLogger *rdDebugLog=0;
boost::logging::rdLogger *rdInfoLog=&coutLogger;
boost::logging::rdLogger *rdErrorLog=&cerrLogger;
boost::logging::rdLogger *rdWarningLog=&cerrLogger;
boost::logging::rdLogger *rdStatusLog=0;

namespace boost {
  namespace logging {

    void enable_logs(const char *arg) { enable_logs(std::string(arg));};
    void enable_logs(const std::string &arg) {
      // Yes... this is extremely crude
      if(arg=="rdApp.debug"||arg=="rdApp.*"){
        if(rdDebugLog) rdDebugLog->df_enabled=true;
      }
      if(arg=="rdApp.info"||arg=="rdApp.*"){
        if(rdInfoLog) rdInfoLog->df_enabled=true;
      }
      if(arg=="rdApp.warning"||arg=="rdApp.*"){
        if(rdWarningLog) rdWarningLog->df_enabled=true;
      }
      if(arg=="rdApp.error"||arg=="rdApp.*"){
        if(rdErrorLog) rdErrorLog->df_enabled=true;
      }
    };
    void disable_logs(const char *arg) {disable_logs(std::string(arg));};
    void disable_logs(const std::string &arg) {
      // Yes... this is extremely crude
      if(arg=="rdApp.debug"||arg=="rdApp.*"){
        if(rdDebugLog) rdDebugLog->df_enabled=false;
      }
      if(arg=="rdApp.info"||arg=="rdApp.*"){
        if(rdInfoLog) rdInfoLog->df_enabled=false;
      }
      if(arg=="rdApp.warning"||arg=="rdApp.*"){
        if(rdWarningLog) rdWarningLog->df_enabled=false;
      }
      if(arg=="rdApp.error"||arg=="rdApp.*"){
        if(rdErrorLog) rdErrorLog->df_enabled=false;
      }
    };
  }
}


namespace RDLog {
  void InitLogs(){
    rdDebugLog=new boost::logging::rdLogger(&std::cerr);
    rdInfoLog=new boost::logging::rdLogger(&std::cout);
    rdWarningLog=new boost::logging::rdLogger(&std::cerr);
    rdErrorLog=new boost::logging::rdLogger(&std::cerr);
  }
  std::ostream &toStream(std::ostream &logstrm) {
    time_t t = time(0); 
    tm details = *localtime( &t);
    logstrm << "["<<std::setw(2)<<std::setfill('0')<<details.tm_hour<<":"<<std::setw(2)<<std::setfill('0')<<details.tm_min<<":"<<std::setw(2)<<std::setfill('0')<<int(details.tm_sec)<<"] ";
    return logstrm;
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
