//
// Copyright (C)  2005-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//

#ifndef _RDLOG_H_29JUNE2005_
#define _RDLOG_H_29JUNE2005_

#if 1
#include <iostream>
extern std::ostream *rdAppLog;
extern std::ostream *rdDebugLog;
extern std::ostream *rdInfoLog;
extern std::ostream *rdErrorLog;
extern std::ostream *rdWarningLog;
extern std::ostream *rdStatusLog;
namespace boost {
  namespace logging {
    void enable_logs(const char *arg);
    void enable_logs(const std::string &arg);
    void disable_logs(const char *arg);
    void disable_logs(const std::string &arg);
  }
}
namespace RDLog {
  std::ostream &toStream(std::ostream &);
}
#define BOOST_LOG(__arg__) if(!__arg__) ; else RDLog::toStream(*__arg__)

#else
#define BOOST_LOG_NO_LIB
#include <boost/log/log.hpp>
BOOST_DECLARE_LOG(rdAppLog)
BOOST_DECLARE_LOG(rdDebugLog)
BOOST_DECLARE_LOG(rdInfoLog)
BOOST_DECLARE_LOG(rdErrorLog)
BOOST_DECLARE_LOG(rdWarningLog)
BOOST_DECLARE_LOG(rdStatusLog)
#endif
namespace RDLog {
  void InitLogs();
}
#endif

