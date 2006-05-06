//
// Copyright (C)  2005 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//

#ifndef _RDLOG_H_29JUNE2005_
#define _RDLOG_H_29JUNE2005_

#define BOOST_LOG_NO_LIB
#include <boost/log/log.hpp>
BOOST_DECLARE_LOG(rdAppLog)
BOOST_DECLARE_LOG(rdDebugLog)
BOOST_DECLARE_LOG(rdInfoLog)
BOOST_DECLARE_LOG(rdErrorLog)
BOOST_DECLARE_LOG(rdWarningLog)
BOOST_DECLARE_LOG(rdStatusLog)

namespace RDLog {
  void InitLogs();
}
#endif

