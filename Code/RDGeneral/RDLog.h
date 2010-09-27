//
// Copyright (C)  2005-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RDLOG_H_29JUNE2005_
#define _RDLOG_H_29JUNE2005_

#if 1
#include <iostream>
namespace boost {
  namespace logging {
    class rdLogger{
    public:
    rdLogger(std::ostream *dest,bool owner=false) : dp_dest(dest), df_owner(owner),
        df_enabled(true) {};
      std::ostream *dp_dest;
      bool df_owner,df_enabled;
      ~rdLogger(){
        if(dp_dest){
          dp_dest->flush();
          if(df_owner){
            delete dp_dest;
          }
        }
      }
    };
    void enable_logs(const char *arg);
    void enable_logs(const std::string &arg);
    void disable_logs(const char *arg);
    void disable_logs(const std::string &arg);
  }
}
namespace RDLog {
  std::ostream &toStream(std::ostream &);
}
#define BOOST_LOG(__arg__) if((!__arg__)||(!__arg__->dp_dest)||!(__arg__->df_enabled)) ; else RDLog::toStream(*(__arg__->dp_dest))


extern boost::logging::rdLogger *rdAppLog;
extern boost::logging::rdLogger *rdDebugLog;
extern boost::logging::rdLogger *rdInfoLog;
extern boost::logging::rdLogger *rdErrorLog;
extern boost::logging::rdLogger *rdWarningLog;
extern boost::logging::rdLogger *rdStatusLog;


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

