//
// Copyright (C)  2005-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef _RDLOG_H_29JUNE2005_
#define _RDLOG_H_29JUNE2005_

#if 1
#include "BoostStartInclude.h"
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include "BoostEndInclude.h"
#include <iostream>
namespace boost {
namespace logging {

typedef boost::iostreams::tee_device<std::ostream, std::ostream> RDTee;
typedef boost::iostreams::stream<RDTee> RDTeeStream;

class RDKIT_RDGENERAL_EXPORT rdLogger {
 public:
  std::ostream *dp_dest;
  bool df_owner, df_enabled;

  RDTee *tee;
  RDTeeStream *teestream;

  rdLogger(std::ostream *dest, bool owner = false)
      : dp_dest(dest),
        df_owner(owner),
        df_enabled(true),
        tee(NULL),
        teestream(NULL){};

  //! Sets a stream to tee the output to.
  void SetTee(std::ostream &stream) {
    if (dp_dest) {
      delete teestream;
      delete tee;
      tee = new RDTee(*dp_dest, stream);
      teestream = new RDTeeStream(*tee);
    }
  }
  //! Remove our tee if it's set.
  void ClearTee() {
    if (dp_dest) {
      delete teestream;
      delete tee;
      tee = NULL;
      teestream = NULL;
    }
  }
  ~rdLogger() {
    if (dp_dest) {
      dp_dest->flush();
      if (df_owner) {
        delete dp_dest;
      }
      dp_dest = NULL;
    }
    delete teestream;
    teestream = NULL;
    delete tee;
    tee = NULL;
  }

 private:
  // disable copy ctor and assignment
  rdLogger(const rdLogger &);
  rdLogger &operator=(const rdLogger &);
};
RDKIT_RDGENERAL_EXPORT void enable_logs(const char *arg);
RDKIT_RDGENERAL_EXPORT void enable_logs(const std::string &arg);
RDKIT_RDGENERAL_EXPORT void disable_logs(const char *arg);
RDKIT_RDGENERAL_EXPORT void disable_logs(const std::string &arg);
}
}
namespace RDLog {
RDKIT_RDGENERAL_EXPORT std::ostream &toStream(std::ostream &);
}
#define BOOST_LOG(__arg__)                                      \
  if ((__arg__) && (__arg__->dp_dest) && (__arg__->df_enabled)) \
  RDLog::toStream((__arg__->teestream) ? *(__arg__->teestream)  \
                                       : *(__arg__->dp_dest))

RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdAppLog;
RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdDebugLog;
RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdInfoLog;
RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdErrorLog;
RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdWarningLog;
RDKIT_RDGENERAL_EXPORT extern std::shared_ptr<boost::logging::rdLogger> rdStatusLog;

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
RDKIT_RDGENERAL_EXPORT void InitLogs();
}
#endif
