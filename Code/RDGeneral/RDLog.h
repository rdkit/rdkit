//
// Copyright (C)  2005-2022 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDLOG_H_29JUNE2005
#define RDLOG_H_29JUNE2005

#if 1
#include "BoostStartInclude.h"
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include "BoostEndInclude.h"
#include <iostream>
#include <vector>

#include <vector>
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
        tee(nullptr),
        teestream(nullptr) {}

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
      tee = nullptr;
      teestream = nullptr;
    }
  }
  ~rdLogger() {
    if (dp_dest) {
      dp_dest->flush();
      if (df_owner) {
        delete dp_dest;
      }
      dp_dest = nullptr;
    }
    delete teestream;
    teestream = nullptr;
    delete tee;
    tee = nullptr;
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
RDKIT_RDGENERAL_EXPORT std::string log_status();
}  // namespace logging
}  // namespace boost
namespace RDLog {
RDKIT_RDGENERAL_EXPORT std::ostream &toStream(std::ostream &);
}
#define BOOST_LOG(__arg__)                                      \
  if ((__arg__) && (__arg__->dp_dest) && (__arg__->df_enabled)) \
  RDLog::toStream((__arg__->teestream) ? *(__arg__->teestream)  \
                                       : *(__arg__->dp_dest))

using RDLogger = std::shared_ptr<boost::logging::rdLogger>;

RDKIT_RDGENERAL_EXPORT extern RDLogger rdAppLog;
RDKIT_RDGENERAL_EXPORT extern RDLogger rdDebugLog;
RDKIT_RDGENERAL_EXPORT extern RDLogger rdInfoLog;
RDKIT_RDGENERAL_EXPORT extern RDLogger rdErrorLog;
RDKIT_RDGENERAL_EXPORT extern RDLogger rdWarningLog;
RDKIT_RDGENERAL_EXPORT extern RDLogger rdStatusLog;

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

using RDLoggerList = std::vector<RDLogger>;
class RDKIT_RDGENERAL_EXPORT LogStateSetter : public boost::noncopyable {
 public:
  //! enables only the logs in the list, the current state will be restored when
  //! this object is destroyed
  LogStateSetter(RDLoggerList toEnable);
  //! disables all logs, the current state will be restored when this object is
  //! destroyed
  LogStateSetter();
  ~LogStateSetter();

 private:
  std::uint64_t d_origState = 0;
};

}  // namespace RDLog
#endif
