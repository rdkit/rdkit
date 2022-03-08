//
//  Copyright (C) 2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <sstream>
#include "catch.hpp"
#include "RDLog.h"

TEST_CASE("LogStateSetter") {
  RDLog::RDLoggerList allLogs({rdErrorLog, rdWarningLog, rdInfoLog});
  SECTION("disable all") {
    for (auto strm : allLogs) {
      // explicitly enable the stream so that we know something is happening
      std::stringstream ostrm;
      strm->df_enabled = true;
      strm->SetTee(ostrm);
      {
        RDLog::LogStateSetter disabler;
        BOOST_LOG(strm) << "should be silent" << std::endl;
        auto txt = ostrm.str();
        CHECK(txt.find("should") == std::string::npos);
      }
      BOOST_LOG(strm) << "should not be silent" << std::endl;
      auto txt = ostrm.str();
      CHECK(txt.find("should") != std::string::npos);
      strm->ClearTee();
    }
  }
  SECTION("enable one") {
    for (auto strm : allLogs) {
      strm->df_enabled = true;
      RDLog::LogStateSetter disabler;
      {
        RDLog::LogStateSetter justone(RDLog::RDLoggerList({strm}));
        for (auto strm2 : allLogs) {
          std::stringstream ostrm;
          strm2->SetTee(ostrm);
          BOOST_LOG(strm2) << "should not be silent" << std::endl;
          auto txt = ostrm.str();
          if (strm == strm2) {
            CHECK(txt.find("should") != std::string::npos);
          } else {
            CHECK(txt.find("should") == std::string::npos);
          }

          strm2->ClearTee();
        }
      }

      // should be disabled again
      std::stringstream ostrm;
      strm->SetTee(ostrm);
      BOOST_LOG(strm) << "should be silent" << std::endl;
      auto txt = ostrm.str();
      CHECK(txt.find("should") == std::string::npos);
      strm->ClearTee();
    }
  }
}