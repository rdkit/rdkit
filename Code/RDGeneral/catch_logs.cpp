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
#include <fstream>
#include <sstream>
#include <catch2/catch_all.hpp>
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

TEST_CASE("GitHub Issue #5172", "[bug][logging]") {
  std::stringstream err_ostrm;
  std::stringstream warn_ostrm;
  rdErrorLog->df_enabled = true;
  rdWarningLog->df_enabled = true;
  rdErrorLog->SetTee(err_ostrm);
  rdWarningLog->SetTee(warn_ostrm);

  {
    RDLog::LogStateSetter disabler;

    BOOST_LOG(rdErrorLog) << "should be silent" << std::endl;
    auto txt = err_ostrm.str();
    CHECK(txt.find("should") == std::string::npos);

    BOOST_LOG(rdWarningLog) << "should be silent" << std::endl;
    txt = warn_ostrm.str();
    CHECK(txt.find("should") == std::string::npos);

    {
      // The second setter overrides the first one:
      RDLog::LogStateSetter disabler2({rdWarningLog});

      BOOST_LOG(rdErrorLog) << "should be silent" << std::endl;
      txt = err_ostrm.str();
      CHECK(txt.find("should") == std::string::npos);

      BOOST_LOG(rdWarningLog) << "should not be silent" << std::endl;
      txt = warn_ostrm.str();
      CHECK(txt.find("should") != std::string::npos);
      warn_ostrm.clear();
    }
  }

  // Both setters are destroyed, and we revert to initial state

  BOOST_LOG(rdErrorLog) << "should not be silent" << std::endl;
  auto txt = err_ostrm.str();
  CHECK(txt.find("should") != std::string::npos);

  BOOST_LOG(rdWarningLog) << "should not be silent" << std::endl;
  txt = warn_ostrm.str();
  CHECK(txt.find("should") != std::string::npos);

  rdErrorLog->ClearTee();
  rdWarningLog->ClearTee();
}

TEST_CASE("CaptureErrorLog") {
  RDLog::InitLogs();

  SECTION("basic capture") {
    RDLog::CaptureErrorLog capture;
    BOOST_LOG(rdErrorLog) << "test error" << std::endl;
    CHECK(capture.messages().find("test error") != std::string::npos);
  }

  SECTION("empty when nothing logged") {
    RDLog::CaptureErrorLog capture;
    CHECK(capture.messages().empty());
  }

  SECTION("enables a disabled log and restores its state") {
    rdErrorLog->df_enabled = false;
    {
      RDLog::CaptureErrorLog capture;
      CHECK(rdErrorLog->df_enabled);
      BOOST_LOG(rdErrorLog) << "captured despite disabled" << std::endl;
      CHECK(capture.messages().find("captured despite disabled") !=
            std::string::npos);
    }
    // enabled state restored to false
    CHECK(!rdErrorLog->df_enabled);
    rdErrorLog->df_enabled = true;
  }

  SECTION("restores original stream destination") {
    {
      RDLog::CaptureErrorLog capture;
      BOOST_LOG(rdErrorLog) << "captured" << std::endl;
    }
    // After capture, logs go to the original destination again
    std::stringstream ss;
    rdErrorLog->SetTee(ss);
    BOOST_LOG(rdErrorLog) << "after restore" << std::endl;
    rdErrorLog->ClearTee();
    CHECK(ss.str().find("after restore") != std::string::npos);
  }

  SECTION("does not capture other logs") {
    RDLog::CaptureErrorLog capture;
    BOOST_LOG(rdWarningLog) << "test warning" << std::endl;
    BOOST_LOG(rdErrorLog) << "test error" << std::endl;
    CHECK(capture.messages().find("test error") != std::string::npos);
    CHECK(capture.messages().find("test warning") == std::string::npos);
  }

  SECTION("restores dp_dest to original stream") {
    auto *original_dest = rdErrorLog->dp_dest;
    {
      RDLog::CaptureErrorLog capture;
      CHECK(rdErrorLog->dp_dest != original_dest);
    }
    CHECK(rdErrorLog->dp_dest == original_dest);
  }

  SECTION("re-enables log even when silenced by LogStateSetter") {
    RDLog::CaptureErrorLog outer;
    RDLog::LogStateSetter silence;

    BOOST_LOG(rdErrorLog) << "silenced" << std::endl;
    CHECK(outer.messages().empty());

    {
      RDLog::CaptureErrorLog inner;
      BOOST_LOG(rdErrorLog) << "captured despite silence" << std::endl;
      CHECK(inner.messages().find("captured despite silence") !=
            std::string::npos);
    }

    BOOST_LOG(rdErrorLog) << "silenced again" << std::endl;
    CHECK(outer.messages().empty());
  }

  SECTION("nested captures") {
    RDLog::CaptureErrorLog outer;
    BOOST_LOG(rdErrorLog) << "outer message" << std::endl;
    {
      RDLog::CaptureErrorLog inner;
      BOOST_LOG(rdErrorLog) << "inner message" << std::endl;
      CHECK(inner.messages().find("inner message") != std::string::npos);
      CHECK(inner.messages().find("outer message") == std::string::npos);
    }
    CHECK(outer.messages().find("outer message") != std::string::npos);
    CHECK(outer.messages().find("inner message") == std::string::npos);
  }
}

TEST_CASE("Tee to file") {
  const std::string filename = "error_log.txt";
  rdErrorLog->SetTee(filename);
  BOOST_LOG(rdErrorLog) << "should not be silent" << std::endl;
  std::ifstream istrm(filename);
  std::string txt;
  CHECK(std::getline(istrm, txt));
  CHECK(txt.find("should") != std::string::npos);
}
