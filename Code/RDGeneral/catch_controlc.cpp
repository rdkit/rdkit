//
//  Copyright (C) 2026 Steven Kearnes
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <atomic>
#include <csignal>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <vector>
#endif
#include <catch2/catch_all.hpp>
#include "ControlCHandler.h"

using namespace RDKit;

namespace {
std::atomic<int> outsideCalls{0};
void outsideHandler(int) { ++outsideCalls; }

// The SIGINT handler installed now.
using SignalHandler = void (*)(int);
SignalHandler currentHandler() {
  auto handler = std::signal(SIGINT, SIG_DFL);
  std::signal(SIGINT, handler);
  return handler;
}
}  // namespace

TEST_CASE("ControlCHandler restores the handler it replaced") {
  auto original = std::signal(SIGINT, outsideHandler);
  SECTION("one handler") {
    {
      ControlCHandler hdlr;
      CHECK(currentHandler() != &outsideHandler);
    }
    CHECK(currentHandler() == &outsideHandler);
  }
  SECTION("nested handlers") {
    {
      ControlCHandler outer;
      { ControlCHandler inner; }
      CHECK(currentHandler() != &outsideHandler);
    }
    CHECK(currentHandler() == &outsideHandler);
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("handlers on several threads") {
    // As in substructure searches run in parallel, each of which has one.
    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < 8; ++i) {
      threads.emplace_back([] {
        for (unsigned int j = 0; j < 10000; ++j) {
          ControlCHandler hdlr;
        }
      });
    }
    for (auto &thread : threads) {
      thread.join();
    }
    CHECK(currentHandler() == &outsideHandler);
    CHECK(!ControlCHandler::getGotSignal());
  }
#endif
  std::signal(SIGINT, original);
}

TEST_CASE("ControlCHandler catches SIGINT") {
  auto original = std::signal(SIGINT, outsideHandler);
  outsideCalls = 0;
  {
    ControlCHandler outer;
    {
      ControlCHandler inner;
      std::raise(SIGINT);
      CHECK(ControlCHandler::getGotSignal());
    }
    // The flag stays set while any handler is alive.
    CHECK(ControlCHandler::getGotSignal());
    CHECK(outsideCalls == 0);
    // Catching the signal restored the outside handler.
    std::raise(SIGINT);
    CHECK(outsideCalls == 1);
    ControlCHandler::reset();
    CHECK(!ControlCHandler::getGotSignal());
    std::raise(SIGINT);
    CHECK(ControlCHandler::getGotSignal());
    CHECK(outsideCalls == 1);
  }
  CHECK(!ControlCHandler::getGotSignal());
  CHECK(currentHandler() == &outsideHandler);
  ControlCHandler::reset();
  CHECK(currentHandler() == &outsideHandler);
  std::signal(SIGINT, original);
}
