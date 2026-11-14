//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "ControlCHandler.h"

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

namespace RDKit {

namespace {
using SignalHandler = void (*)(int);

// The signal handler reads and writes these, so they must be lock-free.
static_assert(std::atomic<bool>::is_always_lock_free);
static_assert(std::atomic<SignalHandler>::is_always_lock_free);

// The state is defined here, once, rather than in the header, so that every
// library shares it.  The mutex orders the handlers' constructors and
// destructors, which may run on several threads at once.
#ifdef RDK_BUILD_THREADSAFE_SSS
std::mutex handlerMutex;
#endif
unsigned int numHandlers = 0;
std::atomic<SignalHandler> prevHandler{SIG_DFL};
}  // namespace

std::atomic<bool> ControlCHandler::d_gotSignal{false};

ControlCHandler::ControlCHandler() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(handlerMutex);
#endif
  if (numHandlers++ == 0) {
    prevHandler = std::signal(SIGINT, signalHandler);
  }
}

ControlCHandler::~ControlCHandler() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(handlerMutex);
#endif
  if (--numHandlers == 0) {
    std::signal(SIGINT, prevHandler);
    d_gotSignal = false;
  }
}

void ControlCHandler::signalHandler(int signalNumber) {
  if (signalNumber == SIGINT) {
    d_gotSignal = true;
    std::signal(SIGINT, prevHandler);
  }
}

void ControlCHandler::reset() {
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::lock_guard<std::mutex> lock(handlerMutex);
#endif
  d_gotSignal = false;
  if (numHandlers) {
    std::signal(SIGINT, signalHandler);
  }
}

}  // namespace RDKit
