//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef CONTROLCHANDLER_H
#define CONTROLCHANDLER_H

#include <atomic>
#include <csignal>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif
#include <stdexcept>

namespace RDKit {

class ControlCCaught : public std::runtime_error {
 public:
  explicit ControlCCaught()
      : std::runtime_error("The process was interrupted with Ctrl+c") {};
};

//! This class catches a control-C/SIGINT and sets the flag d_gotSignal
//! if one is received.  It is intended to be used inside a long
//! C++ calculation called from Python which intercepts the signal
//! handler.  The C++ code must check the value of d_gotSignal
//! periodically and act accordingly.  The destructor resets
//! the signal handler and flag after the last concurrent user exits.
//! Example usage, inside a boost::python wrapper:
//!  ResultsObject results;
//!  {
//!   NOGIL gil;
//!    results = someFunction();
//!  }
//!  if (results.getCancelled()) {
//!    throw_runtime_error("someFunction cancelled");
//!  }
//! It's important that the exception is thrown once the GIL has been
//! released, otherwise a crash is inevitable at some future point.
class ControlCHandler {
 public:
  ControlCHandler() {
#ifdef RDK_BUILD_THREADSAFE_SSS
    std::lock_guard<std::mutex> lock(d_mutex);
#endif
    HandlerTransitionGuard transition;
    if (d_handlerCount++ == 0) {
      installHandler();
    }
  }
  ControlCHandler(const ControlCHandler &) = delete;
  ControlCHandler(ControlCHandler &&) = delete;
  ControlCHandler &operator=(const ControlCHandler &) = delete;
  ControlCHandler &operator=(ControlCHandler &&) = delete;
  ~ControlCHandler() {
#ifdef RDK_BUILD_THREADSAFE_SSS
    std::lock_guard<std::mutex> lock(d_mutex);
#endif
    HandlerTransitionGuard transition;
    if (--d_handlerCount == 0) {
      std::signal(SIGINT, d_prev_handler.load(std::memory_order_relaxed));
      d_gotSignal.clear(std::memory_order_relaxed);
    }
  }
  static bool getGotSignal() {
    return d_gotSignal.test(std::memory_order_relaxed);
  }
  static void signalHandler(int signalNumber) {
    if (signalNumber == SIGINT) {
      d_gotSignal.test_and_set();
      restoreInterruptedHandler();
    }
  }
  static void reset() {
#ifdef RDK_BUILD_THREADSAFE_SSS
    std::lock_guard<std::mutex> lock(d_mutex);
#endif
    HandlerTransitionGuard transition;
    installHandler();
  }

 private:
  using SignalHandler = void (*)(int);
  class HandlerTransitionGuard {
   public:
    HandlerTransitionGuard() {
      while (d_handlerTransition.test_and_set(std::memory_order_acquire)) {
      }
    }
    ~HandlerTransitionGuard() {
      d_handlerTransition.clear(std::memory_order_release);
      // A signal arriving during a transition records cancellation and leaves
      // restoration to the code completing that transition.
      restoreInterruptedHandler();
    }
  };

  static void restoreInterruptedHandler() {
    // Never wait in a signal handler: it may have interrupted a transition
    // on this thread. Ordinary transitions and signal-driven restoration use
    // the same gate so an old handler cannot overwrite a new registration.
    if (d_gotSignal.test() &&
        !d_handlerTransition.test_and_set(std::memory_order_acquire)) {
      if (d_gotSignal.test()) {
        std::signal(SIGINT, d_prev_handler.load(std::memory_order_relaxed));
      }
      d_handlerTransition.clear(std::memory_order_release);
    }
  }

  static void installHandler() {
    d_gotSignal.clear();
    const auto previous = std::signal(SIGINT, signalHandler);
    // reset() may have armed us before the first RAII scope, or may rearm an
    // active scope. In either case retain the original external handler.
    if (previous != signalHandler) {
      d_prev_handler.store(previous, std::memory_order_relaxed);
    }
  }
  static_assert(std::atomic<SignalHandler>::is_always_lock_free);
  inline static std::atomic_flag d_gotSignal{};
  inline static std::atomic_flag d_handlerTransition{};
  // A concurrent scope must neither replace the original handler nor restore
  // it while another calculation is still using our handler.
#ifdef RDK_BUILD_THREADSAFE_SSS
  inline static std::mutex d_mutex;
#endif
  inline static unsigned int d_handlerCount{0};
  inline static std::atomic<SignalHandler> d_prev_handler{SIG_DFL};
};
}  // namespace RDKit
#endif  // CONTROLCHANDLER_H
