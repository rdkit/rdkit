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

#include <RDGeneral/export.h>
#include <atomic>
#include <csignal>
#include <stdexcept>

namespace RDKit {

class ControlCCaught : public std::runtime_error {
 public:
  explicit ControlCCaught()
      : std::runtime_error("The process was interrupted with Ctrl+c") {};
};

//! This class catches a control-C/SIGINT and sets a flag, read with
//! getGotSignal(), if one is received.  It is intended to be used inside a
//! long C++ calculation called from Python which intercepts the signal
//! handler.  The C++ code must check the value of getGotSignal()
//! periodically and act accordingly.
//! Handlers may be alive on several threads at once, and nested on one: the
//! first to be constructed installs the signal handler, and the last to be
//! destroyed restores the one it replaced and clears the flag, so every
//! calculation running while the signal arrives sees it.
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
class RDKIT_RDGENERAL_EXPORT ControlCHandler {
 public:
  ControlCHandler();
  ControlCHandler(const ControlCHandler &) = delete;
  ControlCHandler(ControlCHandler &&) = delete;
  ControlCHandler &operator=(const ControlCHandler &) = delete;
  ControlCHandler &operator=(ControlCHandler &&) = delete;
  ~ControlCHandler();
  static bool getGotSignal() {
    return d_gotSignal.load(std::memory_order_relaxed);
  }
  static void signalHandler(int signalNumber);
  //! clears the flag and, while any handler is alive, installs the signal
  //! handler again
  static void reset();

 private:
  static std::atomic<bool> d_gotSignal;
};
}  // namespace RDKit
#endif  // CONTROLCHANDLER_H
