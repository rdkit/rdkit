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

#include <csignal>

#include <RDGeneral/export.h>

namespace RDKit {
//! This class catches a control-C/SIGINT and sets the flag d_gotSignal
//! if one is received.  It is intended to be used inside a long
//! C++ calculation called from Python which intercepts the signal
//! handler.  The C++ code must check the value of d_gotSignal
//! periodically and act accordingly.  The destructor resets
//! the signal handler and flag for next use, which is essential
//! because it's a static variable.
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
  ControlCHandler() { d_prev_handler = std::signal(SIGINT, signalHandler); }
  ~ControlCHandler() {
    std::signal(SIGINT, d_prev_handler);
    d_gotSignal = false;
  }
  static bool getGotSignal() { return d_gotSignal; }
  static void setGotSignal(bool newVal) { d_gotSignal = newVal; }
  static void signalHandler(int signalNumber) {
    if (signalNumber == SIGINT) {
      d_gotSignal = true;
    }
  }

 private:
  inline static volatile std::atomic<bool> d_gotSignal{false};
  void (*d_prev_handler)(int);
};
}  // namespace RDKit
#endif  // CONTROLCHANDLER_H
