//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#if defined(__CYGWIN__) && !defined(_GNU_SOURCE)
// -std=c++11 turns off recent POSIX features!
#define _GNU_SOURCE
#endif

#include "LocaleSwitcher.h"
#include <string>

// LocaleSwitcher Dependencies
#include <clocale>
#ifdef __APPLE__
#include <xlocale.h>
#endif
#ifdef __CYGWIN__
#include <locale.h>
#endif
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#endif

namespace RDKit {
namespace Utils {
namespace detail {
// Implementation of LocaleSwitcher
//  The locale switcher has state to indicate how many times
//   it has been used (in a thread if thread safe RDKIT)

const static int CurrentState = 0;  // return current state
const static int SwitchLocale = 1;  // indicate we are now switched to "C"
const static int ResetLocale = -1;  // return to previous locale (not "C")

// Returns how many LocaleSwitches we have gone down.
//  Such as
//    void foo() { LocaleSwitcher switch; }
//    void bar() { LocaleSwitcher switch; foo(); }
//
//  Implementations are free to not switch the locale if the locale
//   is already in "C", but this is implementation defined.
//  When Recurse(CurrentState) == 0 we are in the "global" locale,
//   or at least at the top locale that isn't "C".
//
static int recurseLocale(int state) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static thread_local int recursion = 0;
#else
  static int recursion = 0;
#endif
  if (state == SwitchLocale) {
    ++recursion;
  } else if (state == ResetLocale) {
    --recursion;
  }
  return recursion;
}

// allows an RAII-like approach to ensuring the locale is temporarily "C"
// instead of whatever we started in.
class LocaleSwitcherImpl {
 public:
#ifdef _WIN32

  LocaleSwitcherImpl() {
    if (!recurseLocale(CurrentState)) {
#ifdef RDK_BUILD_THREADSAFE_SSS
      _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
      old_locale = ::setlocale(LC_ALL, nullptr);
      ::setlocale(LC_ALL, "C");  // thread safe on windows
#else
      old_locale = std::setlocale(LC_ALL, nullptr);
      std::setlocale(LC_ALL, "C");
#endif  // RDK_BUILD_THREADSAFE_SSS
      recurseLocale(SwitchLocale);
      switched = true;
    }
  }
  ~LocaleSwitcherImpl() {
    if (switched) {
#ifdef RDK_BUILD_THREADSAFE_SSS
      ::setlocale(LC_ALL, old_locale.c_str());
#else
      // back to last (global?) locale
      std::setlocale(LC_ALL, old_locale.c_str());
#endif  // RDK_BUILD_THREADSAFE_SSS
      recurseLocale(ResetLocale);
    }
  }

 public:
  bool switched = false;
  std::string old_locale;
#else  // _WIN32
  LocaleSwitcherImpl() {
    // set locale for this thread

    if (!recurseLocale(CurrentState)) {
      auto loc = newlocale(LC_ALL_MASK, "C", (locale_t) nullptr);
      old_loc = uselocale(loc);
      recurseLocale(SwitchLocale);
      switched = true;
    }
  }
  ~LocaleSwitcherImpl() {
    if (switched) {
      auto loc = uselocale(old_loc);
      freelocale(loc);
      recurseLocale(ResetLocale);
    }
  }

 private:
  bool switched = false;
  locale_t old_loc;  // locale we came from

#endif  // _WIN32
};
}  // namespace detail

// allows an RAII-like approach to ensuring the locale is temporarily "C"
// instead of whatever we started in.

LocaleSwitcher::LocaleSwitcher() : impl(new detail::LocaleSwitcherImpl) {}
LocaleSwitcher::~LocaleSwitcher() { delete impl; }
}  // namespace Utils
}  // namespace RDKit
