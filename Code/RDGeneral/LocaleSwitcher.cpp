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
#include "LocaleSwitcher.h"

// LocaleSwitcher Dependencies
#include <clocale>
#ifndef _WIN32
#include <xlocale.h>
#include <string>
#else
#include <locale.h>
#endif

#ifdef RDK_THREADSAFE_SSS
#include <boost/thread/tss.hpp>
#endif

namespace RDKit {
namespace Utils {
namespace detail {
// Implementation of LocaleSwitcher
//  The locale switcher has state to indicate how many times
//   it has been used (in a thread if thread safe RDKIT)

const static int CurrentState  = 0; // return current state
const static int SwitchLocale  = 1; // indicate we are now switched to "C"
const static int ResetLocale   = -1;// return to previous locale (not "C")

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
#ifndef RDK_THREADSAFE_SSS
  static int recursion = 0;  
  if      (state==SwitchLocale) recursion++;
  else if (state==ResetLocale)  recursion--;
  return recursion;
#else
  static boost::thread_specific_ptr<int> recursion;
  if( ! recursion.get() ) {
    recursion.reset( new int(0));
  }
  if      (state==SwitchLocale) (*recursion)++;
  else if (state==ResetLocale)  (*recursion)--;
  return (*recursion);
#endif
}


// allows an RAII-like approach to ensuring the locale is temporarily "C"
// instead of whatever we started in.
class LocaleSwitcherImpl {
 public:
#ifdef _WIN32

  LocaleSwitcherImpl() {
    if (!recurseLocale(CurrentState)) {
#ifdef RDK_THREADSAFE_SSS
      _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
      ::setlocale(LC_ALL, "C"); // thread safe on windows
#else
      std::setlocale(LC_ALL, "C");
#endif // RDK_THREADSAFE_SSS
      recurseLocale(SwitchLocale);
      switched = true;
    }
  }
  ~LocaleSwitcherImpl() {
    if (switched) {
#ifdef RDK_THREADSAFE_SSS
      ::setlocale(LC_ALL, "C");
#else
      std::setlocale(LC_ALL, ""); // back to last (global?) locale
#endif // RDK_THREADSAFE_SSS
      recurseLocale(ResetLocale);
    }
  }

 public:
  bool switched;
#else // _WIN32
  locale_t loc;     // current "C" locale
  locale_t old_loc; // locale we came frome

  LocaleSwitcherImpl() : old_locale(setlocale(LC_ALL, NULL)) {
    // set locale for this thread

    if (!recurseLocale(CurrentState) && old_locale != "C") {
      recurseLocale(SwitchLocale);
      old_loc = uselocale(0);
      loc = newlocale(LC_ALL_MASK, "C", (locale_t)0);
      uselocale(loc);
      // Don't free "C" or "GLOBAL" Locales
    } else
      old_locale = "C"; // prevents recursion
  }
  ~LocaleSwitcherImpl() {
    if (old_locale != "C") {
      uselocale(old_loc);
      freelocale(loc);
      recurseLocale(ResetLocale);
    }
  }

 public:
  std::string old_locale;

#endif // _WIN32
};
}

// allows an RAII-like approach to ensuring the locale is temporarily "C"
// instead of whatever we started in.

LocaleSwitcher::LocaleSwitcher() : impl(new detail::LocaleSwitcherImpl) {}
LocaleSwitcher::~LocaleSwitcher() { delete impl; }
}
}
