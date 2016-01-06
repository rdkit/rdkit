//
//  Copyright (C) 2012 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#ifndef _RD_LOCALESWITCHER_H
#define _RD_LOCALESWITCHER_H

#include <clocale>
#ifndef _MSC_VER
# include <xlocale.h>
#else
# include <locale.h>
#endif

namespace RDKit {
namespace Utils {
// allows an RAII-like approach to ensuring the locale is temporarily "C"
// instead of whatever we started in.
class LocaleSwitcher {
 public:
#ifdef _MSC_VER  

  LocaleSwitcher() {
#ifdef RDK_THREADSAFE_SSS
    _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
    ::setlocale(LC_ALL, "C");
#else    
    std::setlocale(LC_ALL, "C");
#endif
  }
  ~LocaleSwitcher() {
#ifdef RDK_THREADSAFE_SSS
    ::setlocale(LC_ALL, "C");
#else    
    std::setlocale(LC_ALL, "");
#endif
  }
#else
  LocaleSwitcher() {
    std::locale loc;
    if (loc.name() != std::locale::classic().name()) {
      // set locale for this thread
      uselocale(newlocale(LC_ALL, "C", (locale_t)0));
      switched = true;
    }
  }
  ~LocaleSwitcher() {
    if (switched) {
      uselocale(newlocale(LC_ALL, "", (locale_t)0));
    }
  }
  bool switched;
#endif
  
};
}
}

#endif
