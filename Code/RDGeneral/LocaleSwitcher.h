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

namespace RDKit{
  namespace Utils {
    // allows an RAII-like approach to ensuring the locale is temporarily "C"
    // instead of whatever we started in.
    class LocaleSwitcher {
    public:
      LocaleSwitcher() {
        std::setlocale(LC_ALL,"C");
      }
      ~LocaleSwitcher() {
        std::setlocale(LC_ALL,"");
      }
    };
  }
}

#endif
