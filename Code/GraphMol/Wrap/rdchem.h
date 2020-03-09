//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RDCHEM_INCL_
#define _RDCHEM_INCL_

#define PY_ARRAY_UNIQUE_SYMBOL rdchem_array_API

namespace RDKit {
class ConformerException;
}
void rdExceptionTranslator(RDKit::ConformerException const& x);

#endif
