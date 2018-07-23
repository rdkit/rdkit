//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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

#include <RDGeneral/export.h>
#ifndef __RDKIT_FILTERDATA__
#define __RDKIT_FILTERDATA__

#include "FilterCatalogEntry.h"
#include "FilterCatalog.h"

namespace RDKit {
struct FilterData_t {
  const char *name;
  const char *smarts;
  unsigned int max;
  const char *comment;
};

struct FilterProperty_t {
  const char *key;
  const char *value;
};

unsigned int GetNumEntries(FilterCatalogParams::FilterCatalogs catalog);
const FilterData_t *GetFilterData(FilterCatalogParams::FilterCatalogs catalog);
unsigned int GetNumPropertyEntries(FilterCatalogParams::FilterCatalogs catalog);
const FilterProperty_t *GetFilterProperties(
    FilterCatalogParams::FilterCatalogs catalog);

FilterCatalogEntry *MakeFilterCatalogEntry(const FilterData_t &,
                                           unsigned int num_props = 0,
                                           const FilterProperty_t *props = 0);
}

#endif
