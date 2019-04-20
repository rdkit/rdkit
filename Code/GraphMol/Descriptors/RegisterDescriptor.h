//
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

#include <RDGeneral/export.h>
#ifndef RDKIT_REGISTER_DESCRIPTOR_H
#define RDKIT_REGISTER_DESCRIPTOR_H

#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "Property.h"

namespace RDKit {
class ROMol;
namespace Descriptors {

// these macros create a static class that can call the specified function.
//    using double class::operator(const ROMol&)
//  this class also contains a function pointer with the signature
//      double(*)(const ROMol&)
// these classes are automatically registered with the property registry
/*
#define REGISTER_FULL_DESCRIPTOR( NAME, VERSION, FUNC )         \
double NAME##PropertyFunction(const ROMol&m){return
static_cast<double>(FUNC(mol));}\
struct NAME##PropertyFunctor : public PropertyFunctor{          \
  NAME##PropertyFunctor(bool registerProp=true) : PropertyFunctor(#NAME,
VERSION) { \
    if (registerProp) Properties::registerProperty(new
NAME##PropertyFunctor(false)); \
    d_dataFunc = &NAME##PropertyFunction;  \
  } \
  double operator()(const RDKit::ROMol &mol) const { \
     return NAME##PropertyFunction(mol); } \
}; \
static NAME##PropertyFunctor NAME##PropertyFunctor__;
*/

#define REGISTER_DESCRIPTOR(NAME, FUNC)                                     \
  struct NAME##PropertyFunctor : public PropertyFunctor {                   \
    static double _func(const ROMol &m) {                                   \
      return static_cast<double>(FUNC(m));                                  \
    }                                                                       \
    NAME##PropertyFunctor(bool registerProp = true)                         \
        : PropertyFunctor(#NAME, NAME##Version, _func) {                    \
      if (registerProp)                                                     \
        Properties::registerProperty(new NAME##PropertyFunctor(false));     \
    }                                                                       \
    double operator()(const RDKit::ROMol &mol) const { return _func(mol); } \
  };                                                                        \
  static NAME##PropertyFunctor NAME##PropertyFunctor__;
}  // namespace Descriptors
}  // namespace RDKit

#endif
