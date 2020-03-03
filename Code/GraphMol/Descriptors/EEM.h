//
//  Copyright (c) 2017, Guillaume GODIN
//  "Copyright 2013-2016 Tomas Racek (tom@krab1k.net)"
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

#ifndef EEMRDKIT_H_SEPT2017
#define EEMRDKIT_H_SEPT2017

#ifdef RDK_BUILD_DESCRIPTORS3D
namespace RDKit {
class ROMol;
namespace Descriptors {

namespace {
class EEM_arrays {
 public:
  unsigned int n;
  unsigned int *Atomindex;
  unsigned int *EEMatomtype;

  EEM_arrays() = delete;
  EEM_arrays(const EEM_arrays &) = delete;
  void operator=(const EEM_arrays &) = delete;

  EEM_arrays(const ROMol &mol, unsigned int n);
  ~EEM_arrays();
};
}  // namespace

const std::string EEMVersion = "1.0.0";
void RDKIT_DESCRIPTORS_EXPORT EEM(ROMol &mol, std::vector<double> &res,
                                  int confId);
}  // namespace Descriptors
}  // namespace RDKit
#endif
#endif
