// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <fuzzer/FuzzedDataProvider.h>

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

// fuzz_target.cc
extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
  FuzzedDataProvider fdp(Data, Size);

  const bool sanitize = fdp.ConsumeIntegralInRange(0, 1);
  std::string smiles_string = fdp.ConsumeRemainingBytesAsString();
  try {
    std::shared_ptr<RDKit::ROMol> mol3(
        RDKit::SmilesToMol(smiles_string, 0, sanitize));
  } catch (...) {
  }

  return 0;  // Non-zero return values are reserved for future use.
}