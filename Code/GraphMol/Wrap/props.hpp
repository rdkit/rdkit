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
#ifndef _RD_WRAPPED_PROPS_H_
#define _RD_WRAPPED_PROPS_H_

#include <RDBoost/python.h>
#include <RDBoost/pyint_api.h>
#include <RDBoost/Wrap.h>
#include <RDGeneral/Dict.h>

namespace RDKit
{

template<class T, class U>
bool AddToDict(const U& ob, boost::python::dict &dict, const std::string &key) {
  T res;
  try {
    if (ob.getPropIfPresent(key, res)) {
      dict[key] = res;
    }
  } catch (boost::bad_any_cast &) {
    return false;
  }
  return false;
}

template<class T>
boost::python::dict GetPropsAsDict(const T &obj) {
  boost::python::dict dict;
  // precedence double, int, unsigned, std::vector<double>,
  // std::vector<int>, std::vector<unsigned>, string
  STR_VECT keys = obj.getPropList();
  for(size_t i=0;i<keys.size();++i) {
    if (AddToDict<double>(obj, dict, keys[i])) continue;
    if (AddToDict<int>(obj, dict, keys[i])) continue;
    if (AddToDict<unsigned int>(obj, dict, keys[i])) continue;
    if (AddToDict<bool>(obj, dict, keys[i])) continue;
    if (AddToDict<std::vector<double> >(obj, dict, keys[i])) continue;
    if (AddToDict<std::vector<int> >(obj, dict, keys[i])) continue;
    if (AddToDict<std::vector<unsigned int> >(obj, dict, keys[i])) continue;
    if (AddToDict<std::string>(obj, dict, keys[i])) continue;
  }
  return dict;
}

}

#endif
