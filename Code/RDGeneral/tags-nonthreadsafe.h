/* 
*
*  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef RDKIT_TAGS_NONTHREADSAFE_H
#define RDKIT_TAGS_NONTHREADSAFE_H

#include "types.h"

#include "BoostStartInclude.h"
#include <boost/unordered_map.hpp>
#include "BoostEndInclude.h"


namespace RDKit
{

// RDTags is designed to be used as a static singleton
class RDTags
{
  RDTags(const RDTags&);
  RDTags& operator=(const RDTags&);
public:
  // Load the known keys
  RDTags() {
    for(int i=0;i<=common_properties::MAX;++i) {
      keyIntMap[common_properties::getPropName(i)] = i;
      keys.push_back(common_properties::getPropName(i));
    }
  };

  //! \brief Get the current string key for the integer key.
  //!  throws std::out_of_range when key is out
  //!  of range.  
  const std::string &get(int tag) const {
    return keys[tag]; // raises exception
  }

  //! \brief Return the int key for the specified string key
  //! If the string key has never been used before, it is
  //!  added to the global map and a new integer key is created.
  int get(const std::string &key) const
  {
    boost::unordered_map<std::string, int>::const_iterator it=keyIntMap.find(key);
    
    if (it != keyIntMap.end()) {
      return it->second;
    }

    int res = static_cast<int>(keyIntMap.size());
    std::pair<std::string, int> p(key,res);
    keyIntMap.insert(p);
    keys.push_back(key);
    return res;
  }
  
public:
  mutable boost::unordered_map<std::string, int> keyIntMap;
  mutable std::vector<std::string > keys;
};

}

#endif
