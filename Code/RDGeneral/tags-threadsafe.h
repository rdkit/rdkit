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
#ifndef RDKIT_TAGS_THREADSAFE_H
#define RDKIT_TAGS_THREADSAFE_H

#include "BoostStartInclude.h"
#include <boost/unordered_map.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#include "BoostEndInclude.h"

#include "types.h"

namespace RDKit
{

class RDTags
{
  RDTags(const RDTags&);
  RDTags& operator=(const RDTags&);
public:
  RDTags() {
    // Load the known keys
    for(int i=0;i<=common_properties::MAX;++i) {
      keyIntMap[common_properties::getPropName(i)] = i;
      keys.push_back(common_properties::getPropName(i));
    }
  };

  //! \brief Get the current string key for the integer key.
  //!  throws std::out_of_range when key is out
  //!  of range.
  inline const std::string &get(int tag) const {
    return keys[tag]; // raises exception
  }

  //! \brief Return the int key for the specified string key
  //! If the string key has never been used before, it is
  //!  added to the global map and a new integer key is created.
  inline int get(const std::string &key) const
  {
    // local threadsafe map used for quick non-contested lookups in threads
    static boost::thread_specific_ptr<boost::unordered_map<std::string, int> >
        tls_keyIntMap;

    if ( !tls_keyIntMap.get() ) {
      // enable TLS loaded with current map
      boost::unique_lock<boost::mutex> lock(_globalMapMutex);
      tls_keyIntMap.reset(new boost::unordered_map<std::string, int>(keyIntMap));
    }
    
    boost::unordered_map<std::string, int> &local_keyIntMap= *tls_keyIntMap;
    // check the thread local cache
    {
      boost::unordered_map<std::string, int>::const_iterator it=local_keyIntMap.find(key);
      if (it != local_keyIntMap.end())
        return it->second;
    }
    
    {
      boost::unique_lock<boost::mutex> lock(_globalMapMutex);

      // Check to see if the key has been added to the global map
      boost::unordered_map<std::string, int>::const_iterator it=keyIntMap.find(key);
        
      if (it != keyIntMap.end()) {
        // if so, insert into local map and return
        local_keyIntMap.insert(*it);
        return it->second;
      }

      // otherwise make a new map and add it to the local and global maps
      int res = static_cast<int>(keyIntMap.size());
      std::pair<std::string, int> p(key,res);
      local_keyIntMap.insert(p);
      keyIntMap.insert(p);
      keys.push_back(key);
      return res;
      // lock is released
    }
  }
  
private:
    mutable boost::mutex _globalMapMutex;
public:
    mutable boost::unordered_map<std::string, int> keyIntMap;
    mutable std::vector<std::string > keys;
};

}

#endif
