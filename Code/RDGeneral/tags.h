/* 
*
*  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
#ifndef RDKIT_TAGS_H
#define RDKIT_TAGS_H

#include "BoostStartInclude.h"

#ifdef RDK_TEST_MULTITHREADED
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#endif

#include <boost/unordered_map.hpp>
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
    for(int i=0;i<=common_properties::MAX;++i) {
      m[common_properties::getPropName(i)] = i;
      keys.push_back(common_properties::getPropName(i));
    }
  };

  const std::string &get(int tag) const {
    return keys[tag]; // raises exception
  }
  
  int get(const std::string &k) const
  {
#ifdef RDK_TEST_MULTITHREADED
    static boost::thread_specific_ptr<boost::unordered_map<std::string, int> >
        instance;

    if ( !instance.get() ) {
      boost::unique_lock<boost::mutex> lock(_m);
      // enable TLS loaded with current map
      instance.reset(new boost::unordered_map<std::string, int>(m));
    }
    
    boost::unordered_map<std::string, int> &ptr = *instance;
    // check tls map
    {
      boost::unordered_map<std::string, int>::const_iterator it=ptr.find(k);
      if (it != ptr.end())
        return it->second;
    }
#endif
    
    {
      boost::unique_lock<boost::mutex> lock(_m);
      
      {
        boost::unordered_map<std::string, int>::const_iterator it=m.find(k);
        
        if (it != m.end()) {
#ifdef RDK_TEST_MULTITHREADED          
          // insert into local map        
          instance->insert(*it);
#endif
          return it->second;
        }
      }
      {
        int res = m.size();
        std::pair<std::string, int> p(k,res);
#ifdef RDK_TEST_MULTITHREADED        
        instance->insert(p);
#endif        
        m.insert(p);
        keys.push_back(k);
        return res;
      }
    }
  }
  
private:
    mutable boost::mutex _m;
public:
    mutable boost::unordered_map<std::string, int> m;
    mutable std::vector<std::string > keys;
};

}

#endif
