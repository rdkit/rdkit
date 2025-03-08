// Reimplementation of missing csInitMap

// BSD 3-Clause License
// 
// Copyright (c) 2025, Glysade Inc
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

namespace cs {
  template<typename Key, typename Value>
  class csInitMap {
  public:
    std::map<Key,Value> storage;
    typedef std::pair<Key,Value> pair;
    csInitMap() : storage() {}

    typename std::map<Key,Value>::iterator find(const Key &key) {
      return storage.find(key);
    }
    Value & operator[](const Key&key) {
      return storage[key];
    }
    typename std::map<Key,Value>::iterator end() { return storage.end(); }
    
    csInitMap<Key,Value>& operator<<(const pair &) {
      return *this;
    }
  //private:
  //  friend csInitMap<Key,Value> & operator<<(csInitMap<Key,Value> &map, const csInitMap<Key,Value>::pair& p);
  };

//template<typename Key, typename Value>
//csInitMap<Key,Value>& operator<<(csInitMap<Key,Value>& os, const std::pair<Key, Value>& p) {
//    return os;
//}
}
