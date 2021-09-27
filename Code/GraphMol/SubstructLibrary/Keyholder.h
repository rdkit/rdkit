//
//  Copyright (C) 2021 Brian P Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDK_SUBSTRUCT_KEYHOLDER
#define RDK_SUBSTRUCT_KEYHOLDER
#include <RDGeneral/export.h>
#include <algorithm>
#include <string>
#include <vector>
#include <GraphMol/RDKitBase.h>

namespace RDKit {
class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyHolderBase {
public:
  virtual ~KeyHolderBase() {}

  //! Add a key to the database getting it from the molecule
  virtual unsigned int addMol(const ROMol &m) = 0;
  
  //! Add a key to the database, this needs to be in the same order
  //!  as the molecule, no validation is done
  virtual unsigned int addKey(const std::string &) = 0;

  // !get the key at the requested index
  // implementations should throw IndexError on out of range
  virtual const std::string & getKey(unsigned int) = 0;

  // !get keys from a bunch of indices
  virtual std::vector<std::string> getKeys(const std::vector<unsigned int> &indices) const = 0;
  //! Get the current keeyholder size
  virtual unsigned int size() const = 0;
};

class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyFromPropHolder : public KeyHolderBase {
  std::string propname;
  std::vector<std::string> keys;
  const std::string empty_string = {};
  
public:
  KeyFromPropHolder(const std::string &propname = "_Name") : propname(propname) {
  }

  std::string &getPropName() { return propname; }  
  const std::string &getPropName() const { return propname; }  

  std::vector<std::string> &getKeys() { return keys; }  
  const std::vector<std::string> &getKeys() const { return keys; }  

  unsigned int addMol(const ROMol &m) override {
    std::string key;
    if (m.getPropIfPresent(propname, key)) {
      keys.push_back(std::move(key));
    } else {
      // is this a warning?  Should we push back the numeric index?
      std::cerr << "Molecule at index " << keys.size() << " doesn't have propname " << propname << std::endl;
      keys.emplace_back("");
    }
    return keys.size() - 1u;
  };

  unsigned int addKey(const std::string &key) override {
    keys.push_back(key);
    return keys.size() - 1u;
  }

  const std::string & getKey(unsigned int idx) override {
    if (idx >= keys.size()) throw IndexErrorException(idx);
    return keys[idx];
  }
  
  std::vector<std::string> getKeys(const std::vector<unsigned int> &indices) const override{
    std::vector<std::string> res;
    std::transform(indices.begin(), indices.end(), std::back_inserter(res),
		   [=](unsigned idx){return keys.at(idx);});
    return res;
  }
  unsigned int size() const override {
    return keys.size();
  }
  
};
}
#endif 
