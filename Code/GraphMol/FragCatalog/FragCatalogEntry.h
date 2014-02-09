//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_FRAGCATALOGENTRY_H_
#define _RD_FRAGCATALOGENTRY_H_

#include "FragCatParams.h"
#include <RDGeneral/utils.h>
#include <Catalogs/CatalogEntry.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <map>
#include <sstream>

namespace RDKit {
  
  class FragCatalogEntry : public RDCatalog::CatalogEntry {
    public :

    FragCatalogEntry() : dp_mol(0),d_descrip(""),d_order(0) {
      dp_props = new Dict();
      setBitId(-1); 
    }

    FragCatalogEntry(const ROMol *omol, const PATH_TYPE &path, const MatchVectType &aidToFid);
    FragCatalogEntry(const std::string &pickle);
    
    ~FragCatalogEntry() { 
      delete dp_mol;
      dp_mol=0;
      if(dp_props){
	delete dp_props;
	dp_props=0;
      }
    }
    
    std::string getDescription() const {return d_descrip;}
    
    void setDescription(std::string val) {d_descrip = val;}

    void setDescription(const FragCatParams *params);
    
    // check if this fragment macthes the one specified
    //
    
    bool match(const FragCatalogEntry *other, double tol) const;
         
    Subgraphs::DiscrimTuple getDiscrims() const;
    
    unsigned int getOrder() const {
      return dp_mol->getNumBonds();
    }

    const INT_INT_VECT_MAP &getFuncGroupMap() const {
      return d_aToFmap;
    }

    // REVIEW: this should be removed?
    std::string getSmarts(){return "";}

    // FUnctions on the property dictionary
    template <typename T> void setProp(const char *key, T &val) const {
      dp_props->setVal(key, val);
    }

    template <typename T> void setProp(const std::string key, T &val) const {
      setProp(key.c_str(), val);
    }

    void setProp(const char *key, int val) const {
      dp_props->setVal(key, val);
    }

    void setProp(const std::string key, int val) const {
      setProp(key.c_str(), val);
    }

    void setProp(const char *key, float val) const {
      dp_props->setVal(key, val);
    }

    void setProp(const std::string key, float val) const {
      setProp(key.c_str(), val);
    }

    void setProp(const std::string key, std::string val) const {
      setProp(key.c_str(), val);
    }
    
    template <typename T> 
      void getProp(const char *key, T &res) const {
      dp_props->getVal(key, res);
    }
    template <typename T>
      void getProp(const std::string key, T &res) const {
      getProp(key.c_str(), res);
    }
    
    bool hasProp(const char *key) const {
      if (!dp_props) return false;
      return dp_props->hasVal(key);
    }
    bool hasProp(const std::string key) const {
      return hasProp(key.c_str());
    }
    
    void clearProp(const char *key) const {
      dp_props->clearVal(key);
    }

    void clearProp(const std::string key) const {
      clearProp(key.c_str());
    }

    void toStream(std::ostream &ss) const;
    std::string Serialize() const;
    void initFromStream(std::istream &ss);
    void initFromString(const std::string &text);

    
  private:
    
    ROMol *dp_mol;
    Dict *dp_props;

    std::string d_descrip;
    
    unsigned int d_order;

    // a map between the atom ids in mol that connect to 
    // a functional group and the corresponding functional
    // group ID
    INT_INT_VECT_MAP d_aToFmap;
  };
}

#endif
    
