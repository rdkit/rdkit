//
//  Copyright (C) 2006 Greg Landrum
//
#ifndef _RD_MOLCATALOGENTRY_H_
#define _RD_MOLCATALOGENTRY_H_

#include <RDGeneral/Dict.h>
#include <Catalogs/CatalogEntry.h>
#include <fstream>
#include <string>


namespace RDKit {
  class ROMol;

  class MolCatalogEntry : public RDCatalog::CatalogEntry {
    public :

    MolCatalogEntry() : dp_mol(0),d_descrip("") {
      dp_props = new Dict();
      setBitId(-1); 
    }

    MolCatalogEntry(const MolCatalogEntry &other);
    MolCatalogEntry(const ROMol *omol);
    MolCatalogEntry(const std::string &pickle){
      this->initFromString(pickle);
    }
    
    ~MolCatalogEntry();
    
    std::string getDescription() const {return d_descrip;}
    
    void setDescription(std::string val) {d_descrip = val;}

    unsigned int getOrder() const { return d_order; };
    void setOrder(unsigned int order) { d_order=order; };
    
    const ROMol *getMol() const { return dp_mol; };
    void setMol(const ROMol *molPtr);
    
    // Functions on the property dictionary
    template <typename T> void setProp(const char *key, T &val) const {
      dp_props->setVal(key, val);
    }

    template <typename T> void setProp(const std::string key, T &val) const {
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
    const ROMol *dp_mol;
    Dict *dp_props;

    unsigned int d_order;
    std::string d_descrip;
  };
}

#endif
    
