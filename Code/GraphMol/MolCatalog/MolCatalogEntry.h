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

  //! This class is used to store ROMol objects in a MolCatalog
  class MolCatalogEntry : public RDCatalog::CatalogEntry {
    public :

    MolCatalogEntry() : dp_mol(0),d_descrip("") {
      dp_props = new Dict();
      setBitId(-1); 
    }

    //! copy constructor
    MolCatalogEntry(const MolCatalogEntry &other);

    //! create an entry to hold the provided ROMol
    /*!
      The MolCatalogEntry takes ownership of the pointer
     */
    MolCatalogEntry(const ROMol *omol);

    //! construct from a pickle
    MolCatalogEntry(const std::string &pickle){
      this->initFromString(pickle);
    }
    
    ~MolCatalogEntry();
    
    std::string getDescription() const {return d_descrip;}
    
    void setDescription(std::string val) {d_descrip = val;}

    unsigned int getOrder() const { return d_order; };
    void setOrder(unsigned int order) { d_order=order; };
    
    const ROMol *getMol() const { return dp_mol; };
    //! hold the provided ROMol
    /*!
      The MolCatalogEntry takes ownership of the pointer.
      If the MolCatalogEntry already has a molecule, this one will be deleted.
     */
    void setMol(const ROMol *molPtr);
    
    //! set a named property
    template <typename T> void setProp(const char *key, T &val) const {
      dp_props->setVal(key, val);
    }

    //! \overload
    template <typename T> void setProp(const std::string key, T &val) const {
      setProp(key.c_str(), val);
    }

    //! get the value of a named property
    template <typename T> 
      void getProp(const char *key, T &res) const {
      dp_props->getVal(key, res);
    }
    //! \overload
    template <typename T>
      void getProp(const std::string key, T &res) const {
      getProp(key.c_str(), res);
    }
    
    //! returns true if such a property exists
    bool hasProp(const char *key) const {
      if (!dp_props) return false;
      return dp_props->hasVal(key);
    }
    //! \overload
    bool hasProp(const std::string key) const {
      return hasProp(key.c_str());
    }
    
    //! clears a named property
    void clearProp(const char *key) const {
      dp_props->clearVal(key);
    }
    //! \overload
    void clearProp(const std::string key) const {
      clearProp(key.c_str());
    }

    //! serializes this entry to the stream
    void toStream(std::ostream &ss) const;
    //! returns a serialized (pickled) form of the entry
    std::string Serialize() const;
    //! initialize from a stream containing a pickle
    void initFromStream(std::istream &ss);
    //! initialize from a string containing a pickle
    void initFromString(const std::string &text);
    
  private:
    const ROMol *dp_mol;
    Dict *dp_props;

    unsigned int d_order;
    std::string d_descrip;
  };
}

#endif
    
