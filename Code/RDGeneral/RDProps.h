#ifndef RDKIT_RDPROPS_H
#define RDKIT_RDPROPS_H
#include "Dict.h"
#include "types.h"
#include <boost/foreach.hpp>

namespace RDKit {

class RDProps {
protected:
  // It is a quirk of RDKit history that this is mutable
  //  as the RDKit allows properties to be set
  //  on const obects.

  mutable Dict dp_props;              // Holds the properties
  mutable std::vector<int> compLst;   // Holds the tags that are computed props

 public:
  RDProps() : dp_props() {}
  RDProps(const RDProps &rhs) : dp_props(rhs.dp_props) {}
  RDProps& operator=(const RDProps &rhs) {
    if (this != &rhs) {
      dp_props = rhs.dp_props;
      compLst = rhs.compLst;
    }
    return *this;
  }
  void clear() {
    dp_props.reset();
    compLst.clear();
  }
  // ------------------------------------
  //  Local Property Dict functionality
  //  all setProp functions are const because they
  //     are not meant to change the atom chemically
  // ------------------------------------
  //! returns a list with the names of our \c properties
  STR_VECT getPropList(bool includePrivate = true,
                       bool includeComputed = true) const {
    STR_VECT res;
    for(Dict::const_iterator it=dp_props.begin();
        it != dp_props.end();
        ++it) {
      int key = it->getKey();
      const char *name = common_properties::getPropName(key);
      CHECK_INVARIANT(name, "Bad tag value in dictionary!!!");
      
      if (!includePrivate && name[0] == '_') continue;
      if (!includeComputed && isComputedProp(key))
        continue;
      res.push_back(name);
    }
    return res;
    }

  STR_VECT getComputedPropList() const {
    STR_VECT res;
    for(Dict::const_iterator it=dp_props.begin();
        it != dp_props.end();
        ++it) {
      if (isComputedProp(it->getKey())) {
        res.push_back(common_properties::getPropName(it->getKey()));
      }
    }
    return res;
  }
  
  bool isComputedProp(int tag) const {
    return std::find(compLst.begin(), compLst.end(), tag) != compLst.end();
  }
  
  bool isComputedProp(const std::string &prop) const {
    return isComputedProp(GetKey(prop));
  }
  //! sets a \c property value
  /*!
    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param val the value to be stored
    \param computed (optional) allows the \c property to be flagged
    \c computed.
  */

  //! \overload
  template <typename T>
  void setProp(int tag, T val, bool computed = false) const {
    if (computed && !isComputedProp(tag)) {
        compLst.push_back(tag);
      }

    dp_props.setVal(tag, val);
  }

  template <typename T>
  void setProp(const std::string &key, T val, bool computed = false) const {
    setProp(GetKey(key), val, computed);
  }
  //! allows retrieval of a particular property value
  /*!

    \param key the name under which the \c property should be stored.
    If a \c property is already stored under this name, it will be
    replaced.
    \param res a reference to the storage location for the value.

    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException will be
    thrown.
    - the \c boost::lexical_cast machinery is used to attempt type
    conversions.
    If this fails, a \c boost::bad_lexical_cast exception will be thrown.

  */
  //! \overload
  template <typename T>
  void getProp(int tag, T &res) const {
    dp_props.getVal(tag, res);
  }
  template <typename T>
  void getProp(const std::string &key, T &res) const {
    return getProp(GetKey(key), res);
  }

  //! \overload
  template <typename T>
  T getProp(int tag) const {
    return dp_props.getVal<T>(tag);
  }
  template <typename T>
  T getProp(const std::string &key) const {
    return getProp<T>(GetKey(key));
  }
  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  //! \overload
  template <typename T>
  bool getPropIfPresent(int tag, T &res) const {
    return dp_props.getValIfPresent(tag, res);
  }

  template <typename T>
  bool getPropIfPresent(const std::string &key, T &res) const {
    return getPropIfPresent(GetKey(key), res);
  }
  //! \overload
  bool hasProp(int tag) const {
    return dp_props.hasVal(tag);
  }
  
  bool hasProp(const std::string &key) const {
    return dp_props.hasVal(key);
  };

  //! clears the value of a \c property
  /*!
    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException
    will be thrown.
    - if the \c property is marked as \c computed, it will also be removed
    from our list of \c computedProperties
  */
  //! \overload
  void clearProp(int tag) const {
    std::vector<int>::iterator it = std::find(compLst.begin(),
                                              compLst.end(),
                                              tag);
    if (it != compLst.end())
      compLst.erase(it);
    dp_props.clearVal(tag);
  };

  void clearProp(const std::string &key) const {
    clearProp(GetKey(key));
  };

  //! clears all of our \c computed \c properties
  void clearComputedProps() const {
    for(std::vector<int>::iterator it=compLst.begin(); it!=compLst.end(); ++it)
      dp_props.clearVal(*it);
      compLst.clear();
  }
};
}
#endif
