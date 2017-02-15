#ifndef RDKIT_RDPROPS_H
#define RDKIT_RDPROPS_H
#include "Dict.h"
#include "types.h"
#include <boost/foreach.hpp>

namespace RDKit {

class RDProps {
 protected:
  mutable Dict dp_props;
  // It is a quirk of history that this is mutable
  //  as the RDKit allows properties to be set
  //  on const obects.

 public:
  RDProps() : dp_props() {}
  RDProps(const RDProps &rhs) : dp_props(rhs.dp_props) {}
  RDProps &operator=(const RDProps &rhs) {
    dp_props = rhs.dp_props;
    return *this;
  }
  void clear() { dp_props.reset(); }
  //! gets the underlying Dictionary
  const Dict &getDict() const { return dp_props; }
  Dict &getDict() { return dp_props; }

  // ------------------------------------
  //  Local Property Dict functionality
  //  all setProp functions are const because they
  //     are not meant to change the atom chemically
  // ------------------------------------
  //! returns a list with the names of our \c properties
  STR_VECT getPropList(bool includePrivate = true,
                       bool includeComputed = true) const {
    const STR_VECT &tmp = dp_props.keys();
    STR_VECT res, computed;
    if (!includeComputed &&
        getPropIfPresent(RDKit::detail::computedPropName, computed)) {
      computed.push_back(RDKit::detail::computedPropName);
    }

    STR_VECT::const_iterator pos = tmp.begin();
    while (pos != tmp.end()) {
      if ((includePrivate || (*pos)[0] != '_') &&
          std::find(computed.begin(), computed.end(), *pos) == computed.end()) {
        res.push_back(*pos);
      }
      pos++;
    }
    return res;
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
  void setProp(const std::string &key, T val, bool computed = false) const {
    if (computed) {
      STR_VECT compLst;
      getPropIfPresent(RDKit::detail::computedPropName, compLst);
      if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        compLst.push_back(key);
        dp_props.setVal(RDKit::detail::computedPropName, compLst);
      }
    }
    // setProp(key.c_str(),val);
    dp_props.setVal(key, val);
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
  void getProp(const std::string &key, T &res) const {
    dp_props.getVal(key, res);
  }

  //! \overload
  template <typename T>
  T getProp(const std::string &key) const {
    return dp_props.getVal<T>(key);
  }

  //! returns whether or not we have a \c property with name \c key
  //!  and assigns the value if we do
  //! \overload
  template <typename T>
  bool getPropIfPresent(const std::string &key, T &res) const {
    return dp_props.getValIfPresent(key, res);
  }

  //! \overload
  bool hasProp(const std::string &key) const { return dp_props.hasVal(key); };

  //! clears the value of a \c property
  /*!
    <b>Notes:</b>
    - if no \c property with name \c key exists, a KeyErrorException
    will be thrown.
    - if the \c property is marked as \c computed, it will also be removed
    from our list of \c computedProperties
  */
  //! \overload
  void clearProp(const std::string &key) const {
    STR_VECT compLst;
    if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
      STR_VECT_I svi = std::find(compLst.begin(), compLst.end(), key);
      if (svi != compLst.end()) {
        compLst.erase(svi);
        dp_props.setVal(RDKit::detail::computedPropName, compLst);
      }
    }
    dp_props.clearVal(key);
  };

  //! clears all of our \c computed \c properties
  void clearComputedProps() const {
    STR_VECT compLst;
    if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
      BOOST_FOREACH (const std::string &sv, compLst) { dp_props.clearVal(sv); }
      compLst.clear();
      dp_props.setVal(RDKit::detail::computedPropName, compLst);
    }
  }
};
}
#endif
