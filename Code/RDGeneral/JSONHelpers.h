#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <array>

// Convert a JSON object with boolean keys to a bitwise enum flag set.
// The T enum must have been declared as a BETTER_ENUM.
// The JSON object should have keys corresponding to the enum names, and boolean
// values. Example JSON: { "FlagA": true, "FlagB": false, "FlagC": true,
// "FlagD": false } This would set FlagA and FlagC in the resulting flag set,
// then clear FlagB and FlagD. If keysFound is provided and non-null, it will be
// set to true if any of the enum keys were found in the JSON object, false
// otherwise.
template <typename T>
typename T::_integral flagsFromJson(const boost::property_tree::ptree &pt,
                                    bool *keysFound = nullptr) {
  std::array<typename T::_integral, 2> flagsByType;
  flagsByType.fill(static_cast<T::_integral>(0));
  if (keysFound) {
    *keysFound = false;
  }
  for (const auto *key : T::_names()) {
    const auto it = pt.find(key);
    if (it == pt.not_found()) {
      continue;
    }
    if (keysFound) {
      *keysFound = true;
    }
    auto value = T::_from_string(key)._to_integral();
    auto i = static_cast<unsigned int>(it->second.template get_value<bool>());
    flagsByType[i] |= value;
  }
  return (flagsByType[1] & ~flagsByType[0]);
}
