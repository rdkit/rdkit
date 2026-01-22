//  Copyright (C) 2026, Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_WRAPPED_PROPS_H
#define RD_WRAPPED_PROPS_H

#include <nanobind/nanobind.h>

#include <RDGeneral/Dict.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace nb = nanobind;

namespace RDKit {

template <class T>
inline const char *GetTypeName() {
  // PRECONDITION(0, "Unregistered c++ type");
  return "unregistered C++ type";
}

template <>
inline const char *GetTypeName<double>() {
  return "a double value";
}
template <>
inline const char *GetTypeName<int>() {
  return "an integer value";
}
template <>
inline const char *GetTypeName<unsigned int>() {
  return "an unsigned integer value";
}
template <>
inline const char *GetTypeName<bool>() {
  return "a True or False value";
}

template <class T, class U>
bool AddToDict(const U &ob, nb::dict &dict, const std::string &key) {
  T res;
  try {
    if (ob.getPropIfPresent(key, res)) {
      dict[key] = res;
    }
  } catch (std::bad_any_cast &) {
    return false;
  }
  return true;
}

const std::string getPropsAsDictDocString =
    "Returns a dictionary populated with properties.\n"
    "When possible, string values will be converted to integers or doubles (trimming if necessary)\n"

    " n.b. Some properties are not able to be converted to python "
    "types.\n\n"
    "  ARGUMENTS:\n"
    "    - includePrivate: (optional) toggles inclusion of private "
    "properties in the result set.\n"
    "                      Defaults to False.\n"
    "    - includeComputed: (optional) toggles inclusion of computed "
    "properties in the result set.\n"
    "                      Defaults to False.\n\n"
    "    - autoConvertStrings: (optional) toggles automatic conversion of string "
    "properties to integers or doubles.\n"
    "                      Defaults to True.\n\n"
    "  RETURNS: a dictionary\n";

template <class T>
nb::dict GetPropsAsDict(const T &obj, bool includePrivate, bool includeComputed,
                        bool autoConvertStrings = true) {
  nb::dict dict;
  auto &rd_dict = obj.getDict();
  auto &data = rd_dict.getData();

  STR_VECT keys = obj.getPropList(includePrivate, includeComputed);
  for (auto &rdvalue : data) {
    if (std::find(keys.begin(), keys.end(), rdvalue.key) == keys.end())
      continue;
    try {
      const auto tag = rdvalue.val.getTag();
      switch (tag) {
        case RDTypeTag::IntTag:
          dict[rdvalue.key] = from_rdvalue<int>(rdvalue.val);
          break;
        case RDTypeTag::DoubleTag:
          dict[rdvalue.key] = from_rdvalue<double>(rdvalue.val);
          break;
        case RDTypeTag::StringTag: {
          auto value = from_rdvalue<std::string>(rdvalue.val);
          if (autoConvertStrings) {
            auto trimVal = value;
            boost::trim(trimVal);
            // Auto convert strings to ints and double if possible
            int ivalue;
            if (boost::conversion::try_lexical_convert(trimVal, ivalue)) {
              dict[rdvalue.key] = ivalue;
              break;
            }
            double dvalue;
            if (boost::conversion::try_lexical_convert(trimVal, dvalue)) {
              dict[rdvalue.key] = dvalue;
              break;
            }
          }
          dict[rdvalue.key] = value;
        } break;
        case RDTypeTag::FloatTag:
          dict[rdvalue.key] = from_rdvalue<float>(rdvalue.val);
          break;
        case RDTypeTag::BoolTag:
          dict[rdvalue.key] = from_rdvalue<bool>(rdvalue.val);
          break;
        case RDTypeTag::UnsignedIntTag:
          dict[rdvalue.key] = from_rdvalue<unsigned int>(rdvalue.val);
          break;
        case RDTypeTag::AnyTag:
          // we skip these for now
          break;
        case RDTypeTag::VecDoubleTag:
          dict[rdvalue.key] = from_rdvalue<std::vector<double>>(rdvalue.val);
          break;
        case RDTypeTag::VecFloatTag:
          dict[rdvalue.key] = from_rdvalue<std::vector<float>>(rdvalue.val);
          break;
        case RDTypeTag::VecIntTag:
          dict[rdvalue.key] = from_rdvalue<std::vector<int>>(rdvalue.val);
          break;
        case RDTypeTag::VecUnsignedIntTag:
          dict[rdvalue.key] =
              from_rdvalue<std::vector<unsigned int>>(rdvalue.val);
          break;
        case RDTypeTag::VecStringTag:
          dict[rdvalue.key] =
              from_rdvalue<std::vector<std::string>>(rdvalue.val);
          break;
        case RDTypeTag::EmptyTag:
          dict[rdvalue.key] = nb::none();
          break;
        default:
          std::string message =
              std::string(
                  "Unhandled property type encountered for property: ") +
              rdvalue.key;
          UNDER_CONSTRUCTION(message.c_str());
      }
    } catch (std::bad_any_cast &) {
      // C++ datatypes can really be anything, this just captures mislabelled
      // data, it really shouldn't happen
      std::string message =
          std::string("Unhandled type conversion occured for property: ") +
          rdvalue.key;
      UNDER_CONSTRUCTION(message.c_str());
    }
  }
  return dict;
}
#if 0
static PyObject *rawPy(python::object &&pyobj) {
  Py_INCREF(pyobj.ptr());
  return pyobj.ptr();
}

template <class T>
PyObject *rawPy(T &&thing) {
  return rawPy(python::object(thing));
}
#endif

template <class RDOb, class T>
nb::object GetProp(const RDOb *ob, const std::string &key) {
  T res;
  try {
    ob->getProp(key, res);
  } catch (const KeyErrorException &) {
    auto msg = std::string("key `") + key + "` not found";
    throw nb::key_error(msg.c_str());
  } catch (const std::exception &e) {
    auto msg = std::string("key `") + key + "` exists but does not result in " +
               GetTypeName<T>() + " reason: " + e.what();
    throw nb::value_error(msg.c_str());
  }
  return nb::object(res);
}

template <class RDOb>
python::object autoConvertString(const RDOb *ob, const std::string &key) {
  int ivalue;
  double dvalue;
  std::string svalue;

  if (ob->getPropIfPresent(key, ivalue))
    return python::object(ivalue);
  else if (ob->getPropIfPresent(key, dvalue))
    return python::object(dvalue);
  else if (ob->getPropIfPresent(key, svalue))
    return python::object(svalue);

  return python::object();
}

template <class RDOb>
PyObject *GetPyProp(const RDOb *obj, const std::string &key, bool autoConvert) {
  python::object pobj;
  if (!autoConvert) {
    std::string res;
    if (obj->getPropIfPresent(key, res)) {
      return rawPy(res);
    } else {
      PyErr_SetString(PyExc_KeyError, key.c_str());
      return nullptr;
    }
  } else {
    const auto &data = obj->getDict().getData();
    for (auto &rdvalue : data) {
      if (rdvalue.key == key) {
        try {
          const auto tag = rdvalue.val.getTag();
          switch (tag) {
            case RDTypeTag::IntTag:
              return rawPy(from_rdvalue<int>(rdvalue.val));

            case RDTypeTag::DoubleTag:
              return rawPy(from_rdvalue<double>(rdvalue.val));

            case RDTypeTag::StringTag:
              if (autoConvert) {
                pobj = autoConvertString(obj, rdvalue.key);
              }
              return rawPy(from_rdvalue<std::string>(rdvalue.val));
            case RDTypeTag::FloatTag:
              return rawPy(from_rdvalue<float>(rdvalue.val));
              break;
            case RDTypeTag::BoolTag:
              return rawPy(from_rdvalue<bool>(rdvalue.val));
              break;
            case RDTypeTag::UnsignedIntTag:
              return rawPy(from_rdvalue<unsigned int>(rdvalue.val));
              break;
            case RDTypeTag::AnyTag:
              // we skip these for now
              break;
            case RDTypeTag::VecDoubleTag:
              return rawPy(from_rdvalue<std::vector<double>>(rdvalue.val));
              break;
            case RDTypeTag::VecFloatTag:
              return rawPy(from_rdvalue<std::vector<float>>(rdvalue.val));
              break;
            case RDTypeTag::VecIntTag:
              return rawPy(from_rdvalue<std::vector<int>>(rdvalue.val));
              break;
            case RDTypeTag::VecUnsignedIntTag:
              return rawPy(
                  from_rdvalue<std::vector<unsigned int>>(rdvalue.val));
              break;
            case RDTypeTag::VecStringTag:
              return rawPy(from_rdvalue<std::vector<std::string>>(rdvalue.val));
              break;
            case RDTypeTag::EmptyTag:
              return Py_None;
              break;
            default:
              std::string message =
                  std::string(
                      "Unhandled property type encountered for property: ") +
                  rdvalue.key;
              UNDER_CONSTRUCTION(message.c_str());
              return Py_None;
          }
        } catch (std::bad_any_cast &) {
          // C++ datatypes can really be anything, this just captures
          // mislabelled data, it really shouldn't happen
          std::string message =
              std::string("Unhandled type conversion occured for property: ") +
              rdvalue.key;
          UNDER_CONSTRUCTION(message.c_str());
          return Py_None;
        }
      }
    }
  }
  PyErr_SetString(PyExc_KeyError, key.c_str());
  return nullptr;
}

// Return policy for functions that directly return a PyObject* and
// are fully responsible for setting the Python error state.
struct return_pyobject_passthrough {
  template <class T>
  struct apply {
    struct type {
      static bool convertible() { return true; }

      PyObject *operator()(PyObject *inner) const { return inner; }
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
      PyTypeObject const *get_pytype() const {
        return boost::python::converter::expected_pytype_for_arg<
            T>::get_pytype();
      }
#endif
    };
  };
};

template <class RDOb>
int MolHasProp(const RDOb &mol, const std::string &key) {
  int res = mol.hasProp(key);
  // std::cout << "key: "  << key << ": " << res << std::endl;
  return res;
}

template <class RDOb, class T>
void MolSetProp(const RDOb &mol, const std::string &key, const T &val,
                bool computed = false) {
  mol.setProp(key, val, computed);
}

template <class RDOb>
void MolClearProp(const RDOb &mol, const std::string &key) {
  if (!mol.hasProp(key)) {
    return;
  }
  mol.clearProp(key);
}

template <class RDOb>
void MolClearComputedProps(const RDOb &mol) {
  mol.clearComputedProps();
}
}  // namespace RDKit

#endif
