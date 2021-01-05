//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.


#ifndef UTIL_H_
#define UTIL_H_

#include <cmath>
#include <limits>
#include <cstdio>
#include <ctime>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <boost/optional.hpp>
#include "export.h"

/*
 * Utility methods
 */
namespace GarethUtil {

    using namespace std;

/*
 * Return current time as string
 */
    string currentTime();

/**
 *
 * @param str
 * @param prefix
 * @return true id str starts with prefix
 */
    GA_EXPORT bool startsWith(string str, string prefix);

/**
 *
 * @return current user name
 */
    GA_EXPORT string getUserName();

/**
 * Template function to convert a string to a type.  Can optionally process a substring
 *
 * @param str
 * @param ok optional boolean reference to indicate success or failure
 * @param start optional start position
 * @param end optional end position
 * @return
 */
    template<typename T>
    T convertString(const string &str, bool *const ok =
    nullptr, size_t start = 0, size_t noChars = 0) {
        if (ok != nullptr)
            *ok = false;
        T rtn;
        if (str.empty())
            return rtn;
        if (start >= str.length())
            return rtn;
        if (noChars == 0)
            noChars = str.length() - start;
        string field = str.substr(start, noChars);
        stringstream ss(field);
        ss >> rtn;
        if (!ss.fail()) {
            if (ok != nullptr)
                *ok = true;
        }
        return rtn;
    }

/**
 * Template function to write all elements of a vector to a string
 *
 * @param vec
 * @param seperator
 * @return
 */
    template<typename T>
    string collectionToString(const T &vec,
                              const string &seperator = ", ") {
        stringstream ss;
        std::copy(vec.begin(), vec.end(),
                  ostream_iterator<typename T::value_type>(ss, seperator.c_str()));
        // erase trailing separator
        string rtn = ss.str();
        if (!rtn.empty())
            rtn.erase(rtn.length() - seperator.length(), seperator.length());
        return rtn;
    }

/**
 * Removes any trailing linefeed (\r) from a string
 * @param line
 * @return
 */
    GA_EXPORT string &removeTrailingLF(string &line);

/**
 * Determines if a key is present in a map
 *
 * @param key
 * @param map
 * @param value  optionally set the value for the key
 * @return
 */
    template<typename K, typename V>
    bool isKeyPresent(const K key, const map<K, V> &map,
                      V *const value = nullptr) {
        auto iter = map.find(key);
        if (iter != map.end()) {
            if (value != nullptr)
                *value = iter->second;
            return true;
        }
        return false;
    }

/**
 * Bs whitespace from beginning and end of string.  Original argument is modified and returned.
 * @param str
 * @return
 */
    GA_EXPORT string &trim(string &str);

/**
 * Converts string to upper case.  Original argument is modified and returned.
 * @param str
 * @return
 */
    GA_EXPORT string &toUpperCase(string &str);

/**
 * Converts string to lower case.  Original argument is modified and returned.
 * @param str
 * @return
 */
    GA_EXPORT string &toLowerCase(string &str);

/**
 * @param str1
 * @param str2
 * @return true if the two strings are equal
 */
    GA_EXPORT bool equals(const string &str1, const string &str2);

/**
 * @param str1
 * @param str2
 * @return true if the two strings are equal (case insensitive)
 */
    GA_EXPORT bool equalsIgnoreCase(const string &str1, const string &str2);

/**
 *
 * @param str
 * @param suffix
 * @return true if str ends with suffix
 */
    GA_EXPORT bool endsWith(const string &str, const string &suffix);

/**
 *
 * @param vector
 * @param value
 * @return true if vector contains value
 */
    template<typename T>
    bool contains(std::vector<T> vector, T value) {
        return std::find(vector.begin(), vector.end(), value) != vector.end();
    }

/**
 * Determines if two double numbers are equal within multiples of machine precision
 *
 * @param d1
 * @param d2
 * @param ulp machine precision (units in the last place)
 * @return
 */
    GA_EXPORT bool equals(const double d1, const double d2, const int ulp);

/**
 * Returns true if the two numbers are within epsilon of each other
 *
 * @param d1
 * @param d2
 * @param epsilon
 * @return
 */
    GA_EXPORT bool equals(const double d1, const double d2, const double epsilon);

/**
 * Determines if two double numbers are within one unit of the last place.
 *
 * @param d1
 * @param d2
 * @return
 */
    GA_EXPORT bool equals(const double d1, const double d2);

/**
 * Finds the first occurrence of an matching item in a list
 *
 * @param values
 * @param matcher
 * @return
 */
    template<typename T>
    boost::optional<T> findFirstInList(const std::vector<T> &values,
                                       const function<bool(T)> &matcher) {
        auto iter = std::find_if(values.cbegin(), values.cend(), matcher);
        if (iter != values.cend()) {
            return *iter;
        }
        return boost::none;
    }

/**
 * Finds the first occurrence of an matching item in a list of unique pointers
 *
 * @param values
 * @param matcher
 * @return
 */
    template<typename T>
    boost::optional<T *> findFirstInUniquePtrList(
            const std::vector<unique_ptr<T>> &values,
            const function<bool(const unique_ptr<T> &)> &matcher) {
        auto iter = std::find_if(values.cbegin(), values.cend(), matcher);
        if (iter != values.cend()) {
            return iter->get();
        }
        return boost::none;
    }

/**
 * Removes all items from a list that DON'T match the filter.  Alters and returns the original list.
 *
 * @param values
 * @param filter
 * @return
 */
    template<typename T>
    std::vector<T> &filterList(std::vector<T> &values,
                               const function<bool(T &)> &filter) {
        auto iter = remove_if(values.begin(), values.end(), not1(filter));
        values.erase(iter, values.end());
        return values;
    }

/**
 * Creates a new vector of items in the input vector that match the filter
 * @param values
 * @param filter
 * @return
 */
    template<typename T>
    std::vector<T> filterListToNewList(const std::vector<T> &values,
                                       const function<bool(const T &)> &filter) {
        std::vector<T> newValues;
        newValues.reserve(values.size());
        copy_if(values.begin(), values.end(), back_inserter(newValues), filter);
        return newValues;
    }

/**
 * Creates a new vector of items in the input vector that match the filter.
 * This method is specialized for an input list of unique pointers.
 * A list of raw pointers is returned.  It is assumed that the input list will
 * outlast the filtered list.
 *
 * @param values
 * @param filter
 * @return
 */
    template<typename T>
    std::vector<T *> filterUniquePtrListToRawPtrList(
            const std::vector<unique_ptr<T>> &values,
            const function<bool(T &)> &filter) {
        std::vector<T *> newValues;
        newValues.reserve(values.size());
        for (const auto &value : values) {
            if (filter(*value)) {
                newValues.push_back(value.get());
            }
        }
        return newValues;
    }

/**
 * Specialization of the above for filtering non-const unique_ptrs to const raw pointer filtered list
 * @param values
 * @param filter
 * @return
 */
    template<typename T>
    std::vector<const T *> filterUniquePtrListToConstRawPtrList(
            const std::vector<unique_ptr<T>> &values,
            const function<bool(T &)> &filter) {
        std::vector<const T *> newValues;
        newValues.reserve(values.size());
        for (const auto &value : values) {
            if (filter(*value)) {
                newValues.push_back(value.get());
            }
        }
        return newValues;
    }

/**
 * Creates a new list from an existing list by applying the mapping function.
 *
 * @param values
 * @param map
 * @return
 */
    template<typename T, typename V>
    std::vector<V> mapToNewList(const std::vector<T> &values,
                                const function<V(const T &)> &map) {
        std::vector<V> mappedValues;
        mappedValues.reserve(values.size());
        std::transform(values.cbegin(), values.cend(), back_inserter(mappedValues),
                       map);
        return mappedValues;
    }

/**
 * Creates a new list from an existing list by applying the mapping function.
 *
 * @param values
 * @param map
 * @return
 */
    template<typename T>
    std::vector<T> mapToNewList(const std::vector<T> &values,
                                const function<T(const T &)> &map) {
        return mapToNewList<T, T>(values, map);
    }

/**
 * Applies a function to every value in a list.
 *
 * @param values
 * @param func
 */
    template<typename T>
    void forEach(const std::vector<T> &values, const function<void(T &)> &func) {
        for (const auto &v : values) {
            func(v);
        }
    }

/**
 * Reduces a vector to a single value of a different type using the accumulator function.
 *
 * @param startingValue
 * @param accumulator function takes accumulating value as first argument and current list value as second.
 * @return
 */
    template<typename T, typename R>
    R reduce(const std::vector<T> &values, const R &startingValue,
             const function<R(R &, const T &)> accumulator) {
        R currentValue = startingValue;
        for (auto iter = values.cbegin(); iter != values.cend(); ++iter) {
            currentValue = accumulator(currentValue, *iter);
        }
        return currentValue;
    }

/**
 * Writes an object to a string using the stream << operator.
 *
 * @param obj
 * @return
 */
    template<typename T>
    string printToString(T obj) {
        stringstream ss;
        ss << obj;
        return ss.str();
    }

/**
 * Retrieve an environment variable value, or none if it is not set.
 *
 * @param name
 * @return
 */
    GA_EXPORT boost::optional<string> getEnv(const string &name);

}

#endif /* UTIL_H_ */
