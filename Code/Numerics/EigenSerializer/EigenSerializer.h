//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
#ifndef EIGEN_CONFIG_H_
#define EIGEN_CONFIG_H_

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>

#include "portable_binary_iarchive.hpp"
#include "portable_binary_oarchive.hpp"
#define EIGEN_DENSEBASE_PLUGIN "Numerics/EigenSerializer/EigenBaseAddons.h"
#include <Eigen/Dense>
#endif

#include <fstream>

namespace RDNumeric {
namespace EigenSerializer {

template <typename T>
bool serialize(const T& data, const std::string& filename);

template <typename T>
bool serialize(const T& data, std::ofstream& ofs);

template <typename T>
bool deserialize(T& data, const std::string& filename);

template <typename T>
bool deserialize(T& data, std::ifstream& ifs);

template <typename T>
bool deserializeAll(std::vector<T>& data, std::vector<std::string>& labels,
                    const std::string& filename);

template <typename T>
bool deserializeAll(std::vector<T>& data, std::vector<std::string>& labels,
                    std::ifstream& ifs);

template <typename T>
bool serializeAll(const std::vector<T>& data,
                  const std::vector<std::string>& labels,
                  const std::string& filename);

template <typename T>
bool serializeAll(const std::vector<T>& data,
                  const std::vector<std::string>& labels,
                  std::ofstream& ofs);
}  // namespace EigenSerializer
}  // namespace RDNumeric