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

#include <boost/serialization/array.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
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
bool serialize(const T& data, const std::ofstream& ofs);

template <typename T>
bool deserialize(T& data, const std::string& filename);

template <typename T>
bool deserialize(T& data, const std::ifstream& ifs);

template <typename T>
bool deserializeAll(std::vector<T>* weights, std::vector<T>* biases,
                    std::string& filename, std::string atomType);

template <typename T>
bool deserializeAll(std::vector<T>* weights, std::vector<T>* biases,
                    std::ifstream& ifs, std::string atomType);

template <typename T>
bool serializeAll(
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, T>>>>*
        weightsAndBiasesForEachAtomType,
    std::ofstream& ofs);

/*!
    Stores boost serialized eigen matrix in "fileName"
    \param weightsAndBiasesForEachAtomType formatted as follows
        H : "weight"    -> ArrayXXd/ArrayXXf
            "bias"      -> ArrayXXd/ArrayXXf
            "weight"    -> ArrayXXd/ArrayXXf
            "bias"      -> ArrayXXd/ArrayXXf
                            (in order of layers of NN)
        O : "weight"    -> ArrayXXd/ArrayXXf
            "bias"      -> ArrayXXd/ArrayXXf
            "weight"    -> ArrayXXd/ArrayXXf
            "bias"      -> ArrayXXd/ArrayXXf
        and so on for different atom types
    \param fileName     File in which the the first argument is stored
    \return true/false if values were stored or not
*/

template <typename T>
bool serializeAll(
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, T>>>>*
        weightsAndBiasesForEachAtomType,
    const std::string& fileName);
}  // namespace EigenSerializer
}  // namespace RDNumeric