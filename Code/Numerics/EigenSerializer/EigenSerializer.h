#ifndef EIGEN_CONFIG_H_
#define EIGEN_CONFIG_H_

#include <boost/serialization/array.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#define EIGEN_DENSEBASE_PLUGIN "Numerics/EigenSerializer/EigenBaseAddons.h"
#include <Eigen/Dense>
#endif

#include <fstream>

namespace RDNumeric {
namespace EigenSerializer {

template <typename T>
bool serialize(const T& data, const std::string& filename);

template <typename T>
bool deSerialize(T& data, const std::string& filename);

}  // namespace EigenSerializer
}  // namespace RDNumeric