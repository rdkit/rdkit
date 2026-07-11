// Boost python numpy available in Boost 1.63+
// Boost python numeric removed in Boost 1.65+
#include <RDGeneral/export.h>
#if BOOST_VERSION < 106500
#include <boost/python/numeric.hpp>
using NumpyArrayType = boost::python::numeric::array;
#else
#include <boost/python/numpy.hpp>
using NumpyArrayType = boost::python::numpy::ndarray;
#endif
