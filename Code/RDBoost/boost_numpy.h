// Boost python numpy available in Boost 1.63+
// Boost python numeric removed in Boost 1.65+
#if BOOST_VERSION < 106500
#include <boost/python/numeric.hpp>
typedef typename boost::python::numeric::array NumpyArrayType;
#else
#include <boost/python/numpy.hpp>
typedef typename boost::python::numpy::ndarray NumpyArrayType;
#endif
