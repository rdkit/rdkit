#include "rdDeprotect.h"
#include <RDBoost/python.h>

namespace python = boost::python;

void wrap_deprotect();

BOOST_PYTHON_MODULE(rdDeprotect) { wrap_deprotect(); }
