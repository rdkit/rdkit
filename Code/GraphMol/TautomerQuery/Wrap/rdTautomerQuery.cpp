//
// Created by Gareth Jones on 5/30/2020.
//
// Copyright 2020 Schrodinger, Inc
//

#include "rdTautomerQuery.h"
#include <RDBoost/python.h>

namespace python = boost::python;

void wrap_TautomerQuery();

BOOST_PYTHON_MODULE(rdTautomerQuery) { wrap_TautomerQuery(); }
