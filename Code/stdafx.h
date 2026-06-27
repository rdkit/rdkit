//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Aggregate header used as the precompiled header when
// RDK_USE_PRECOMPILED_HEADERS is enabled (see Code/CMakeLists.txt). It bundles
// the STL and boost headers that dominate parse time across the tree, plus the
// core RDKit headers pulled in by almost every translation unit. It is not
// included directly by any source file.
//
#pragma once

#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
