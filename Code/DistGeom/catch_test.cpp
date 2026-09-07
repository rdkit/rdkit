//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#ifdef RDK_TEST_MULTITHREADED
#include <csignal>
#include <thread>
#include <chrono>
#endif

// using namespace RDKit;