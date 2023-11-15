//
//  Copyright (C) 2021 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <RDGeneral/RDLog.h>

int main(int argc, char* argv[]) {
  RDLog::InitLogs();

  int result = Catch::Session().run(argc, argv);

  return result;
}