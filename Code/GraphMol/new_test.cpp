//
//  Copyright (C) 2001-2018 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>

using namespace std;
using namespace RDKit;

// -------------------------------------------------------------------
void testGetProp() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Get Property" << std::endl;
  RWMol m2;
  bool caughtKeyError = false;

  try {
    std::string val = m2.getProp<std::string>("prop1");
  }
  catch (const KeyErrorException& e) {
    caughtKeyError = true;
  }
  TEST_ASSERT(caughtKeyError);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

struct GetPropDefaultTestCase {
  std::string propName;
  std::string defaultValue;
};

void testGetPropDefault() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Get Property with Default" << std::endl;
  RWMol m2;

  std::vector<GetPropDefaultTestCase> testCases = {
      {"prop1", "default_value"},
      {"prop1", "42"},
      {"prop1", "3.14"}
  };
  for (const auto& testCase : testCases) {
    std::string val = m2.getProp<std::string>(testCase.propName, testCase.defaultValue);
    TEST_ASSERT(val == testCase.defaultValue);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();
  boost::logging::enable_logs("rdApp.info");
  testGetProp();
  testGetPropDefault();
  return 0;
}
