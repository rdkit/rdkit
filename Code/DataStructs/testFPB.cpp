//
//  Copyright (C) 2015 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/FPBReader.h>

using namespace RDKit;

void test1FPBReaderBasics() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing FPBReader basics "
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {  // operator[] version
      std::pair<ExplicitBitVect *, std::string> tpl = fps[0];
      ExplicitBitVect *fp = tpl.first;
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      delete fp;
      TEST_ASSERT(tpl.second == "ZINC00902219");
    }
    {  // get* version
      ExplicitBitVect *fp = fps.getFP(0);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      delete fp;
      std::string nm = fps.getId(0);
      TEST_ASSERT(nm == "ZINC00902219");
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();

  test1FPBReaderBasics();

  return 0;
}
