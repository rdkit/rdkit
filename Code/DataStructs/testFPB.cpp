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
    {  // get* version
      std::string nm = fps.getId(0);
      TEST_ASSERT(nm == "ZINC00902219");
      ExplicitBitVect *fp = fps.getFP(0);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      delete fp;
    }
    {  // operator[] version
      std::pair<ExplicitBitVect *, std::string> tpl = fps[0];
      ExplicitBitVect *fp = tpl.first;
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      delete fp;
      TEST_ASSERT(tpl.second == "ZINC00902219");
    }
    {  // test another fp
      ExplicitBitVect *fp = fps.getFP(3);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 20);
      unsigned int obs[20] = {1,   8,    80,   95,   222,  227, 457,
                              482, 650,  680,  715,  807,  831, 845,
                              888, 1226, 1556, 1711, 1917, 1982};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      delete fp;
      std::string nm = fps.getId(3);
      TEST_ASSERT(nm == "ZINC04803506");
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();

  test1FPBReaderBasics();

  return 0;
}
