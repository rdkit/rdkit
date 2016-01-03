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
#include <RDGeneral/utils.h>

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
    TEST_ASSERT(fps.nBits() == 2048);

    {  // pop counts
      std::pair<unsigned int, unsigned int> offsets;
      offsets = fps.getFPIdsInCountRange(17, 17);
      TEST_ASSERT(offsets.first == 0);
      TEST_ASSERT(offsets.second == 1);
      offsets = fps.getFPIdsInCountRange(60, 65);
      TEST_ASSERT(offsets.first == 96);
      TEST_ASSERT(offsets.second == 100);
      offsets = fps.getFPIdsInCountRange(160, 165);
      TEST_ASSERT(offsets.first == 100);
      TEST_ASSERT(offsets.second == 100);
    }
    {  // get* version
      std::string nm = fps.getId(0);
      TEST_ASSERT(nm == "ZINC00902219");
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(0);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
    }
    {  // operator[] version
      std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> tpl = fps[0];
      boost::shared_ptr<ExplicitBitVect> fp = tpl.first;
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      TEST_ASSERT(tpl.second == "ZINC00902219");
    }
    {  // test another fp
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(3);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 20);
      unsigned int obs[20] = {1,   8,    80,   95,   222,  227, 457,
                              482, 650,  680,  715,  807,  831, 845,
                              888, 1226, 1556, 1711, 1917, 1982};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      std::string nm = fps.getId(3);
      TEST_ASSERT(nm == "ZINC04803506");
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test2FPBReaderTanimoto() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing FPBReader tanimoto " << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(1);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 0.3703));
      TEST_ASSERT(feq(fps.getTanimoto(2, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(5, bytes), 0.2903));
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test3FPBReaderTanimotoNeighbors() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing FPBReader tanimoto neighbors"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 95);
      TEST_ASSERT(feq(nbrs[1].first, 0.4125));
      TEST_ASSERT(nbrs[1].second == 89);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test4LazyFPBReaderBasics() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Lazy FPBReader basics "
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    TEST_ASSERT(fps.nBits() == 2048);

    {  // pop counts
      std::pair<unsigned int, unsigned int> offsets;
      offsets = fps.getFPIdsInCountRange(17, 17);
      TEST_ASSERT(offsets.first == 0);
      TEST_ASSERT(offsets.second == 1);
      offsets = fps.getFPIdsInCountRange(60, 65);
      TEST_ASSERT(offsets.first == 96);
      TEST_ASSERT(offsets.second == 100);
      offsets = fps.getFPIdsInCountRange(160, 165);
      TEST_ASSERT(offsets.first == 100);
      TEST_ASSERT(offsets.second == 100);
    }
    {  // get* version
      std::string nm = fps.getId(0);
      TEST_ASSERT(nm == "ZINC00902219");
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(0);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
    }
    {  // operator[] version
      std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> tpl = fps[0];
      boost::shared_ptr<ExplicitBitVect> fp = tpl.first;
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      TEST_ASSERT(tpl.second == "ZINC00902219");
    }
    {  // test another fp
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(3);
      TEST_ASSERT(fp);
      TEST_ASSERT(fp->getNumBits() == 2048);
      TEST_ASSERT(fp->getNumOnBits() == 20);
      unsigned int obs[20] = {1,   8,    80,   95,   222,  227, 457,
                              482, 650,  680,  715,  807,  831, 845,
                              888, 1226, 1556, 1711, 1917, 1982};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        TEST_ASSERT(fp->getBit(obs[i]));
      }
      std::string nm = fps.getId(3);
      TEST_ASSERT(nm == "ZINC04803506");
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test5LazyFPBReaderTanimoto() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Lazy FPBReader tanimoto "
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(1);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 0.3703));
      TEST_ASSERT(feq(fps.getTanimoto(2, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(5, bytes), 0.2903));
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test6LazyFPBReaderTanimotoNeighbors() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Lazy FPBReader tanimoto neighbors"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes = fps.getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int> > nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 95);
      TEST_ASSERT(feq(nbrs[1].first, 0.4125));
      TEST_ASSERT(nbrs[1].second == 89);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();

  test1FPBReaderBasics();
  test2FPBReaderTanimoto();
  test3FPBReaderTanimotoNeighbors();
  test4LazyFPBReaderBasics();
  test5LazyFPBReaderTanimoto();
  test6LazyFPBReaderTanimotoNeighbors();

  return 0;
}
