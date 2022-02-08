//
//  Copyright (C) 2016 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/utils.h>

#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/FPBReader.h>

using namespace RDKit;

void _basicsTest(FPBReader &fps) {
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
    // std::cerr << " nm: >" << nm << "<" << std::endl;
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

void test1FPBReaderBasics() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing FPBReader basics "
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    _basicsTest(fps);
    fps.cleanup();  // make sure this doesn't cause problems
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test9LazyFPBReaderBasics() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Lazy FPBReader basics "
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    _basicsTest(fps);
    fps.cleanup();  // make sure this doesn't cause problems
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
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
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
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold, no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30, false);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
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
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTanimoto(0, bytes), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(0);
      TEST_ASSERT(ebv);
      TEST_ASSERT(feq(fps.getTanimoto(0, *ebv.get()), 1.0));
      TEST_ASSERT(feq(fps.getTanimoto(1, *ebv.get()), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
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
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 95);
      TEST_ASSERT(feq(nbrs[1].first, 0.4125));
      TEST_ASSERT(nbrs[1].second == 89);
    }
    {  // ebv with a threshold
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(95);
      TEST_ASSERT(ebv);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(*ebv.get(), 0.30);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 95);
      TEST_ASSERT(feq(nbrs[1].first, 0.4125));
      TEST_ASSERT(nbrs[1].second == 89);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

// need forward declarations of some of the detail functions (doesn't make sense
// to include these in a public interface and there aren't enough of them to
// make it worth adding a detail header)
namespace RDKit {
namespace detail {
boost::dynamic_bitset<> *bytesToBitset(const std::uint8_t *fpData,
                                       std::uint32_t nBits);
std::uint8_t *bitsetToBytes(const boost::dynamic_bitset<> &bitset);
}  // namespace detail
}  // namespace RDKit
void test7BitsetDetails() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing some internal bitset details"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {  // test round-tripping bytes -> dynamic_bitest -> bytes
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);

      // -------------------
      // start with bytes -> a bitset
      boost::dynamic_bitset<> *dbs =
          RDKit::detail::bytesToBitset(bytes.get(), fps.nBits());
      TEST_ASSERT(dbs);
      TEST_ASSERT(dbs->size() == fps.nBits());
      TEST_ASSERT(dbs->count() == 17);

      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < dbs->count(); ++i) {
        TEST_ASSERT((*dbs)[obs[i]]);
      }

      // -------------------
      // and now go the other way
      std::uint8_t *newBytes = RDKit::detail::bitsetToBytes(*dbs);
      TEST_ASSERT(newBytes);
      for (unsigned int i = 0; i < fps.nBits() / 8; ++i) {
        TEST_ASSERT(newBytes[i] == bytes[i]);
      }

      delete dbs;
      delete[] newBytes;
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test8FPBReaderContains() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing FPBReader contains search"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(nbrs[0] == 0);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      TEST_ASSERT(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 4);
      TEST_ASSERT(nbrs[0] == 1);
      TEST_ASSERT(nbrs[1] == 2);
      TEST_ASSERT(nbrs[2] == 3);
      TEST_ASSERT(nbrs[3] == 4);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(16);
      TEST_ASSERT(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(nbrs[0] == 16);
      TEST_ASSERT(nbrs[1] == 17);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(87);
      TEST_ASSERT(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 4);
      TEST_ASSERT(nbrs[0] == 85);
      TEST_ASSERT(nbrs[1] == 86);
      TEST_ASSERT(nbrs[2] == 87);
      TEST_ASSERT(nbrs[3] == 88);
    }
    {
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(87);
      TEST_ASSERT(ebv);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(*ebv.get());
      TEST_ASSERT(nbrs.size() == 4);
      TEST_ASSERT(nbrs[0] == 85);
      TEST_ASSERT(nbrs[1] == 86);
      TEST_ASSERT(nbrs[2] == 87);
      TEST_ASSERT(nbrs[3] == 88);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test9FPBReaderTversky() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing FPBReader Tversky "
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTversky(0, bytes, 1., 1.), 1.0));
      TEST_ASSERT(feq(fps.getTversky(1, bytes, 1., 1.), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      TEST_ASSERT(bytes);
      TEST_ASSERT(feq(fps.getTversky(1, bytes, 1., 1.), 1.0));
      TEST_ASSERT(feq(fps.getTversky(0, bytes, 1., 1.), 0.3703));
      TEST_ASSERT(feq(fps.getTversky(2, bytes, 1, 1.), 1.0));
      TEST_ASSERT(feq(fps.getTversky(5, bytes, 1., 1.), 0.2903));
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test10FPBReaderTverskyNeighbors() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing FPBReader Tversky neighbors"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    TEST_ASSERT(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1.);
      TEST_ASSERT(nbrs.size() == 1);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1., 0.3);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold, asymmetric
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 0.5, 0.3);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.4255));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold,  no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1., 0.3, false);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.3703));
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {  // with a threshold, asymmetric, no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 0.5, 0.3, false);
      TEST_ASSERT(nbrs.size() == 5);
      TEST_ASSERT(feq(nbrs[0].first, 1.));
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(feq(nbrs[1].first, 0.4255));
      TEST_ASSERT(nbrs[1].second == 1);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub1118() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github issue "
                          "1118: deleting non-initialized FPBReader causes seg "
                          "fault"
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    auto *fps = new FPBReader(filename);
    delete fps;
  }
  {
    auto *fps = new FPBReader();
    delete fps;
  }
  {
    std::string filename = pathName + "zim.head100.fpb";
    std::ifstream ifs(filename.c_str());
    auto *fps = new FPBReader(&ifs, false);
    delete fps;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();

  test1FPBReaderBasics();
  test9LazyFPBReaderBasics();

  test2FPBReaderTanimoto();
  test3FPBReaderTanimotoNeighbors();
  test8FPBReaderContains();
  test4LazyFPBReaderBasics();
  test5LazyFPBReaderTanimoto();
  test6LazyFPBReaderTanimotoNeighbors();
  test7BitsetDetails();

  test9FPBReaderTversky();
  test10FPBReaderTverskyNeighbors();
  testGithub1118();
  return 0;
}
