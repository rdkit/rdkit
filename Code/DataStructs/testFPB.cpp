//
//  Copyright (C) 2016-2025 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <RDGeneral/utils.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/FPBReader.h>

using namespace RDKit;

void _basicsTest(FPBReader &fps) {
  fps.init();
  REQUIRE(fps.length() == 100);
  REQUIRE(fps.nBits() == 2048);

  {  // pop counts
    std::pair<unsigned int, unsigned int> offsets;
    offsets = fps.getFPIdsInCountRange(17, 17);
    REQUIRE(offsets.first == 0);
    REQUIRE(offsets.second == 1);
    offsets = fps.getFPIdsInCountRange(60, 65);
    REQUIRE(offsets.first == 96);
    REQUIRE(offsets.second == 100);
    offsets = fps.getFPIdsInCountRange(160, 165);
    REQUIRE(offsets.first == 100);
    REQUIRE(offsets.second == 100);
  }
  {  // get* version
    std::string nm = fps.getId(0);
    REQUIRE(nm == "ZINC00902219");
    boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(0);
    REQUIRE(fp);
    REQUIRE(fp->getNumBits() == 2048);
    REQUIRE(fp->getNumOnBits() == 17);
    unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                            811, 831, 888, 1335, 1411, 1664, 1820, 1917};
    for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
      REQUIRE(fp->getBit(obs[i]));
    }
  }
  {  // operator[] version
    std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> tpl = fps[0];
    boost::shared_ptr<ExplicitBitVect> fp = tpl.first;
    REQUIRE(fp);
    REQUIRE(fp->getNumBits() == 2048);
    REQUIRE(fp->getNumOnBits() == 17);
    unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                            811, 831, 888, 1335, 1411, 1664, 1820, 1917};
    for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
      REQUIRE(fp->getBit(obs[i]));
    }
    REQUIRE(tpl.second == "ZINC00902219");
  }
  {  // test another fp
    boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(3);
    REQUIRE(fp);
    REQUIRE(fp->getNumBits() == 2048);
    REQUIRE(fp->getNumOnBits() == 20);
    unsigned int obs[20] = {1,   8,    80,   95,   222,  227, 457,
                            482, 650,  680,  715,  807,  831, 845,
                            888, 1226, 1556, 1711, 1917, 1982};
    for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
      REQUIRE(fp->getBit(obs[i]));
    }
    std::string nm = fps.getId(3);
    REQUIRE(nm == "ZINC04803506");
  }
}

TEST_CASE("FPBReader Basics") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    _basicsTest(fps);
    fps.cleanup();  // make sure this doesn't cause problems
  }
}

TEST_CASE("Lazy FPBReader Basics") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    _basicsTest(fps);
    fps.cleanup();  // make sure this doesn't cause problems
  }
}

TEST_CASE("FPBReader Tanimoto") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTanimoto(0, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTanimoto(1, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(0, bytes), 0.3703));
      REQUIRE(feq(fps.getTanimoto(2, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(5, bytes), 0.2903));
    }
  }
}

TEST_CASE("FPBReader Tanimoto Neighbors") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes);
      REQUIRE(nbrs.size() == 1);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.3703));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold, no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30, false);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.3703));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(95);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 95);
      REQUIRE(feq(nbrs[1].first, 0.4125));
      REQUIRE(nbrs[1].second == 89);
    }
  }
}

TEST_CASE("Lazy FPBReader Basics 2") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    REQUIRE(fps.length() == 100);
    REQUIRE(fps.nBits() == 2048);

    {  // pop counts
      std::pair<unsigned int, unsigned int> offsets;
      offsets = fps.getFPIdsInCountRange(17, 17);
      REQUIRE(offsets.first == 0);
      REQUIRE(offsets.second == 1);
      offsets = fps.getFPIdsInCountRange(60, 65);
      REQUIRE(offsets.first == 96);
      REQUIRE(offsets.second == 100);
      offsets = fps.getFPIdsInCountRange(160, 165);
      REQUIRE(offsets.first == 100);
      REQUIRE(offsets.second == 100);
    }
    {  // get* version
      std::string nm = fps.getId(0);
      REQUIRE(nm == "ZINC00902219");
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(0);
      REQUIRE(fp);
      REQUIRE(fp->getNumBits() == 2048);
      REQUIRE(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        REQUIRE(fp->getBit(obs[i]));
      }
    }
    {  // operator[] version
      std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> tpl = fps[0];
      boost::shared_ptr<ExplicitBitVect> fp = tpl.first;
      REQUIRE(fp);
      REQUIRE(fp->getNumBits() == 2048);
      REQUIRE(fp->getNumOnBits() == 17);
      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        REQUIRE(fp->getBit(obs[i]));
      }
      REQUIRE(tpl.second == "ZINC00902219");
    }
    {  // test another fp
      boost::shared_ptr<ExplicitBitVect> fp = fps.getFP(3);
      REQUIRE(fp);
      REQUIRE(fp->getNumBits() == 2048);
      REQUIRE(fp->getNumOnBits() == 20);
      unsigned int obs[20] = {1,   8,    80,   95,   222,  227, 457,
                              482, 650,  680,  715,  807,  831, 845,
                              888, 1226, 1556, 1711, 1917, 1982};
      for (unsigned int i = 0; i < fp->getNumOnBits(); ++i) {
        REQUIRE(fp->getBit(obs[i]));
      }
      std::string nm = fps.getId(3);
      REQUIRE(nm == "ZINC04803506");
    }
  }
}

TEST_CASE("Lazy FPBReader Tanimoto") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTanimoto(0, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(1, bytes), 0.3703));
    }
    {
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(0);
      REQUIRE(ebv);
      REQUIRE(feq(fps.getTanimoto(0, *ebv.get()), 1.0));
      REQUIRE(feq(fps.getTanimoto(1, *ebv.get()), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTanimoto(1, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(0, bytes), 0.3703));
      REQUIRE(feq(fps.getTanimoto(2, bytes), 1.0));
      REQUIRE(feq(fps.getTanimoto(5, bytes), 0.2903));
    }
  }
}

TEST_CASE("Lazy FPBReader Tanimoto Neighbors") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes);
      REQUIRE(nbrs.size() == 1);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.3703));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(95);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 95);
      REQUIRE(feq(nbrs[1].first, 0.4125));
      REQUIRE(nbrs[1].second == 89);
    }
    {  // ebv with a threshold
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(95);
      REQUIRE(ebv);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTanimotoNeighbors(*ebv.get(), 0.30);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 95);
      REQUIRE(feq(nbrs[1].first, 0.4125));
      REQUIRE(nbrs[1].second == 89);
    }
  }
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
TEST_CASE("Bitset Details") {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing some internal bitset details"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {  // test round-tripping bytes -> dynamic_bitest -> bytes
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename, true);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);

      // -------------------
      // start with bytes -> a bitset
      boost::dynamic_bitset<> *dbs =
          RDKit::detail::bytesToBitset(bytes.get(), fps.nBits());
      REQUIRE(dbs);
      REQUIRE(dbs->size() == fps.nBits());
      REQUIRE(dbs->count() == 17);

      unsigned int obs[17] = {1,   80,  183, 222,  227,  231,  482,  650, 807,
                              811, 831, 888, 1335, 1411, 1664, 1820, 1917};
      for (unsigned int i = 0; i < dbs->count(); ++i) {
        REQUIRE((*dbs)[obs[i]]);
      }

      // -------------------
      // and now go the other way
      std::uint8_t *newBytes = RDKit::detail::bitsetToBytes(*dbs);
      REQUIRE(newBytes);
      for (unsigned int i = 0; i < fps.nBits() / 8; ++i) {
        REQUIRE(newBytes[i] == bytes[i]);
      }

      delete dbs;
      delete[] newBytes;
    }
  }
}

TEST_CASE("FPBReader Contains") {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing FPBReader contains search"
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 1);
      REQUIRE(nbrs[0] == 0);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      REQUIRE(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 4);
      REQUIRE(nbrs[0] == 1);
      REQUIRE(nbrs[1] == 2);
      REQUIRE(nbrs[2] == 3);
      REQUIRE(nbrs[3] == 4);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(16);
      REQUIRE(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(nbrs[0] == 16);
      REQUIRE(nbrs[1] == 17);
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(87);
      REQUIRE(bytes);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 4);
      REQUIRE(nbrs[0] == 85);
      REQUIRE(nbrs[1] == 86);
      REQUIRE(nbrs[2] == 87);
      REQUIRE(nbrs[3] == 88);
    }
    {
      boost::shared_ptr<ExplicitBitVect> ebv = fps.getFP(87);
      REQUIRE(ebv);
      std::vector<unsigned int> nbrs = fps.getContainingNeighbors(*ebv.get());
      REQUIRE(nbrs.size() == 4);
      REQUIRE(nbrs[0] == 85);
      REQUIRE(nbrs[1] == 86);
      REQUIRE(nbrs[2] == 87);
      REQUIRE(nbrs[3] == 88);
    }
  }
}

TEST_CASE("FPBReader Tversky") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTversky(0, bytes, 1., 1.), 1.0));
      REQUIRE(feq(fps.getTversky(1, bytes, 1., 1.), 0.3703));
    }
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(1);
      REQUIRE(bytes);
      REQUIRE(feq(fps.getTversky(1, bytes, 1., 1.), 1.0));
      REQUIRE(feq(fps.getTversky(0, bytes, 1., 1.), 0.3703));
      REQUIRE(feq(fps.getTversky(2, bytes, 1, 1.), 1.0));
      REQUIRE(feq(fps.getTversky(5, bytes, 1., 1.), 0.2903));
    }
  }
}

TEST_CASE("FPBReader Tversky Neighbors") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps(filename);
    fps.init();
    REQUIRE(fps.length() == 100);
    {
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1.);
      REQUIRE(nbrs.size() == 1);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1., 0.3);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.3703));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold, asymmetric
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 0.5, 0.3);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.4255));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold,  no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 1., 0.3, false);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.3703));
      REQUIRE(nbrs[1].second == 1);
    }
    {  // with a threshold, asymmetric, no screen
      boost::shared_array<std::uint8_t> bytes = fps.getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<double, unsigned int>> nbrs =
          fps.getTverskyNeighbors(bytes, 1., 0.5, 0.3, false);
      REQUIRE(nbrs.size() == 5);
      REQUIRE(feq(nbrs[0].first, 1.));
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(feq(nbrs[1].first, 0.4255));
      REQUIRE(nbrs[1].second == 1);
    }
  }
}

TEST_CASE("Github Issue 1118") {
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

int main(int argc, char *argv[]) { return Catch::Session().run(argc, argv); }
