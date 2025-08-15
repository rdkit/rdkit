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
#include <DataStructs/MultiFPBReader.h>
#include <DataStructs/BitOps.h>

using namespace RDKit;

TEST_CASE("MultiFPBReader Basics") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 2);
    REQUIRE(mfps.getReader(0));
    REQUIRE(mfps.getReader(0)->nBits() == mfps.nBits());
  }
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    MultiFPBReader mfps;
    REQUIRE(mfps.addReader(&fps1) == 1);
    REQUIRE(mfps.addReader(&fps2) == 2);
    mfps.init();
    REQUIRE(mfps.length() == 2);
    REQUIRE(mfps.getReader(0));
    REQUIRE(mfps.getReader(0)->nBits() == mfps.nBits());
    REQUIRE(mfps.getReader(0));
    REQUIRE(mfps.getReader(0)->nBits() == mfps.nBits());
  }
}

TEST_CASE("MultiFPBReader Tanimoto") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 2);

    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 0);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 10);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 0);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.3703));
      REQUIRE(std::get<1>(nbrs[2]) == 1);
      REQUIRE(std::get<2>(nbrs[2]) == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(95);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      REQUIRE(nbrs.size() == 4);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 95);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 95);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.4125));
      REQUIRE(std::get<1>(nbrs[2]) == 89);
      REQUIRE(std::get<2>(nbrs[2]) == 0);
    }
  }
}

TEST_CASE("MultiFPBReader Tversky") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 2);

    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1.);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 0);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1., 0.30);
      REQUIRE(nbrs.size() == 10);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 0);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.3703));
      REQUIRE(std::get<1>(nbrs[2]) == 1);
      REQUIRE(std::get<2>(nbrs[2]) == 0);
    }
    {  // with a threshold, asymmetric
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 0.5, 0.30);
      REQUIRE(nbrs.size() == 10);
      REQUIRE(feq(std::get<0>(nbrs[0]), 1.));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[1]), 1.));
      REQUIRE(std::get<1>(nbrs[1]) == 0);
      REQUIRE(std::get<2>(nbrs[1]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.4255));
      REQUIRE(std::get<1>(nbrs[2]) == 1);
      REQUIRE(std::get<2>(nbrs[2]) == 0);
    }
  }
}

TEST_CASE("MultiFPBReader Contains") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 2);
    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      REQUIRE(bytes);
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 2);
      REQUIRE(nbrs[0].first == 0);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 0);
      REQUIRE(nbrs[1].second == 1);
    }
    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(1);
      REQUIRE(bytes);
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(bytes);
      REQUIRE(nbrs.size() == 8);
      REQUIRE(nbrs[0].first == 1);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[2].first == 2);
      REQUIRE(nbrs[4].first == 3);
      REQUIRE(nbrs[6].first == 4);
      REQUIRE(nbrs[1].first == 1);
      REQUIRE(nbrs[1].second == 1);
      REQUIRE(nbrs[3].first == 2);
      REQUIRE(nbrs[5].first == 3);
      REQUIRE(nbrs[7].first == 4);
    }
  }
}

TEST_CASE("MultiFPBReader Similarity Threaded") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    FPBReader fps1(pathName + "zinc_random200.1.patt.fpb");
    FPBReader fps2(pathName + "zinc_random200.2.patt.fpb");
    FPBReader fps3(pathName + "zinc_random200.3.patt.fpb");
    FPBReader fps4(pathName + "zinc_random200.4.patt.fpb");

    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    rdrs.push_back(&fps3);
    rdrs.push_back(&fps4);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 4);
    std::string fps =
        "0000000000404000100000001000040000300040222000002004000240000020000000"
        "8200010200000090000024040860070044003214820000220401054008018000226000"
        "4800800140000042000080008008020482400000200410800000300430200800400000"
        "0000080a0000800400010c800200648818100010880040";
    ExplicitBitVect qbv(1024);
    UpdateBitVectFromFPSText(qbv, fps);

    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6);
      REQUIRE(nbrs.size() == 6);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << std::get<0>(nbrs[i]) << " " <<
      //   std::get<1>(nbrs[i])
      //             << " " << std::get<2>(nbrs[i]) << " " << std::endl;
      // }
      REQUIRE(feq(std::get<0>(nbrs[0]), 0.66412));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 3);
      REQUIRE(feq(std::get<0>(nbrs[1]), 0.65289));
      REQUIRE(std::get<1>(nbrs[1]) == 1);
      REQUIRE(std::get<2>(nbrs[1]) == 2);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.64341));
      REQUIRE(std::get<1>(nbrs[2]) == 2);
      REQUIRE(std::get<2>(nbrs[2]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[3]), 0.61940));
      REQUIRE(std::get<1>(nbrs[3]) == 1);
      REQUIRE(std::get<2>(nbrs[3]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[4]), 0.61905));
      REQUIRE(std::get<1>(nbrs[4]) == 0);
      REQUIRE(std::get<2>(nbrs[4]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[5]), 0.61344));
      REQUIRE(std::get<1>(nbrs[5]) == 0);
      REQUIRE(std::get<2>(nbrs[5]) == 1);
    }

#ifdef RDK_TEST_MULTITHREADED
    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 4);
      REQUIRE(nbrs.size() == 6);
      REQUIRE(feq(std::get<0>(nbrs[0]), 0.66412));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 3);
      REQUIRE(feq(std::get<0>(nbrs[1]), 0.65289));
      REQUIRE(std::get<1>(nbrs[1]) == 1);
      REQUIRE(std::get<2>(nbrs[1]) == 2);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.64341));
      REQUIRE(std::get<1>(nbrs[2]) == 2);
      REQUIRE(std::get<2>(nbrs[2]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[3]), 0.61940));
      REQUIRE(std::get<1>(nbrs[3]) == 1);
      REQUIRE(std::get<2>(nbrs[3]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[4]), 0.61905));
      REQUIRE(std::get<1>(nbrs[4]) == 0);
      REQUIRE(std::get<2>(nbrs[4]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[5]), 0.61344));
      REQUIRE(std::get<1>(nbrs[5]) == 0);
      REQUIRE(std::get<2>(nbrs[5]) == 1);
    }
    {  // request more threads than we have readers, this shouldn't be a problem
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 8);
      REQUIRE(nbrs.size() == 6);
      REQUIRE(feq(std::get<0>(nbrs[0]), 0.66412));
      REQUIRE(std::get<1>(nbrs[0]) == 0);
      REQUIRE(std::get<2>(nbrs[0]) == 3);
      REQUIRE(feq(std::get<0>(nbrs[1]), 0.65289));
      REQUIRE(std::get<1>(nbrs[1]) == 1);
      REQUIRE(std::get<2>(nbrs[1]) == 2);
      REQUIRE(feq(std::get<0>(nbrs[2]), 0.64341));
      REQUIRE(std::get<1>(nbrs[2]) == 2);
      REQUIRE(std::get<2>(nbrs[2]) == 1);
      REQUIRE(feq(std::get<0>(nbrs[3]), 0.61940));
      REQUIRE(std::get<1>(nbrs[3]) == 1);
      REQUIRE(std::get<2>(nbrs[3]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[4]), 0.61905));
      REQUIRE(std::get<1>(nbrs[4]) == 0);
      REQUIRE(std::get<2>(nbrs[4]) == 0);
      REQUIRE(feq(std::get<0>(nbrs[5]), 0.61344));
      REQUIRE(std::get<1>(nbrs[5]) == 0);
      REQUIRE(std::get<2>(nbrs[5]) == 1);
    }

#endif
  }
}

TEST_CASE("MultiFPBReader Contains Threaded") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    FPBReader fps1(pathName + "zinc_random200.1.patt.fpb");
    FPBReader fps2(pathName + "zinc_random200.2.patt.fpb");
    FPBReader fps3(pathName + "zinc_random200.3.patt.fpb");
    FPBReader fps4(pathName + "zinc_random200.4.patt.fpb");

    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    rdrs.push_back(&fps3);
    rdrs.push_back(&fps4);
    MultiFPBReader mfps(rdrs);
    mfps.init();
    REQUIRE(mfps.length() == 4);
    std::string fps =
        "40081010824820021000500010110410003000402b20285000a4040240010030050000"
        "080001420040009000003d04086007080c03b31d920004220400074008098010206080"
        "00488001080000c64002a00080000200024c2000602410049200340820200002400010"
        "02200106090401056801080182006088101000088a0048";
    ExplicitBitVect qbv(1024);
    UpdateBitVectFromFPSText(qbv, fps);

    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv);
      REQUIRE(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      REQUIRE(nbrs[0].first == 160);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 163);
      REQUIRE(nbrs[1].second == 0);
      REQUIRE(nbrs[2].first == 170);
      REQUIRE(nbrs[2].second == 0);
      REQUIRE(nbrs[3].first == 180);
      REQUIRE(nbrs[3].second == 2);
      REQUIRE(nbrs[4].first == 182);
      REQUIRE(nbrs[4].second == 3);
      REQUIRE(nbrs[5].first == 185);
      REQUIRE(nbrs[5].second == 0);
      REQUIRE(nbrs[6].first == 189);
      REQUIRE(nbrs[6].second == 0);
      REQUIRE(nbrs[7].first == 192);
      REQUIRE(nbrs[7].second == 3);
      REQUIRE(nbrs[8].first == 193);
      REQUIRE(nbrs[8].second == 0);
    }
#ifdef RDK_TEST_MULTITHREADED
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 4);
      REQUIRE(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      REQUIRE(nbrs[0].first == 160);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 163);
      REQUIRE(nbrs[1].second == 0);
      REQUIRE(nbrs[2].first == 170);
      REQUIRE(nbrs[2].second == 0);
      REQUIRE(nbrs[3].first == 180);
      REQUIRE(nbrs[3].second == 2);
      REQUIRE(nbrs[4].first == 182);
      REQUIRE(nbrs[4].second == 3);
      REQUIRE(nbrs[5].first == 185);
      REQUIRE(nbrs[5].second == 0);
      REQUIRE(nbrs[6].first == 189);
      REQUIRE(nbrs[6].second == 0);
      REQUIRE(nbrs[7].first == 192);
      REQUIRE(nbrs[7].second == 3);
      REQUIRE(nbrs[8].first == 193);
      REQUIRE(nbrs[8].second == 0);
    }
    {  // request more threads than we have readers, this shouldn't be a problem
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 8);
      REQUIRE(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      REQUIRE(nbrs[0].first == 160);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 163);
      REQUIRE(nbrs[1].second == 0);
      REQUIRE(nbrs[2].first == 170);
      REQUIRE(nbrs[2].second == 0);
      REQUIRE(nbrs[3].first == 180);
      REQUIRE(nbrs[3].second == 2);
      REQUIRE(nbrs[4].first == 182);
      REQUIRE(nbrs[4].second == 3);
      REQUIRE(nbrs[5].first == 185);
      REQUIRE(nbrs[5].second == 0);
      REQUIRE(nbrs[6].first == 189);
      REQUIRE(nbrs[6].second == 0);
      REQUIRE(nbrs[7].first == 192);
      REQUIRE(nbrs[7].second == 3);
      REQUIRE(nbrs[8].first == 193);
      REQUIRE(nbrs[8].second == 0);
    }
#endif
  }

  {  // test with initOnSearch set
    FPBReader fps1(pathName + "zinc_random200.1.patt.fpb");
    FPBReader fps2(pathName + "zinc_random200.2.patt.fpb");
    FPBReader fps3(pathName + "zinc_random200.3.patt.fpb");
    FPBReader fps4(pathName + "zinc_random200.4.patt.fpb");

    std::vector<FPBReader *> rdrs;
    rdrs.push_back(&fps1);
    rdrs.push_back(&fps2);
    rdrs.push_back(&fps3);
    rdrs.push_back(&fps4);
    MultiFPBReader mfps(rdrs, false, true);
    std::string fps =
        "40081010824820021000500010110410003000402b20285000a4040240010030050000"
        "080001420040009000003d04086007080c03b31d920004220400074008098010206080"
        "00488001080000c64002a00080000200024c2000602410049200340820200002400010"
        "02200106090401056801080182006088101000088a0048";
    ExplicitBitVect qbv(1024);
    UpdateBitVectFromFPSText(qbv, fps);
#ifndef RDK_TEST_MULTITHREADED
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv);
      REQUIRE(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      REQUIRE(nbrs[0].first == 160);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 163);
      REQUIRE(nbrs[1].second == 0);
      REQUIRE(nbrs[2].first == 170);
      REQUIRE(nbrs[2].second == 0);
      REQUIRE(nbrs[3].first == 180);
      REQUIRE(nbrs[3].second == 2);
      REQUIRE(nbrs[4].first == 182);
      REQUIRE(nbrs[4].second == 3);
      REQUIRE(nbrs[5].first == 185);
      REQUIRE(nbrs[5].second == 0);
      REQUIRE(nbrs[6].first == 189);
      REQUIRE(nbrs[6].second == 0);
      REQUIRE(nbrs[7].first == 192);
      REQUIRE(nbrs[7].second == 3);
      REQUIRE(nbrs[8].first == 193);
      REQUIRE(nbrs[8].second == 0);
    }
#else
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 4);
      REQUIRE(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      REQUIRE(nbrs[0].first == 160);
      REQUIRE(nbrs[0].second == 0);
      REQUIRE(nbrs[1].first == 163);
      REQUIRE(nbrs[1].second == 0);
      REQUIRE(nbrs[2].first == 170);
      REQUIRE(nbrs[2].second == 0);
      REQUIRE(nbrs[3].first == 180);
      REQUIRE(nbrs[3].second == 2);
      REQUIRE(nbrs[4].first == 182);
      REQUIRE(nbrs[4].second == 3);
      REQUIRE(nbrs[5].first == 185);
      REQUIRE(nbrs[5].second == 0);
      REQUIRE(nbrs[6].first == 189);
      REQUIRE(nbrs[6].second == 0);
      REQUIRE(nbrs[7].first == 192);
      REQUIRE(nbrs[7].second == 3);
      REQUIRE(nbrs[8].first == 193);
      REQUIRE(nbrs[8].second == 0);
    }
#endif
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

TEST_CASE("MultiFPBReader edge cases") {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader edge cases "
      << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/DataStructs/testData/";
  {
    MultiFPBReader mfps;
    mfps.init();
    std::string fps =
        "40081010824820021000500010110410003000402b20285000a4040240010030050000"
        "080001420040009000003d04086007080c03b31d920004220400074008098010206080"
        "00488001080000c64002a00080000200024c2000602410049200340820200002400010"
        "02200106090401056801080182006088101000088a0048";
    ExplicitBitVect qbv(1024);
    UpdateBitVectFromFPSText(qbv, fps);
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv);
      REQUIRE(nbrs.size() == 0);
    }
    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.01);
      REQUIRE(nbrs.size() == 0);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main(int argc, char *argv[]) { return Catch::Session().run(argc, argv); }
