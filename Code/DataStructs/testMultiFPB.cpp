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
#include <DataStructs/MultiFPBReader.h>
#include <DataStructs/BitOps.h>

using namespace RDKit;

void test1MultiFPBReaderBasics() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader basics "
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 2);
    TEST_ASSERT(mfps.getReader(0));
    TEST_ASSERT(mfps.getReader(0)->nBits() == mfps.nBits());
  }
  {
    std::string filename = pathName + "zim.head100.fpb";
    FPBReader fps1(filename), fps2(filename);
    MultiFPBReader mfps;
    TEST_ASSERT(mfps.addReader(&fps1) == 1);
    TEST_ASSERT(mfps.addReader(&fps2) == 2);
    mfps.init();
    TEST_ASSERT(mfps.length() == 2);
    TEST_ASSERT(mfps.getReader(0));
    TEST_ASSERT(mfps.getReader(0)->nBits() == mfps.nBits());
    TEST_ASSERT(mfps.getReader(0));
    TEST_ASSERT(mfps.getReader(0)->nBits() == mfps.nBits());
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
void test2MultiFPBReaderTanimoto() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader Tanimoto "
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 2);

    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.3703));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 0);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 4);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 95);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 95);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.4125));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 89);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 0);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test3MultiFPBReaderTversky() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader Tversky "
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 2);

    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1.);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
    }
    {  // with a threshold
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1., 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.3703));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 0);
    }
    {  // with a threshold, asymmetric
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 0.5, 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 1.));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.4255));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 0);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test4MultiFPBReaderContains() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader contains search"
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 2);
    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(nbrs[0].first == 0);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 0);
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {
      boost::shared_array<std::uint8_t> bytes = mfps.getReader(0)->getBytes(1);
      TEST_ASSERT(bytes);
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 8);
      TEST_ASSERT(nbrs[0].first == 1);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[2].first == 2);
      TEST_ASSERT(nbrs[4].first == 3);
      TEST_ASSERT(nbrs[6].first == 4);
      TEST_ASSERT(nbrs[1].first == 1);
      TEST_ASSERT(nbrs[1].second == 1);
      TEST_ASSERT(nbrs[3].first == 2);
      TEST_ASSERT(nbrs[5].first == 3);
      TEST_ASSERT(nbrs[7].first == 4);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test5MultiFPBReaderThreaded() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader Similarity Threaded"
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 4);
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
      TEST_ASSERT(nbrs.size() == 6);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << std::get<0>(nbrs[i]) << " " <<
      //   std::get<1>(nbrs[i])
      //             << " " << std::get<2>(nbrs[i]) << " " << std::endl;
      // }
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 0.66412));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 3);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 0.65289));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 2);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.64341));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 2);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[3]), 0.61940));
      TEST_ASSERT(std::get<1>(nbrs[3]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[3]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[4]), 0.61905));
      TEST_ASSERT(std::get<1>(nbrs[4]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[4]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[5]), 0.61344));
      TEST_ASSERT(std::get<1>(nbrs[5]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[5]) == 1);
    }

#ifdef RDK_TEST_MULTITHREADED
    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 4);
      TEST_ASSERT(nbrs.size() == 6);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 0.66412));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 3);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 0.65289));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 2);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.64341));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 2);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[3]), 0.61940));
      TEST_ASSERT(std::get<1>(nbrs[3]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[3]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[4]), 0.61905));
      TEST_ASSERT(std::get<1>(nbrs[4]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[4]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[5]), 0.61344));
      TEST_ASSERT(std::get<1>(nbrs[5]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[5]) == 1);
    }
    {  // request more threads than we have readers, this shouldn't be a problem
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 8);
      TEST_ASSERT(nbrs.size() == 6);
      TEST_ASSERT(feq(std::get<0>(nbrs[0]), 0.66412));
      TEST_ASSERT(std::get<1>(nbrs[0]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[0]) == 3);
      TEST_ASSERT(feq(std::get<0>(nbrs[1]), 0.65289));
      TEST_ASSERT(std::get<1>(nbrs[1]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[1]) == 2);
      TEST_ASSERT(feq(std::get<0>(nbrs[2]), 0.64341));
      TEST_ASSERT(std::get<1>(nbrs[2]) == 2);
      TEST_ASSERT(std::get<2>(nbrs[2]) == 1);
      TEST_ASSERT(feq(std::get<0>(nbrs[3]), 0.61940));
      TEST_ASSERT(std::get<1>(nbrs[3]) == 1);
      TEST_ASSERT(std::get<2>(nbrs[3]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[4]), 0.61905));
      TEST_ASSERT(std::get<1>(nbrs[4]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[4]) == 0);
      TEST_ASSERT(feq(std::get<0>(nbrs[5]), 0.61344));
      TEST_ASSERT(std::get<1>(nbrs[5]) == 0);
      TEST_ASSERT(std::get<2>(nbrs[5]) == 1);
    }

#endif
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test6MultiFPBReaderContainsThreaded() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing MultiFPBReader Contains Threaded"
      << std::endl;
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
    TEST_ASSERT(mfps.length() == 4);
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
      TEST_ASSERT(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      TEST_ASSERT(nbrs[0].first == 160);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 163);
      TEST_ASSERT(nbrs[1].second == 0);
      TEST_ASSERT(nbrs[2].first == 170);
      TEST_ASSERT(nbrs[2].second == 0);
      TEST_ASSERT(nbrs[3].first == 180);
      TEST_ASSERT(nbrs[3].second == 2);
      TEST_ASSERT(nbrs[4].first == 182);
      TEST_ASSERT(nbrs[4].second == 3);
      TEST_ASSERT(nbrs[5].first == 185);
      TEST_ASSERT(nbrs[5].second == 0);
      TEST_ASSERT(nbrs[6].first == 189);
      TEST_ASSERT(nbrs[6].second == 0);
      TEST_ASSERT(nbrs[7].first == 192);
      TEST_ASSERT(nbrs[7].second == 3);
      TEST_ASSERT(nbrs[8].first == 193);
      TEST_ASSERT(nbrs[8].second == 0);
    }
#ifdef RDK_TEST_MULTITHREADED
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 4);
      TEST_ASSERT(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      TEST_ASSERT(nbrs[0].first == 160);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 163);
      TEST_ASSERT(nbrs[1].second == 0);
      TEST_ASSERT(nbrs[2].first == 170);
      TEST_ASSERT(nbrs[2].second == 0);
      TEST_ASSERT(nbrs[3].first == 180);
      TEST_ASSERT(nbrs[3].second == 2);
      TEST_ASSERT(nbrs[4].first == 182);
      TEST_ASSERT(nbrs[4].second == 3);
      TEST_ASSERT(nbrs[5].first == 185);
      TEST_ASSERT(nbrs[5].second == 0);
      TEST_ASSERT(nbrs[6].first == 189);
      TEST_ASSERT(nbrs[6].second == 0);
      TEST_ASSERT(nbrs[7].first == 192);
      TEST_ASSERT(nbrs[7].second == 3);
      TEST_ASSERT(nbrs[8].first == 193);
      TEST_ASSERT(nbrs[8].second == 0);
    }
    {  // request more threads than we have readers, this shouldn't be a problem
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 8);
      TEST_ASSERT(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      TEST_ASSERT(nbrs[0].first == 160);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 163);
      TEST_ASSERT(nbrs[1].second == 0);
      TEST_ASSERT(nbrs[2].first == 170);
      TEST_ASSERT(nbrs[2].second == 0);
      TEST_ASSERT(nbrs[3].first == 180);
      TEST_ASSERT(nbrs[3].second == 2);
      TEST_ASSERT(nbrs[4].first == 182);
      TEST_ASSERT(nbrs[4].second == 3);
      TEST_ASSERT(nbrs[5].first == 185);
      TEST_ASSERT(nbrs[5].second == 0);
      TEST_ASSERT(nbrs[6].first == 189);
      TEST_ASSERT(nbrs[6].second == 0);
      TEST_ASSERT(nbrs[7].first == 192);
      TEST_ASSERT(nbrs[7].second == 3);
      TEST_ASSERT(nbrs[8].first == 193);
      TEST_ASSERT(nbrs[8].second == 0);
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
      TEST_ASSERT(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      TEST_ASSERT(nbrs[0].first == 160);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 163);
      TEST_ASSERT(nbrs[1].second == 0);
      TEST_ASSERT(nbrs[2].first == 170);
      TEST_ASSERT(nbrs[2].second == 0);
      TEST_ASSERT(nbrs[3].first == 180);
      TEST_ASSERT(nbrs[3].second == 2);
      TEST_ASSERT(nbrs[4].first == 182);
      TEST_ASSERT(nbrs[4].second == 3);
      TEST_ASSERT(nbrs[5].first == 185);
      TEST_ASSERT(nbrs[5].second == 0);
      TEST_ASSERT(nbrs[6].first == 189);
      TEST_ASSERT(nbrs[6].second == 0);
      TEST_ASSERT(nbrs[7].first == 192);
      TEST_ASSERT(nbrs[7].second == 3);
      TEST_ASSERT(nbrs[8].first == 193);
      TEST_ASSERT(nbrs[8].second == 0);
    }
#else
    {
      std::vector<std::pair<unsigned int, unsigned int>> nbrs =
          mfps.getContainingNeighbors(qbv, 4);
      TEST_ASSERT(nbrs.size() == 9);
      // for (unsigned int i = 0; i < nbrs.size(); ++i) {
      //   std::cerr << i << ": " << nbrs[i].first << " " << nbrs[i].second
      //             << std::endl;
      // }
      TEST_ASSERT(nbrs[0].first == 160);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 163);
      TEST_ASSERT(nbrs[1].second == 0);
      TEST_ASSERT(nbrs[2].first == 170);
      TEST_ASSERT(nbrs[2].second == 0);
      TEST_ASSERT(nbrs[3].first == 180);
      TEST_ASSERT(nbrs[3].second == 2);
      TEST_ASSERT(nbrs[4].first == 182);
      TEST_ASSERT(nbrs[4].second == 3);
      TEST_ASSERT(nbrs[5].first == 185);
      TEST_ASSERT(nbrs[5].second == 0);
      TEST_ASSERT(nbrs[6].first == 189);
      TEST_ASSERT(nbrs[6].second == 0);
      TEST_ASSERT(nbrs[7].first == 192);
      TEST_ASSERT(nbrs[7].second == 3);
      TEST_ASSERT(nbrs[8].first == 193);
      TEST_ASSERT(nbrs[8].second == 0);
    }
#endif
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test7MultiFPBReaderEdges() {
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
      TEST_ASSERT(nbrs.size() == 0);
    }
    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.01);
      TEST_ASSERT(nbrs.size() == 0);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  RDLog::InitLogs();

  test1MultiFPBReaderBasics();
  test2MultiFPBReaderTanimoto();
  test3MultiFPBReaderTversky();
  test4MultiFPBReaderContains();
  test5MultiFPBReaderThreaded();
  test6MultiFPBReaderContainsThreaded();
  test7MultiFPBReaderEdges();

  return 0;
}
