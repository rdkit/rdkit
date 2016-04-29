//
//  Copyright (C) 2016 Greg Landrum
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
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 0);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 0);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.3703));
      TEST_ASSERT(nbrs[2].get<1>() == 1);
      TEST_ASSERT(nbrs[2].get<2>() == 0);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(95);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(bytes, 0.30);
      TEST_ASSERT(nbrs.size() == 4);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 95);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 95);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.4125));
      TEST_ASSERT(nbrs[2].get<1>() == 89);
      TEST_ASSERT(nbrs[2].get<2>() == 0);
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
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1.);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 0);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
    }
    {  // with a threshold
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 1., 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 0);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.3703));
      TEST_ASSERT(nbrs[2].get<1>() == 1);
      TEST_ASSERT(nbrs[2].get<2>() == 0);
    }
    {  // with a threshold, asymmetric
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTverskyNeighbors(bytes, 1., 0.5, 0.30);
      TEST_ASSERT(nbrs.size() == 10);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 1.));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 1.));
      TEST_ASSERT(nbrs[1].get<1>() == 0);
      TEST_ASSERT(nbrs[1].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.4255));
      TEST_ASSERT(nbrs[2].get<1>() == 1);
      TEST_ASSERT(nbrs[2].get<2>() == 0);
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
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(0);
      TEST_ASSERT(bytes);
      std::vector<std::pair<unsigned int, unsigned int> > nbrs =
          mfps.getContainingNeighbors(bytes);
      TEST_ASSERT(nbrs.size() == 2);
      TEST_ASSERT(nbrs[0].first == 0);
      TEST_ASSERT(nbrs[0].second == 0);
      TEST_ASSERT(nbrs[1].first == 0);
      TEST_ASSERT(nbrs[1].second == 1);
    }
    {
      boost::shared_array<boost::uint8_t> bytes =
          mfps.getReader(0)->getBytes(1);
      TEST_ASSERT(bytes);
      std::vector<std::pair<unsigned int, unsigned int> > nbrs =
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
      << "-----------------------\n Testing MultiFPBReader Tanimoto Threaded"
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
      //   std::cerr << i << ": " << nbrs[i].get<0>() << " " << nbrs[i].get<1>()
      //             << " " << nbrs[i].get<2>() << " " << std::endl;
      // }
      TEST_ASSERT(feq(nbrs[0].get<0>(), 0.66412));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 3);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 0.65289));
      TEST_ASSERT(nbrs[1].get<1>() == 1);
      TEST_ASSERT(nbrs[1].get<2>() == 2);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.64341));
      TEST_ASSERT(nbrs[2].get<1>() == 2);
      TEST_ASSERT(nbrs[2].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[3].get<0>(), 0.61940));
      TEST_ASSERT(nbrs[3].get<1>() == 1);
      TEST_ASSERT(nbrs[3].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[4].get<0>(), 0.61905));
      TEST_ASSERT(nbrs[4].get<1>() == 0);
      TEST_ASSERT(nbrs[4].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[5].get<0>(), 0.61344));
      TEST_ASSERT(nbrs[5].get<1>() == 0);
      TEST_ASSERT(nbrs[5].get<2>() == 1);
    }

#ifdef RDK_TEST_MULTITHREADED
    {
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 4);
      TEST_ASSERT(nbrs.size() == 6);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 0.66412));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 3);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 0.65289));
      TEST_ASSERT(nbrs[1].get<1>() == 1);
      TEST_ASSERT(nbrs[1].get<2>() == 2);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.64341));
      TEST_ASSERT(nbrs[2].get<1>() == 2);
      TEST_ASSERT(nbrs[2].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[3].get<0>(), 0.61940));
      TEST_ASSERT(nbrs[3].get<1>() == 1);
      TEST_ASSERT(nbrs[3].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[4].get<0>(), 0.61905));
      TEST_ASSERT(nbrs[4].get<1>() == 0);
      TEST_ASSERT(nbrs[4].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[5].get<0>(), 0.61344));
      TEST_ASSERT(nbrs[5].get<1>() == 0);
      TEST_ASSERT(nbrs[5].get<2>() == 1);
    }
    {  // request more threads than we have readers, this shouldn't be a problem
      std::vector<MultiFPBReader::ResultTuple> nbrs =
          mfps.getTanimotoNeighbors(qbv, 0.6, 8);
      TEST_ASSERT(nbrs.size() == 6);
      TEST_ASSERT(feq(nbrs[0].get<0>(), 0.66412));
      TEST_ASSERT(nbrs[0].get<1>() == 0);
      TEST_ASSERT(nbrs[0].get<2>() == 3);
      TEST_ASSERT(feq(nbrs[1].get<0>(), 0.65289));
      TEST_ASSERT(nbrs[1].get<1>() == 1);
      TEST_ASSERT(nbrs[1].get<2>() == 2);
      TEST_ASSERT(feq(nbrs[2].get<0>(), 0.64341));
      TEST_ASSERT(nbrs[2].get<1>() == 2);
      TEST_ASSERT(nbrs[2].get<2>() == 1);
      TEST_ASSERT(feq(nbrs[3].get<0>(), 0.61940));
      TEST_ASSERT(nbrs[3].get<1>() == 1);
      TEST_ASSERT(nbrs[3].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[4].get<0>(), 0.61905));
      TEST_ASSERT(nbrs[4].get<1>() == 0);
      TEST_ASSERT(nbrs[4].get<2>() == 0);
      TEST_ASSERT(feq(nbrs[5].get<0>(), 0.61344));
      TEST_ASSERT(nbrs[5].get<1>() == 0);
      TEST_ASSERT(nbrs[5].get<2>() == 1);
    }

#endif
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
  return 0;
}
