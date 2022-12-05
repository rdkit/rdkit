//
//  Copyright (C) 2020 Brian P. Kelley
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/RDLog.h>
#include <vector>

#include <GraphMol/Descriptors/BCUT.h>
using namespace RDKit;

// from knuth.
bool feq(float a, float b, float epsilon = 1e-6) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic BCUT tests." << std::endl;

  std::vector<std::string> vec = {"c1ccccc1S", "c1cscc1",    "CC(=S)C",
                                  "CSC",       "CS(=O)C",    "CP(C)C",
                                  "CP=O",      "CP(C)(C)=O", "C[PH](C)=O"};
  std::vector<std::vector<double>> expected = {
      {32.11691025659743, 10.3711255714102, 1.7258015589384423,
       -1.813025067747632, 2.0032623406582895, -1.5648237280483932,
       7.798895708674262, 1.480463460412681},
      {32.133385673916855, 10.690498814349603, 1.5812862454321042,
       -1.3924158848795094, 1.9109090827066813, -1.1624011856503906,
       7.073072952311218, 2.0294988143496093},
      {32.09203038492705, 10.585642854940236, 1.5531725323483636,
       -1.6095571920399787, 1.6095355058255407, -1.6264025009240035,
       7.798379459554676, 1.8522004089765778},
      {32.16623472793675, 11.912765272063261, 1.3229624149945334,
       -1.5095023144073105, 1.699644712375421, -1.2539447123754206,
       7.974132068991219, 2.370867931008784},
      {32.19696040763111, 11.912196087574864, 1.571634336328628,
       -1.5986639964786287, 1.4708186703193278, -1.7011978973115978,
       7.829987043024966, 0.6900118228151252},
      {31.130920829765557, 11.856079170234391, 1.652795743173261,
       -1.8120831465729346, 2.1415709560249625, -1.4818709560249614,
       7.546153955380794, 2.12884604461921},
      {31.059700934072303, 11.958086319851176, 1.271920282698931,
       -1.2044337292760101, 1.6734005705339297, -0.9875941256320955,
       7.220415332975225, 0.601201124030627},
      {31.163638790374048, 11.855143433013291, 1.8876168729214116,
       -1.862911092419536, 2.279553713445213, -1.6070814330824321,
       7.610231826377116, 0.5859042224902173},
      {31.111807418763107, 11.906469225329653, 1.5848247701762241,
       -1.5909887821270712, 2.0038391834918317, -1.3277290609201042,
       7.422700514979134, 0.5941465179281992}};

  for (size_t i = 0; i < vec.size(); ++i) {
    std::unique_ptr<ROMol> m(SmilesToMol(vec[i]));
    auto bcut = Descriptors::BCUT2D(*m);
    TEST_ASSERT(bcut.size() == expected[i].size());
    for (size_t j = 0; j < bcut.size(); ++j) {
      TEST_ASSERT(feq(bcut[j], expected[i][j]));
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Larger BCUT tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/bace_200_bcut.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);
  while (!reader.atEnd()) {
    std::unique_ptr<RDKit::ROMol> m(reader.next());
    TEST_ASSERT(m.get());
    std::vector<double> bcuts = Descriptors::BCUT2D(*m);
    TEST_ASSERT(bcuts.size() == 8);
    TEST_ASSERT(feq(bcuts[0], m->getProp<double>("bcut1")));
    TEST_ASSERT(feq(bcuts[1], m->getProp<double>("bcut2")));
    TEST_ASSERT(feq(bcuts[2], m->getProp<double>("bcut3")));
    TEST_ASSERT(feq(bcuts[3], m->getProp<double>("bcut4")));
    TEST_ASSERT(feq(bcuts[4], m->getProp<double>("bcut5")));
    TEST_ASSERT(feq(bcuts[5], m->getProp<double>("bcut6")));
    TEST_ASSERT(feq(bcuts[6], m->getProp<double>("bcut7")));
    TEST_ASSERT(feq(bcuts[7], m->getProp<double>("bcut8")));
  }
}

int main() {
  RDLog::InitLogs();
  test1();
  test2();
}
