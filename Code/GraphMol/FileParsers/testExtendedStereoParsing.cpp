//
//  Copyright (C) 2018 T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <clocale>
#include <cstdlib>

#include <string>
#include <fstream>
#include <memory>

using namespace RDKit;

std::unique_ptr<RWMol> readTestFile(const std::string &baseName) {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/" + baseName;
  auto m = MolFileToMol(fName);
  return std::unique_ptr<RWMol>(m);
}

void testOr() {
  BOOST_LOG(rdInfoLog) << "testing extended stereo parsing with an OR block"
                       << std::endl;

  auto m = readTestFile("two_centers_or.mol");
  TEST_ASSERT(m.get());
  TEST_ASSERT(m->getNumAtoms() == 8);

  auto stereo_groups = m->getStereoGroups();
  TEST_ASSERT(stereo_groups.size() == 2);
  TEST_ASSERT(stereo_groups[0].getGroupType() ==
              RDKit::StereoGroupType::STEREO_ABSOLUTE);
  TEST_ASSERT(stereo_groups[0].getAtoms().size() == 1u);
  TEST_ASSERT(stereo_groups[1].getGroupType() ==
              RDKit::StereoGroupType::STEREO_OR);
  TEST_ASSERT(stereo_groups[1].getAtoms().size() == 2u);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAnd() {
  BOOST_LOG(rdInfoLog) << "testing extended stereo parsing with an AND block"
                       << std::endl;

  auto m = readTestFile("two_centers_and.mol");
  TEST_ASSERT(m.get());
  TEST_ASSERT(m->getNumAtoms() == 8);

  auto stereo_groups = m->getStereoGroups();
  TEST_ASSERT(stereo_groups.size() == 2);
  TEST_ASSERT(stereo_groups[0].getGroupType() ==
              RDKit::StereoGroupType::STEREO_ABSOLUTE);
  TEST_ASSERT(stereo_groups[1].getGroupType() ==
              RDKit::StereoGroupType::STEREO_AND);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testWrite() {
  BOOST_LOG(rdInfoLog) << "testing extended stereo file writing" << std::endl;

  auto m0 = readTestFile("two_centers_and.mol");
  TEST_ASSERT(m0.get());
  std::string block = RDKit::MolToMolBlock(*m0);
  auto m1 = RDKit::MolBlockToMol(block);

  // Check that the extended stereo information has the same extended stereo
  // types and same atoms marked for extended stereo.
  auto stereo_groups0 = m0->getStereoGroups();
  auto stereo_groups1 = m1->getStereoGroups();
  TEST_ASSERT(stereo_groups0.size() == stereo_groups1.size());

  for (unsigned i = 0u; i < 2; ++i) {
    TEST_ASSERT(stereo_groups0[i].getGroupType() ==
                stereo_groups1[i].getGroupType());
    TEST_ASSERT(stereo_groups0[i].getAtoms().size() ==
                stereo_groups1[i].getAtoms().size());
    for (auto &&atom0 = stereo_groups0[i].getAtoms().begin(),
              atom1 = stereo_groups1[i].getAtoms().begin();
         atom0 != stereo_groups0[i].getAtoms().end(); ++atom0, ++atom1) {
      TEST_ASSERT((*atom0)->getIdx() == (*atom1)->getIdx());
    }
  }
  delete (m1);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  testOr();
  testAnd();
  testWrite();

  return 0;
}
