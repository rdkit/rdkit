//  Created by Guillaume GODIN
//  Copyright (C) 2012-2016 Greg Landrum
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
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>

#include <GraphMol/Descriptors/AUTOCORR2D.h>
using namespace RDKit;

void testautocorrelation() {
  BOOST_LOG(rdErrorLog)
      << "----------------------------\nstart test autocorr2D\n";

  std::string pathName = getenv("RDBASE");

  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  // std::string sdfName ="PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName, true, false);

  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/auto2D.out";
  // std::string fName = "auto2D.out";

  std::ifstream instrm(fName.c_str());

  std::string line;
  std::vector<std::vector<std::string>> data;

  while (std::getline(instrm, line)) {
    std::string phrase;
    std::vector<std::string> row;
    std::stringstream ss(line);
    while (std::getline(ss, phrase, '\t')) {
      row.push_back(phrase);
    }
    data.push_back(row);
  }
  int nDone = 0;
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);

    std::vector<double> res2d;

    RDKit::Descriptors::AUTOCORR2D(*m, res2d);

    std::vector<std::string> myrow = data[nDone];
    std::string inm = myrow[0];

    TEST_ASSERT(inm == nm);

    for (int i = 0; i < 192; i++) {
      double ref = atof(myrow[i + 1].c_str());

      if (fabs(ref - res2d[i]) > 0.05) {
        std::cout << "value mismatch: pos" << i << " " << inm << " " << ref
                  << " " << res2d[i] << std::endl;
      }
      TEST_ASSERT(fabs(ref - res2d[i]) < 0.05);
    }

    res2d.clear();
    res2d.resize(48);
    delete m;
    ++nDone;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub3806() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github #3806: NaNs from AUTOCORR2D descriptor" << std::endl;

  {
    auto m = "CC"_smiles;
    std::vector<double> vs;
    RDKit::Descriptors::AUTOCORR2D(*m, vs);
    TEST_ASSERT(vs.size() == 192);
    for (unsigned i = 0; i < vs.size(); ++i) {
      TEST_ASSERT(!std::isnan(vs[i]));
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub9887() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github #9887: AUTOCORR2D values depend on atom order"
      << std::endl;

  {
    std::unique_ptr<ROMol> m("CC(O)C(=O)OCCCCCO"_smiles);
    TEST_ASSERT(m);
    std::vector<double> ref;
    RDKit::Descriptors::AUTOCORR2D(*m, ref);
    TEST_ASSERT(ref.size() == 192);

    // The SMILES from the issue report is a renumbering of the same molecule
    // that used to give different values at positions 144 and 146.
    std::unique_ptr<ROMol> m2("C(COC(C(C)O)=O)CCCO"_smiles);
    TEST_ASSERT(m2);
    std::vector<double> vs;
    RDKit::Descriptors::AUTOCORR2D(*m2, vs);
    TEST_ASSERT(vs.size() == ref.size());
    for (unsigned i = 0; i < ref.size(); ++i) {
      TEST_ASSERT(vs[i] == ref[i]);
    }

    // Several deterministic renumberings of the same molecule.
    const unsigned int numAtoms = m->getNumAtoms();
    std::vector<unsigned int> order(numAtoms);
    std::iota(order.begin(), order.end(), 0u);
    std::mt19937 rng(42);
    for (unsigned trial = 0; trial < 20; ++trial) {
      std::shuffle(order.begin(), order.end(), rng);
      std::unique_ptr<ROMol> shuffled(MolOps::renumberAtoms(*m, order));
      TEST_ASSERT(shuffled);
      std::vector<double> shuffledVs;
      RDKit::Descriptors::AUTOCORR2D(*shuffled, shuffledVs);
      TEST_ASSERT(shuffledVs.size() == ref.size());
      for (unsigned i = 0; i < ref.size(); ++i) {
        TEST_ASSERT(shuffledVs[i] == ref[i]);
      }
    }

    // The custom-property variant has the same order dependence; check it
    // too. The property is a per-atom value (atomic number), set after each
    // renumbering so it follows the atoms.
    auto setCustomProp = [](ROMol *mol) {
      for (auto atom : mol->atoms()) {
        atom->setProp("customProp",
                      static_cast<double>(atom->getAtomicNum()));
      }
    };
    setCustomProp(m.get());
    std::vector<double> refCustom;
    RDKit::Descriptors::AUTOCORR2D(*m, refCustom, "customProp");
    TEST_ASSERT(refCustom.size() == 32);

    std::iota(order.begin(), order.end(), 0u);
    std::mt19937 rng2(123);
    for (unsigned trial = 0; trial < 20; ++trial) {
      std::shuffle(order.begin(), order.end(), rng2);
      std::unique_ptr<ROMol> shuffled(MolOps::renumberAtoms(*m, order));
      TEST_ASSERT(shuffled);
      setCustomProp(shuffled.get());
      std::vector<double> shuffledVs;
      RDKit::Descriptors::AUTOCORR2D(*shuffled, shuffledVs, "customProp");
      TEST_ASSERT(shuffledVs.size() == refCustom.size());
      for (unsigned i = 0; i < refCustom.size(); ++i) {
        TEST_ASSERT(shuffledVs[i] == refCustom[i]);
      }
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  testautocorrelation();
  testGithub3806();
  testGithub9887();
}
