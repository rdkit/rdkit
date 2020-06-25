//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/test.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "GeneralFileReader.h"

namespace io = boost::iostreams;
using namespace RDKit;
using namespace GeneralMolSupplier;
void testSdf() {
  //! Open uncompressed SDF file format
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  struct SupplierOptions opt;
  MolSupplier* sdsup = getSupplier(fname, opt);
  unsigned int i = 0;
  while (!sdsup->atEnd()) {
    ROMol* nmol = sdsup->next();
    TEST_ASSERT(nmol || sdsup->atEnd());
    if (nmol) {
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));
      TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 16);

  //! Open compressed SDF file format
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf.gz";
  opt.takeOwnership = false;
  MolSupplier* sdsup2 = getSupplier(fname, opt);
  i = 0;
  while (!sdsup2->atEnd()) {
    ROMol* nmol = sdsup2->next();
    if (nmol) {
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));
      TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 16);

  delete sdsup;
  delete sdsup2;
}

void testSmi() {
  //! Open uncompressed SMI file format, try .csv formats
  std::string mname;
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  struct SupplierOptions opt_smi;
  opt_smi.takeOwnership = false;
  opt_smi.delimiter = ",";
  opt_smi.smilesColumn = 1;
  opt_smi.nameColumn = 0;
  opt_smi.titleLine = true;
  MolSupplier* sup = getSupplier(fname, opt_smi);
  unsigned int i = 0;
  while (!sup->atEnd()) {
    ROMol* mol = sup->next();
    if (i == 3) {
      mol->getProp(common_properties::_Name, mname);
      CHECK_INVARIANT(mname == "4", "");
      mol->getProp("TPSA", mname);
      CHECK_INVARIANT(mname == "82.78", "");
    }
    delete mol;
    i++;
  }
  TEST_ASSERT(i == 10);
  delete sup;
}

void testMae() {
  //! Open uncompressed MAE file format
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/props_test.mae";
  struct SupplierOptions opt;
  MolSupplier* maesup = getSupplier(fname, opt);
  std::unique_ptr<ROMol> nmol(maesup->next());
  TEST_ASSERT(nmol);

  //! Test mol properties
  TEST_ASSERT(nmol->hasProp(common_properties::_Name));
  TEST_ASSERT(nmol->hasProp("b_sd_chiral_flag"));
  TEST_ASSERT(nmol->getProp<bool>("b_sd_chiral_flag") == false);
  TEST_ASSERT(nmol->hasProp("i_sd_NSC"));
  TEST_ASSERT(nmol->getProp<int>("i_sd_NSC") == 48);
  TEST_ASSERT(nmol->hasProp("s_m_entry_name"));
  TEST_ASSERT(nmol->getProp<std::string>("s_m_entry_name") == "NCI_aids_few.1");
  TEST_ASSERT(nmol->hasProp("r_f3d_dummy"));
  TEST_ASSERT(abs(nmol->getProp<double>("r_f3d_dummy") - 42.123) < 0.0001);

  //! Test atom properties
  TEST_ASSERT(nmol->getNumAtoms() == 19);
  for (int i = 0; i < 19; ++i) {
    const auto* atom = nmol->getAtomWithIdx(i);

    //! The integer property is present for all atoms
    TEST_ASSERT(atom->hasProp("i_m_minimize_atom_index"));
    TEST_ASSERT(atom->getProp<int>("i_m_minimize_atom_index") == 1 + i);

    //! The bool property is only defined for i < 10
    if (i < 10) {
      TEST_ASSERT(atom->hasProp("b_m_dummy"));
      TEST_ASSERT(atom->getProp<bool>("b_m_dummy") == static_cast<bool>(i % 2));
    } else {
      TEST_ASSERT(!atom->hasProp("b_m_dummy"));
    }

    //! The real property is only defined for i >= 10
    if (i >= 10) {
      TEST_ASSERT(atom->hasProp("r_f3d_dummy"));
      TEST_ASSERT(abs(atom->getProp<double>("r_f3d_dummy") - (19.1 - i)) <
                  0.0001);
    } else {
      TEST_ASSERT(!atom->hasProp("r_f3d_dummy"));
    }

    //! All atoms have the string prop
    TEST_ASSERT(atom->hasProp("s_m_dummy"));
    TEST_ASSERT(atom->getProp<std::string>("s_m_dummy") ==
                std::to_string(19 - i));
  }
  TEST_ASSERT(maesup->atEnd());

  //! Open compressed MAE file, .maegz format
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1kv1.maegz";
  MolSupplier* cmaesup = getSupplier(fname, opt);

  std::shared_ptr<ROMol> nmol2;
  nmol2.reset(cmaesup->next());
  const Atom* atom = nmol2->getAtomWithIdx(0);
  auto* info = (AtomPDBResidueInfo*)(atom->getMonomerInfo());
  TEST_ASSERT(info->getResidueName() == "ARG ");
  TEST_ASSERT(info->getChainId() == "A");
  TEST_ASSERT(info->getResidueNumber() == 5);

  delete maesup;
  delete cmaesup;
}

void testTdt() {
  //! Open uncompressed TDT file
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  struct SupplierOptions opt;
  opt.nameRecord = "PN";
  MolSupplier* suppl = getSupplier(fname, opt);

  unsigned int i = 0;
  while (!suppl->atEnd()) {
    ROMol* nmol = suppl->next();
    if (nmol) {
      std::string prop1, prop2;
      TEST_ASSERT(nmol->getNumAtoms() > 0);
      TEST_ASSERT(nmol->hasProp("PN"));
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));
      TEST_ASSERT(nmol->hasProp("MFCD"));

      nmol->getProp("PN", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      TEST_ASSERT(prop1 == prop2);

      //! we didn't ask for 2D conformers, so there should be a property 2D:
      TEST_ASSERT(nmol->hasProp(common_properties::TWOD));
      //! and no conformer:
      TEST_ASSERT(!nmol->getNumConformers());

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 10);
  delete suppl;
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSdf();
  BOOST_LOG(rdErrorLog) << "Finished: testSdf()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSmi();
  BOOST_LOG(rdErrorLog) << "Finished: testSmi()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testMae();
  BOOST_LOG(rdErrorLog) << "Finished: testMae()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testTdt();
  BOOST_LOG(rdErrorLog) << "Finished: testTdt()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  return 0;
}
