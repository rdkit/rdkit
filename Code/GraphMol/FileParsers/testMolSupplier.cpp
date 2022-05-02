//
//   Copyright (C) 2002-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>

#include "MolSupplier.h"
#include "MolWriters.h"
#include "FileParsers.h"
#include "FileParserUtils.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/RDLog.h>
#include <RDStreams/streams.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
namespace io = boost::iostreams;

using namespace RDKit;

int testMolSup() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";

  {
    SDMolSupplier sdsup(fname);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i == 16);
  }
  {
    SDMolSupplier sdsup(fname);
    for (unsigned int i = 0; i < 16; ++i) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
    }
    // test issue 3524949:
    TEST_ASSERT(sdsup.atEnd());
    bool ok = false;
    try {
      sdsup.next();
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    std::ifstream strm(fname.c_str());
    SDMolSupplier sdsup(&strm, false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i == 16);
  }
  {
    auto *strm = new std::ifstream(fname.c_str());
    SDMolSupplier sdsup(strm, true);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    TEST_ASSERT(i == 16);
  }
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  {  // Test reading properties
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/props_test.mae";

    MaeMolSupplier maesup(fname);
    std::unique_ptr<ROMol> nmol(maesup.next());
    TEST_ASSERT(nmol);

    // Test mol properties
    TEST_ASSERT(nmol->hasProp(common_properties::_Name));
    TEST_ASSERT(nmol->hasProp("b_sd_chiral_flag"));
    TEST_ASSERT(nmol->getProp<bool>("b_sd_chiral_flag") == false);
    TEST_ASSERT(nmol->hasProp("i_sd_NSC"));
    TEST_ASSERT(nmol->getProp<int>("i_sd_NSC") == 48);
    TEST_ASSERT(nmol->hasProp("s_m_entry_name"));
    TEST_ASSERT(nmol->getProp<std::string>("s_m_entry_name") ==
                "NCI_aids_few.1");
    TEST_ASSERT(nmol->hasProp("r_f3d_dummy"));
    TEST_ASSERT(std::abs(nmol->getProp<double>("r_f3d_dummy") - 42.123) <
                0.0001);

    // Test atom properties
    TEST_ASSERT(nmol->getNumAtoms() == 19);
    for (int i = 0; i < 19; ++i) {
      const auto *atom = nmol->getAtomWithIdx(i);

      // The integer property is present for all atoms
      TEST_ASSERT(atom->hasProp("i_m_minimize_atom_index"));
      TEST_ASSERT(atom->getProp<int>("i_m_minimize_atom_index") == 1 + i);

      // The bool property is only defined for i < 10
      if (i < 10) {
        TEST_ASSERT(atom->hasProp("b_m_dummy"));
        TEST_ASSERT(atom->getProp<bool>("b_m_dummy") ==
                    static_cast<bool>(i % 2));
      } else {
        TEST_ASSERT(!atom->hasProp("b_m_dummy"));
      }

      // The real property is only defined for i >= 10
      if (i >= 10) {
        TEST_ASSERT(atom->hasProp("r_f3d_dummy"));
        TEST_ASSERT(std::abs(atom->getProp<double>("r_f3d_dummy") -
                             (19.1 - i)) < 0.0001);
      } else {
        TEST_ASSERT(!atom->hasProp("r_f3d_dummy"));
      }

      // All atoms have the string prop
      TEST_ASSERT(atom->hasProp("s_m_dummy"));
      TEST_ASSERT(atom->getProp<std::string>("s_m_dummy") ==
                  std::to_string(19 - i));
    }

    TEST_ASSERT(maesup.atEnd());
  }
  {  // Test parsing stereo properties. Mol is 2D and has stereo labels.
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/stereochem.mae";
    MaeMolSupplier maesup(fname);

    {  // Stereo bonds. These get overwritten by the double bond detection.
      std::unique_ptr<ROMol> nmol(maesup.next());
      TEST_ASSERT(nmol);
      {
        Bond *bnd = nmol->getBondWithIdx(1);
        TEST_ASSERT(bnd);
        TEST_ASSERT(bnd->getStereoAtoms() == INT_VECT({0, 3}));
        TEST_ASSERT(bnd->getStereo() == Bond::STEREOTRANS);
      }
      {
        Bond *bnd = nmol->getBondWithIdx(3);
        TEST_ASSERT(bnd);
        TEST_ASSERT(bnd->getStereoAtoms() == INT_VECT({2, 5}));
        TEST_ASSERT(bnd->getStereo() == Bond::STEREOCIS);
      }
    }
    {  // Chiralities (these get CIP codes)
      std::unique_ptr<ROMol> nmol(maesup.next());
      TEST_ASSERT(nmol);
      {
        Atom *at = nmol->getAtomWithIdx(1);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
        TEST_ASSERT(at->getProp<std::string>(common_properties::_CIPCode) ==
                    "R");
      }
      {
        Atom *at = nmol->getAtomWithIdx(3);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
        TEST_ASSERT(at->getProp<std::string>(common_properties::_CIPCode) ==
                    "S");
      }
    }
    {  // Pseudochiralities (no CIP codes)
      std::unique_ptr<ROMol> nmol(maesup.next());
      TEST_ASSERT(nmol);
      {
        Atom *at = nmol->getAtomWithIdx(2);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
        TEST_ASSERT(!at->hasProp(common_properties::_CIPCode));
      }
      {
        Atom *at = nmol->getAtomWithIdx(5);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
        TEST_ASSERT(!at->hasProp(common_properties::_CIPCode));
      }
    }
    {  // intentionally bad chirality label, intended to
      // make sure we can step over parse errors
      std::unique_ptr<ROMol> nmol;
      try {
        nmol.reset(maesup.next());
      } catch (const Invar::Invariant &) {
        // just ignore this failure
      }
      TEST_ASSERT(!nmol);
    }
    {  // "Undefined" chirality label
      std::unique_ptr<ROMol> nmol(maesup.next());
      TEST_ASSERT(nmol);
      {
        Atom *at = nmol->getAtomWithIdx(2);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_UNSPECIFIED);
        TEST_ASSERT(!at->hasProp(common_properties::_CIPCode));
      }
      {
        Atom *at = nmol->getAtomWithIdx(5);
        TEST_ASSERT(at);
        TEST_ASSERT(at->getChiralTag() == Atom::CHI_UNSPECIFIED);
        TEST_ASSERT(!at->hasProp(common_properties::_CIPCode));
      }
    }
    TEST_ASSERT(maesup.atEnd());
  }
  {  // Test loop reading
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.mae";
    MaeMolSupplier maesup(fname);
    std::shared_ptr<ROMol> nmol;
    for (unsigned int i = 0; i < 16; ++i) {
      nmol.reset(maesup.next());
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->getNumAtoms() > 0);
        if (i == 0) {
          auto smiles = MolToSmiles(*nmol);
          TEST_ASSERT(smiles ==
                      "CCC1=[O+][Cu@]2([O+]=C(CC)CC(CC)=[O+]2)[O+]=C(CC)C1");
        }
      }
    }
    TEST_ASSERT(maesup.atEnd());
    bool ok = false;
    try {
      maesup.next();
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/bad_ppty.mae";
    const std::string err_msg_substr = "Bad format for property";

    bool ok = false;
    std::unique_ptr<ROMol> mol;
    MaeMolSupplier maesup(fname);

    // This is in excess: there are only 3 mols in the file, and the second one
    // has an invalid property name, so it won't be read
    for (unsigned int i = 0; i < 5; ++i) {
      try {
        mol.reset(maesup.next());
      } catch (const FileParseException &e) {
        const std::string err_msg(e.what());
        TEST_ASSERT(i == 1);
        TEST_ASSERT(err_msg.find(err_msg_substr) != std::string::npos);
        ok = true;
        break;
      }
      TEST_ASSERT(mol);
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->getNumAtoms() == 1);
      TEST_ASSERT(!maesup.atEnd());
    }
    TEST_ASSERT(!maesup.atEnd());
    TEST_ASSERT(ok);
  }

#if RDK_USE_BOOST_IOSTREAMS
  {  // Test Maestro PDB property reading
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1kv1.maegz";
    auto *strm = new gzstream(fname);
    MaeMolSupplier maesup(strm);

    std::shared_ptr<ROMol> nmol;
    nmol.reset(maesup.next());
    const Atom *atom = nmol->getAtomWithIdx(0);
    auto *info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
    TEST_ASSERT(info->getResidueName() == "ARG ");
    TEST_ASSERT(info->getChainId() == "A");
    TEST_ASSERT(info->getResidueNumber() == 5);
  }
#endif
#endif  // RDK_BUILD_MAEPARSER_SUPPORT
  return 1;
}

void testRandMolSup() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  // std::string fname("../test_data/NCI_aids_few.sdf");
  SDMolSupplier sdsup(fname);

  ROMol *tmol = sdsup[7];
  delete tmol;

  CHECK_INVARIANT(sdsup.length() == 16, "");

  STR_VECT names;
  names.push_back(std::string("48"));
  names.push_back(std::string("128"));
  names.push_back(std::string("164"));
  names.push_back(std::string("180"));
  names.push_back(std::string("192"));
  names.push_back(std::string("210"));
  names.push_back(std::string("213"));
  names.push_back(std::string("229"));

  int i;
  for (i = 0; i < 8; i++) {
    ROMol *mol = sdsup[2 * i];
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == names[i], "");
    delete mol;
  }

  // get a random molecule
  ROMol *mol = sdsup[5];
  TEST_ASSERT(mol);
  std::string mname;
  mol->getProp(common_properties::_Name, mname);
  delete mol;
  CHECK_INVARIANT(mname == "170", "");

  // get the last molecule:
  mol = sdsup[15];
  TEST_ASSERT(mol);
  delete mol;

  // and make sure we're at the end:
  TEST_ASSERT(sdsup.atEnd());
  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  mol = sdsup[0];
  TEST_ASSERT(mol);
  delete mol;

  // Issue 113: calling length before grabbing a molecule results in crashes:
  SDMolSupplier sdsup2(fname);
  CHECK_INVARIANT(sdsup2.length() == 16, "");
}

void testSmilesSup() {
  std::string mname;
  std::string fname;
  ROMol *mol;

  std::string rdbase = getenv("RDBASE");
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);
    TEST_ASSERT(nSup2.length() == 10);
  }
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);

    mol = nSup2[3];
    TEST_ASSERT(!nSup2.atEnd())
    TEST_ASSERT(nSup2.length() == 10);

    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == "4", "");
    mol->getProp("TPSA", mname);
    CHECK_INVARIANT(mname == "82.78", "");
    delete mol;

    mol = nSup2[9];
    TEST_ASSERT(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    TEST_ASSERT(mol);
    delete mol;
  }
  {
    std::ifstream strm(fname.c_str(), std::ios_base::binary);
    SmilesMolSupplier nSup2(&strm, false, ",", 1, 0, true);

    mol = nSup2[3];
    CHECK_INVARIANT(nSup2.length() == 10, "");

    mol->getProp(common_properties::_Name, mname);
    CHECK_INVARIANT(mname == "4", "");
    mol->getProp("TPSA", mname);
    CHECK_INVARIANT(mname == "82.78", "");
    delete mol;

    mol = nSup2[9];
    TEST_ASSERT(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    TEST_ASSERT(mol);
    delete mol;
  }

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  SmilesMolSupplier smiSup(fname, ",", 0, -1);

  mol = smiSup[16];

  mol->getProp("TPSA", mname);
  CHECK_INVARIANT(mname == "46.25", "");
  delete mol;

  mol = smiSup[8];
  mol->getProp("TPSA", mname);
  CHECK_INVARIANT(mname == "65.18", "");
  delete mol;

  int len = smiSup.length();
  CHECK_INVARIANT(len == 200, "");

  smiSup.reset();
  int i = 0;
  mol = smiSup.next();
  while (1) {
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    i++;
    delete mol;
    try {
      mol = smiSup.next();
    } catch (FileParseException &) {
      break;
    }
  }

  CHECK_INVARIANT(i == 200, "");

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);

  // check the length before we read anything out...
  //  this was a problem at one point (Issue 113)
  CHECK_INVARIANT(nSup->length() == 10, "");
  mol = (*nSup)[3];

  mol->getProp(common_properties::_Name, mname);
  CHECK_INVARIANT(mname == "4", "");
  mol->getProp("Column_2", mname);
  CHECK_INVARIANT(mname == "82.78", "");

  delete nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  unsigned int nRead = 0;
  while (!nSup->atEnd()) {
    delete mol;
    mol = nSup->next();
    TEST_ASSERT(mol);
    nRead++;
  }
  TEST_ASSERT(nSup->length() == 10);
  TEST_ASSERT(nRead == 10);

  delete nSup;
  delete mol;
}

void testSmilesSupFromText() {
  std::string mname;
  std::string fname;
  ROMol *mol;

  SmilesMolSupplier nSup2;
  std::string text;
  bool failed;
  int nAts;

  // this was a delightful boundary condition:
  BOOST_LOG(rdErrorLog)
      << "------------------------------------------------------" << std::endl;
  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  {
    nSup2.setData(text, " ", 0, -1, false, true);
    //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
    mol = nSup2.next();
    nAts = mol->getNumAtoms();
    delete mol;
    TEST_ASSERT(nAts == 2);

    mol = nSup2[3];
    nAts = mol->getNumAtoms();
    delete mol;
    TEST_ASSERT(nAts == 6);
    TEST_ASSERT(nSup2.length() == 4);

    failed = false;
    try {
      mol = nSup2[4];
      delete mol;
    } catch (FileParseException &) {
      failed = true;
    }
    TEST_ASSERT(failed);
    mol = nSup2[2];
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts == 4);
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    TEST_ASSERT(mname == "2");
    delete mol;
  }
  {
    nSup2.setData(text, " ", 0, -1, false, true);
    mol = nSup2[2];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts == 4);
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    TEST_ASSERT(mname == "2");
    delete mol;

    mol = nSup2[3];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts == 6);
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    TEST_ASSERT(mname == "3");
    delete mol;
  }
  {
    nSup2.setData(text, " ", 0, -1, false, true);
    mol = nSup2[3];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts == 6);
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    TEST_ASSERT(mname == "3");

    delete mol;
    mol = nSup2[2];
    TEST_ASSERT(mol);
    nAts = mol->getNumAtoms();
    TEST_ASSERT(nAts == 4);
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    TEST_ASSERT(mname == "2");
    delete mol;
  }
  // --------------
  // basics:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-3 CCC 9.0\n"
      "mol-4 CCCC 16.0\n";
#if 1
  nSup2.setData(text, " ", 1, 0, true, true);
  mol = nSup2[3];
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  CHECK_INVARIANT(nSup2.length() == 4, "");
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-4");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "16.0");
  delete mol;

  // ensure that we can call setData a second time:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-3 CCC 9.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);
  CHECK_INVARIANT(nSup2.length() == 3, "");
  mol = nSup2[2];
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-3");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "9.0");
  delete mol;

  // now test for failure handling:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-3 fail 9.0\n"
      "mol-4 CCCC 16.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);
  mol = nSup2[3];
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 4);
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-4");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "16.0");
  delete mol;

  // failures should give null molecules:
  mol = nSup2[2];
  TEST_ASSERT(!mol);
  delete mol;
#endif

  // issue 114, no \n at EOF:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-4");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "16.0");
  TEST_ASSERT(nSup2.atEnd());
  delete mol;

  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0";
  nSup2.setData(text, " ", 1, 0, true, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-4");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "16.0");
  TEST_ASSERT(nSup2.atEnd());
  delete mol;

  try {
    mol = nSup2[3];
    delete mol;
  } catch (FileParseException &) {
    failed = true;
  }
  TEST_ASSERT(failed);

  text =
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0";
  nSup2.setData(text, " ", 1, 0, false, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-4");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "16.0");
  delete mol;

  text =
      "C\n"
      "CC\n"
      "CCCC";
  nSup2.setData(text, " ", 0, -1, false, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  TEST_ASSERT(nSup2.length() == 3);
  mol = nSup2[2];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  delete mol;

  // this was a delightful boundary condition:
  BOOST_LOG(rdErrorLog)
      << "------------------------------------------------------" << std::endl;
  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  nSup2.setData(text, " ", 0, -1, false, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  mol = nSup2.next();
  delete mol;

  mol = nSup2[3];
  TEST_ASSERT(nSup2.length() == 4);
  delete mol;

  failed = false;
  try {
    mol = nSup2[4];
    delete mol;
  } catch (FileParseException &) {
    failed = true;
  }
  TEST_ASSERT(failed);

  BOOST_LOG(rdErrorLog)
      << "------------------------------------------------------" << std::endl;
  // this was a delightful boundary condition:
  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  nSup2.setData(text, " ", 0, -1, false, true);
  //  BOOST_LOG(rdErrorLog) << "SIZE: " << nSup2.length() << std::endl;
  failed = false;
  try {
    mol = nSup2[4];
    delete mol;
  } catch (FileParseException &) {
    failed = true;
  }
  TEST_ASSERT(failed);
  BOOST_LOG(rdErrorLog) << ">>> This may result in an infinite loop.  It "
                           "should finish almost immediately:"
                        << std::endl;
  TEST_ASSERT(nSup2.length() == 4);
  BOOST_LOG(rdErrorLog) << "<<< done." << std::endl;

  nSup2.reset();
  unsigned int nDone = 0;
  while (!nSup2.atEnd()) {
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  TEST_ASSERT(nDone == nSup2.length());

  // ensure that we can call setData a second time:
  text =
      "Id SMILES Column_2\n"
      "# comment, ignore\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-3 CCC 9.0\n"
      "mol-4 CCCC 16.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);
  mol = nSup2[2];
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-3");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "9.0");
  delete mol;

  mol = nSup2[1];
  mol->getProp(common_properties::_Name, mname);
  TEST_ASSERT(mname == "mol-2");
  mol->getProp("Column_2", mname);
  TEST_ASSERT(mname == "4.0");
  delete mol;

  // this was a delightful boundary condition:
  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC\n"
      "\n"
      "\n";
  nSup2.setData(text, " ", 0, -1, false, true);
  TEST_ASSERT(nSup2.length() == 4);
  nSup2.reset();
  nDone = 0;
  while (!nSup2.atEnd()) {
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  TEST_ASSERT(nDone == nSup2.length());
};

void testSmilesWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  // std::string fname = "../test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles_molsupplier.csv";
  // std::string oname = "../test_data/outSmiles.csv";

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();
  // BOOST_LOG(rdErrorLog) << "WRITING" << std::endl;
  while (mol) {
    // BOOST_LOG(rdErrorLog) << "MOL: " << MolToSmiles(*mol) << std::endl;
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    names.push_back(mname);
    props.push_back(pval);
    writer->write(*mol);
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  writer->flush();
  delete nSup;

  // now read the molecules back in a check if we have the same properties etc
  nSup = new SmilesMolSupplier(oname);
  int i = 0;
  mol = nSup->next();
  while (mol) {
    std::string mname, pval;
    mol->getProp(common_properties::_Name, mname);
    mol->getProp("Column_2", pval);
    CHECK_INVARIANT(mname == names[i], "");
    CHECK_INVARIANT(pval == props[i], "");
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  TEST_ASSERT(nSup->length() == writer->numMols());
  writer->close();
  delete writer;
  delete nSup;
}

void testSDWriter() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);

  std::string ofile =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/outNCI_few_molsupplier.sdf";

  auto *writer = new SDWriter(ofile);

  STR_VECT names;

  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    names.push_back(mname);

    writer->write(*mol);
    delete mol;
  }
  writer->flush();
  CHECK_INVARIANT(writer->numMols() == 16, "");
  writer->close();
  delete writer;

  // now read in the file we just finished writing

  SDMolSupplier reader(ofile);
  int i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    BOOST_LOG(rdInfoLog) << mname << "\n";
    // CHECK_INVARIANT(mname == names[i], "");

    delete mol;
    i++;
  }

  BOOST_LOG(rdInfoLog) << i << "\n";
  /*
  // now read in a file with aromatic information on the bonds
  std::string infile = rdbase +
  "/Code/GraphMol/FileParsers/test_data/outNCI_arom.sdf";
  SDMolSupplier nreader(infile);
  i = 0;
  while (!nreader.atEnd()) {
    ROMol *mol = nreader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    BOOST_LOG(rdInfoLog) << mname << "\n";
    //CHECK_INVARIANT(mname == names[i], "");
    i++;

    delete mol;
    }*/
}

void testSDSupplierEnding() {
  std::string rdbase = getenv("RDBASE");
  // test the SD supplier to check if it properly handle the end of sd file
  // conditions
  // should work fine if the sd file end with  a $$$$ follwed by blank line or
  // no
  // no blank lines
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/esters_end.sdf";
  int i = 0;
  SDMolSupplier reader(infile);
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    i++;
    delete mol;
  }
  CHECK_INVARIANT(i == 6, "");
}

void testSuppliersEmptyFile() {
  std::string rdbase = getenv("RDBASE");
  {  // contains no records
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty.sdf";
    SDMolSupplier reader(infile);
    TEST_ASSERT(reader.atEnd());
  }
  {
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty.smi";
    SmilesMolSupplier smiSup(infile, ",", 0, -1);
    TEST_ASSERT(smiSup.atEnd());
  }
  // tests for GitHub issue 19:
  {  // actually an empty file, throws an exception:
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty2.sdf";
    bool failed = false;
    try {
      SDMolSupplier reader(infile);
    } catch (BadFileException &) {
      failed = true;
    }
    TEST_ASSERT(failed);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    TEST_ASSERT(reader.atEnd());
    bool failed = false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed = true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length() == 0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    bool failed = false;
    try {
      reader[0];
    } catch (FileParseException &) {
      failed = true;
    }
    TEST_ASSERT(failed);
    TEST_ASSERT(reader.length() == 0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    TEST_ASSERT(reader.length() == 0);
  }
}

void testCisTrans() {
  std::string text;
  text =
      "mol-1 ClC(C)=C(Br)C\n"
      "mol-2 C1=COC=CC1C(Cl)=C(Br)C\n"
      "mol-3 C1=COC=CC1\\C(Cl)=C(Br)\\C";
  SmilesMolSupplier smiSup;
  smiSup.setData(text, " ", 1, 0, false, true);

  std::string ofile = "cisTrans_molsupplier.sdf";
  SDWriter writer(ofile);
  while (!smiSup.atEnd()) {
    ROMol *mol = smiSup.next();
    TEST_ASSERT(mol);
    RDDepict::compute2DCoords(*mol);
    writer.write(*mol);
    delete mol;
  }
  writer.close();
  // do the round t;est
  // parse the sd file and write it out to smiles

  SDMolSupplier *reader;
  try {
    reader = new SDMolSupplier("cisTrans_molsupplier.sdf");
  } catch (FileParseException &) {
    reader = nullptr;
  }
  TEST_ASSERT(reader);
  while (!reader->atEnd()) {
    ROMol *mol = reader->next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    BOOST_LOG(rdInfoLog) << mname << " ";
    BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, 1) << "\n";
    delete mol;
  }
  delete reader;
}

void testStereoRound() {
  // - we will read ina bunch of cdk2 smiles with stereo on them
  // - generate the canonical smiles for each one
  // - generate 2D coordinates, write to an sdf file
  // - read the sdf file back in and compare the canonical smiles
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/cdk2_stereo.csv";
  SmilesMolSupplier *smiSup;
  try {
    smiSup = new SmilesMolSupplier(infile, ",", 0, 1, false, true);
  } catch (FileParseException &) {
    smiSup = nullptr;
  }
  TEST_ASSERT(smiSup)
  std::map<std::string, std::string> nameSmi;
  std::string ofile =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/cdk2_stereo_molsupplier.sdf";
  auto *writer = new SDWriter(ofile);
  int count = 0;

  while (!smiSup->atEnd()) {
    ROMol *mol = smiSup->next();
    // mol->debugMol(std::cout);
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    nameSmi[mname] = MolToSmiles(*mol, 1);

    ROMol *nmol = SmilesToMol(nameSmi[mname]);
    // nmol->debugMol(std::cout);

    std::string nsmi = MolToSmiles(*nmol, 1);
    // BOOST_LOG(rdErrorLog) << mname << "\n";
    if (nameSmi[mname] != nsmi) {
      BOOST_LOG(rdInfoLog) << mname << " " << nameSmi[mname] << " " << nsmi
                           << "\n";
    }
    RDDepict::compute2DCoords(*mol);
    writer->write(*mol);
    count++;
    delete mol;
    delete nmol;

    if (count % 50 == 0) {
      BOOST_LOG(rdInfoLog) << count << " " << mname << "\n";
    }
  }
  writer->close();
  delete smiSup;
  delete writer;

  // now read the SD file back in check if the canonical smiles are the same
  SDMolSupplier *reader;
  try {
    reader = new SDMolSupplier(ofile);
  } catch (FileParseException &) {
    reader = nullptr;
  }
  TEST_ASSERT(reader);
  count = 0;

  while (!reader->atEnd()) {
    ROMol *mol = reader->next();
    // mol->debugMol(std::cout);
    std::string smiles = MolToSmiles(*mol, 1);
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    if (nameSmi[mname] != smiles) {
      BOOST_LOG(rdInfoLog) << mname << " " << nameSmi[mname] << " " << smiles
                           << "\n";
    }
    delete mol;
    count++;
  }
  delete reader;
}

void testIssue226() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue226.sdf";
  SDMolSupplier sdsup(fname);

  ROMol *mol;

  mol = sdsup.next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("E1"));
  TEST_ASSERT(mol->hasProp("E2"));
  delete mol;

  mol = sdsup.next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("E1"));
  TEST_ASSERT(mol->hasProp("E2"));
  delete mol;
}

int testTDTSupplier1() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  {
    TDTMolSupplier suppl(fname, "PN");
    unsigned int i = 0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1, prop2;
        TEST_ASSERT(nmol->getNumAtoms() > 0);
        TEST_ASSERT(nmol->hasProp("PN"));
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("MFCD"));

        nmol->getProp("PN", prop1);
        nmol->getProp(common_properties::_Name, prop2);
        TEST_ASSERT(prop1 == prop2);

        // we didn't ask for 2D conformers, so there should be a property 2D:
        TEST_ASSERT(nmol->hasProp(common_properties::TWOD));
        // and no conformer:
        TEST_ASSERT(!nmol->getNumConformers());

        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i == 10);
  }
  {
    std::ifstream strm(fname.c_str(), std::ios_base::binary);
    TDTMolSupplier suppl(&strm, false, "PN");
    unsigned int i = 0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1, prop2;
        TEST_ASSERT(nmol->getNumAtoms() > 0);
        TEST_ASSERT(nmol->hasProp("PN"));
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("MFCD"));

        nmol->getProp("PN", prop1);
        nmol->getProp(common_properties::_Name, prop2);
        TEST_ASSERT(prop1 == prop2);

        // we didn't ask for 2D conformers, so there should be a property 2D:
        TEST_ASSERT(nmol->hasProp(common_properties::TWOD));
        // and no conformer:
        TEST_ASSERT(!nmol->getNumConformers());

        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i == 10);
  }
  return 1;
}
int testTDTSupplier2() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  int i;
  std::string prop1, prop2;

  TDTMolSupplier suppl(fname, "PN", 2);
  i = 0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms() > 0);
      TEST_ASSERT(nmol->hasProp("PN"));
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));
      TEST_ASSERT(nmol->hasProp("MFCD"));

      nmol->getProp("PN", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      TEST_ASSERT(prop1 == prop2);

      // we asked for 2D conformers, so there should be no property 2D:
      TEST_ASSERT(!nmol->hasProp(common_properties::TWOD));
      // and a conformer:
      TEST_ASSERT(nmol->getNumConformers() == 1);
      // with id "2":
      TEST_ASSERT(nmol->beginConformers()->get()->getId() == 2);

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 10);
  return 1;
}
int testTDTSupplier3() {
  std::string data;
  int i;
  std::string prop1, prop2;

  TDTMolSupplier suppl;

  data =
      "$SMI<Cc1nnc(N)nc1C>\n"
      "CAS<17584-12-2>\n"
      "|\n"
      "$SMI<Cc1n[nH]c(=O)nc1N>\n"
      "CAS<~>\n"
      "|\n"
      "$SMI<Cc1n[nH]c(=O)[nH]c1=O>\n"
      "CAS<932-53-6>\n"
      "|\n"
      "$SMI<Cc1nnc(NN)nc1O>\n"
      "CAS<~>\n"
      "|\n";
  suppl.setData(data, "CAS");

  i = 0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms() > 0);
      TEST_ASSERT(nmol->hasProp("CAS"));
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));

      nmol->getProp("CAS", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      TEST_ASSERT(prop1 == prop2);

      // no conformers should have been read:
      TEST_ASSERT(nmol->getNumConformers() == 0);

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 4);
  TEST_ASSERT(suppl.length() == 4);

  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  ROMol *mol = suppl[0];
  TEST_ASSERT(mol);
  delete mol;

  // make sure we can reset the supplier and still process it properly;
  suppl.setData(data, "CAS");

  i = 0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      TEST_ASSERT(nmol->getNumAtoms() > 0);
      TEST_ASSERT(nmol->hasProp("CAS"));
      TEST_ASSERT(nmol->hasProp(common_properties::_Name));

      nmol->getProp("CAS", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      TEST_ASSERT(prop1 == prop2);

      // no conformers should have been read:
      TEST_ASSERT(nmol->getNumConformers() == 0);

      delete nmol;
      i++;
    }
  }
  TEST_ASSERT(i == 4);

  return 1;
}

void testSDSupplierFromText() {
  std::string text;
  int i = 0;
  SDMolSupplier reader;

  text =
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "M  END\n"
      ">  <ID> (3)\n"
      "Lig1\n"
      "\n"
      "$$$$\n"
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.8477    0.6988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "  3  6  1  0  0  0  0\n"
      "M  END\n"
      ">  <ID> (4)\n"
      "Lig2\n"
      "\n"
      "$$$$\n";
  reader.setData(text);

  i = 0;
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    std::string mname;
    TEST_ASSERT(mol->hasProp(common_properties::_Name));
    TEST_ASSERT(mol->hasProp("ID"));
    i++;
    delete mol;
  }
  TEST_ASSERT(i == 2);
}

void testSDSupplierFromTextStrLax1() {
  std::string text;
  text =
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "M  END\n"
      "blah\n"
      "\n"
      "blah after blank line\n"
      ">  <ID> (3)\n"
      "Lig1\n"
      "\n"
      "This will be ignored\n"
      ">  <ANOTHER_PROPERTY> (4)\n"
      "Value\n"
      "\n"
      "$$$$\n"
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.8477    0.6988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "  3  6  1  0  0  0  0\n"
      "M  END\n"
      ">  <ID> (4)\n"
      "Lig2\n"
      "\n"
      "This will be ignored\n"
      "\n"
      ">  <ANOTHER_PROPERTY> (4)\n"
      "Value\n"
      "\n"
      "This will be ignored\n"
      "\n"
      "$$$$\n";

  // strict
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      if (i == 0) {
        TEST_ASSERT(!mol->hasProp("ID"));
      }
      TEST_ASSERT(!mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 2);
  }
  // lax
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 2);
  }
}

void testSDSupplierFromTextStrLax2() {
  std::string text;
  text =
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "M  END\n"
      ">  <ID> (3)\n"
      "Lig1\n"
      "\n"
      ">  <ANOTHER_PROPERTY> (4)\n"
      "No blank line before dollars\n"
      "$$$$\n"
      "Structure1\n"
      "csChFnd70/05230312262D\n"
      "\n"
      "  6  5  0  0  0  0  0  0  0  0999 V2000\n"
      "    1.2124    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6373    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4249    2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0000    0.7000    0.0000 Y   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.8477    0.6988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  2  0  0  0  0\n"
      "  1  5  1  0  0  0  0\n"
      "  3  6  1  0  0  0  0\n"
      "M  END\n"
      ">  <ID> (3)\n"
      "Lig2\n"
      "\n"
      ">  <ANOTHER_PROPERTY> (4)\n"
      "Value2\n"
      "\n"
      "$$$$\n";

  // strict
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      TEST_ASSERT(s == "Lig1");
      mol->getProp("ANOTHER_PROPERTY", s);
      TEST_ASSERT(s ==
                  "No blank line before dollars\n"
                  "$$$$\n"
                  "Structure1\n"
                  "csChFnd70/05230312262D");
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 1);
  }
  // lax
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      TEST_ASSERT(s == "Lig2");
      mol->getProp("ANOTHER_PROPERTY", s);
      TEST_ASSERT(s == "Value2");
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 1);
  }
}

void testSDSupplierStrLax1() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/strictLax1.sdf";
  // strict
  {
    SDMolSupplier reader(fname, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      if (i == 0) {
        TEST_ASSERT(!mol->hasProp("ID"));
      }
      TEST_ASSERT(!mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 2);
  }
  // lax
  {
    SDMolSupplier reader(fname, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 2);
  }
}

void testSDSupplierStrLax2() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/strictLax2.sdf";
  // strict
  {
    SDMolSupplier reader(fname, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      TEST_ASSERT(s == "Lig1");
      mol->getProp("ANOTHER_PROPERTY", s);
      TEST_ASSERT(s ==
                  "No blank line before dollars\n"
                  "$$$$\n"
                  "Structure1\n"
                  "csChFnd70/05230312262D");
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 1);
  }
  // lax
  {
    SDMolSupplier reader(fname, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      TEST_ASSERT(mol->hasProp("ID"));
      TEST_ASSERT(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      TEST_ASSERT(s == "Lig2");
      mol->getProp("ANOTHER_PROPERTY", s);
      TEST_ASSERT(s == "Value2");
      i++;
      delete mol;
    }
    TEST_ASSERT(i == 1);
  }
}

void testIssue265() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NotThere.sdf";
  bool ok = false;
  try {
    SDMolSupplier reader(fname);
    ok = false;
  } catch (BadFileException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  try {
    SmilesMolSupplier reader(fname);
    ok = false;
  } catch (BadFileException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  try {
    TDTMolSupplier reader(fname);
    ok = false;
  } catch (BadFileException &) {
    ok = true;
  }
  TEST_ASSERT(ok);
}

void testSDErrorHandling() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors1.sdf";
  SDMolSupplier *sdsup;
  ROMol *nmol = nullptr;

  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(!nmol->hasProp("ID"));
  delete sdsup;
  delete nmol;

  // case 2: can't be sanitized
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors2.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;
  delete nmol;

  // entry 3: bad number of atoms
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors3.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;
  delete nmol;

  // entry 4: bad number of bonds
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors4.sdf";
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  nmol = sdsup->next();
  TEST_ASSERT(!nmol);
  TEST_ASSERT(sdsup->atEnd());
  delete sdsup;
  delete nmol;
}

void testIssue381() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue381.sdf";
  SDMolSupplier *sdsup;

  ROMol *nmol = nullptr;
  int count;

  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(!sdsup->atEnd());
  count = 0;
  while (!sdsup->atEnd()) {
    nmol = sdsup->next();
    if (nmol) {
      delete nmol;
    }
    count++;
  }
  TEST_ASSERT(sdsup->atEnd());
  TEST_ASSERT(count == 9);

  TEST_ASSERT(sdsup->length() == 9);

  delete sdsup;
}

void testSetStreamIndices() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  std::ifstream ifs(fname.c_str(), std::ios_base::binary);
  std::vector<std::streampos> indices;
  bool addIndex = true;
  bool notEof = true;
  std::streampos pos = 0;
  std::string line;
  while (notEof) {
    if (addIndex) {
      pos = ifs.tellg();
    }
    notEof = !std::getline(ifs, line).fail();
    if (notEof) {
      if (addIndex) {
        indices.push_back(pos);
      }
      addIndex = (line.substr(0, 4) == "$$$$");
    }
  }
  ifs.close();
  SDMolSupplier *sdsup;

  ROMol *nmol = nullptr;
  int count;

  sdsup = new SDMolSupplier(fname);
  sdsup->setStreamIndices(indices);
  TEST_ASSERT(!sdsup->atEnd());
  TEST_ASSERT(sdsup->length() == 16);

  count = 0;
  while (!sdsup->atEnd()) {
    nmol = sdsup->next();
    if (nmol) {
      delete nmol;
    }
    count++;
  }
  TEST_ASSERT(sdsup->atEnd());
  TEST_ASSERT(count == 16);

  TEST_ASSERT(sdsup->length() == 16);

  delete sdsup;
}

int testMixIterAndRandom() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/esters.sdf";
  bool ok;

  SDMolSupplier *sdsup;
  ROMol *mol;
  std::string name;

  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(sdsup);
  unsigned int i = 0;
  while (!sdsup->atEnd()) {
    mol = sdsup->next();
    if (mol) {
      TEST_ASSERT(mol->hasProp("ID"));
      delete mol;
    }
    i++;
  }
  TEST_ASSERT(i == 6);
  TEST_ASSERT(sdsup->length() == 6);

  delete sdsup;
  sdsup = new SDMolSupplier(fname);
  TEST_ASSERT(sdsup);
  TEST_ASSERT(sdsup->length() == 6);

  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID", name);
  TEST_ASSERT(name == "Lig1");
  delete mol;

  mol = (*sdsup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID", name);
  TEST_ASSERT(name == "Lig1");
  delete mol;

  sdsup->reset();
  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID", name);
  TEST_ASSERT(name == "Lig1");
  delete mol;
  mol = sdsup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->hasProp("ID"));
  mol->getProp("ID", name);
  TEST_ASSERT(name == "Lig2");
  delete mol;
  delete sdsup;

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  TEST_ASSERT(nSup->length() == 10);
  mol = (*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol = (*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(nSup->length() == 10);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol = nSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(nSup->length() == 10);
  delete mol;
  mol = (*nSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(nSup->length() == 10);
  delete mol;
  mol = nSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 20);
  delete nSup;
  delete mol;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  TEST_ASSERT(nSup);
  mol = nullptr;
  try {
    mol = (*nSup)[20];
    ok = false;
  } catch (FileParseException &) {
    ok = true;
  }
  TEST_ASSERT(ok);
  delete nSup;

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  TDTMolSupplier *tSup;
  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  TEST_ASSERT(tSup->length() == 10);
  mol = (*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol = (*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(tSup->length() == 10);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol = tSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(tSup->length() == 10);
  delete mol;

  mol = (*tSup)[0];
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(tSup->length() == 10);
  delete mol;

  mol = tSup->next();
  TEST_ASSERT(mol);
  delete mol;

  mol = tSup->next();
  TEST_ASSERT(mol);
  delete mol;

  mol = tSup->next();
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 10);
  delete tSup;
  delete mol;

  tSup = new TDTMolSupplier(fname);
  TEST_ASSERT(tSup);
  mol = nullptr;
  try {
    mol = (*tSup)[20];
    delete mol;
    ok = false;
  } catch (FileParseException &) {
    ok = true;
  }
  TEST_ASSERT(ok);
  delete tSup;

  return 1;
}

int testRemoveHs() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/withHs.sdf";

  SDMolSupplier sdsup(fname);
  ROMol *nmol;

  nmol = sdsup.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 23);
  delete nmol;
  nmol = sdsup.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 28);
  delete nmol;

  std::cerr << "build:" << std::endl;
  SDMolSupplier sdsup2(fname, true, false);
  nmol = sdsup2.next();
  TEST_ASSERT(nmol);
  // std::cerr<<" count: "<<nmol->getNumAtoms()<<std::endl;
  TEST_ASSERT(nmol->getNumAtoms() == 39);
  delete nmol;
  nmol = sdsup2.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 30);
  delete nmol;

  return 1;
}

void testGetItemText() {
  std::string rdbase = getenv("RDBASE");
  std::string fname;

  ROMol *mol1, *mol2;
  std::string molB, smiles;
  bool ok;

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);
    TEST_ASSERT(sdsup.length() == 16);

    molB = sdsup.getItemText(0);
    mol1 = sdsup[0];
    TEST_ASSERT(mol1);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    // make sure getItemText() doesn't screw up the current position:
    molB = sdsup.getItemText(10);
    mol1 = sdsup.next();
    molB = sdsup.getItemText(1);
    TEST_ASSERT(mol1);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = sdsup.getItemText(15);
    mol1 = sdsup[15];
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    try {
      molB = sdsup.getItemText(16);
      ok = false;
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);

    try {
      molB = sdsup.getItemText(20);
      ok = false;
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);

    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = sdsup.getItemText(0);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    mol1 = sdsup[0];
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    molB = sdsup.getItemText(5);
    mol2 = MolBlockToMol(molB);
    TEST_ASSERT(mol2);
    TEST_ASSERT(mol2->getNumAtoms() == 16);
    mol1 = sdsup[5];
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);
    TEST_ASSERT(smisup.length() == 10);

    molB = smisup.getItemText(0);
    TEST_ASSERT(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");
    mol1 = smisup[0];
    TEST_ASSERT(mol1);
    delete mol1;

    molB = smisup.getItemText(5);
    TEST_ASSERT(
        molB ==
        "6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");
    mol1 = smisup.next();
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol1->getNumAtoms() == 20);
    delete mol1;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = smisup.getItemText(8);
    TEST_ASSERT(molB == "9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");

    mol1 = smisup[0];
    TEST_ASSERT(mol1);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles == "CC1=CC(=O)C=CC1=O");
    TEST_ASSERT(mol1->getNumAtoms() == 9);
    delete mol1;
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);

    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = smisup.getItemText(0);
    TEST_ASSERT(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");

    molB = smisup.getItemText(5);
    TEST_ASSERT(
        molB ==
        "6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");

    molB = smisup.getItemText(8);
    TEST_ASSERT(molB == "9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);

    // make sure getItemText() flags EOF
    // (this was sf.net issue 3299878)
    molB = smisup.getItemText(0);
    TEST_ASSERT(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");

    ROMol *m = smisup[9];
    TEST_ASSERT(m);
    delete m;
    TEST_ASSERT(smisup.atEnd());
    molB = smisup.getItemText(9);
    TEST_ASSERT(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
    TEST_ASSERT(smisup.atEnd());
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = tdtsup.getItemText(0);
    TEST_ASSERT(molB != "");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    TEST_ASSERT(tdtsup.length() == 10);

    molB = tdtsup.getItemText(0);
    TEST_ASSERT(molB != "");

    mol1 = tdtsup[0];
    TEST_ASSERT(mol1);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles == "Cc1nnc(N)nc1C");
    TEST_ASSERT(mol1->getNumAtoms() == 9);
    delete mol1;

    // make sure getItemText doesn't screw up next()
    molB = tdtsup.getItemText(5);
    mol1 = tdtsup.next();
    TEST_ASSERT(mol1);
    TEST_ASSERT(mol1->getNumAtoms() == 9);
    smiles = MolToSmiles(*mol1, 1);
    TEST_ASSERT(smiles == "Cc1n[nH]c(=O)nc1N");
    delete mol1;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = tdtsup.getItemText(9);
    TEST_ASSERT(molB != "");
    TEST_ASSERT(molB.substr(0, 12) == "$SMI<Cc1n[nH");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    TEST_ASSERT(tdtsup.length() == 10);

    ROMol *mol = tdtsup[9];
    TEST_ASSERT(mol);
    delete mol;
    TEST_ASSERT(tdtsup.atEnd());

    // (this was sf.net issue 3299878
    molB = tdtsup.getItemText(9);
    TEST_ASSERT(molB != "");
    TEST_ASSERT(molB.substr(0, 12) == "$SMI<Cc1n[nH");
    TEST_ASSERT(tdtsup.atEnd());
  }
}

int testForwardSDSupplier() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  std::string fname2 =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf.gz";

  {
    std::ifstream strm(fname.c_str());
    ForwardSDMolSupplier sdsup(&strm, false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      TEST_ASSERT(nmol || sdsup.atEnd());
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i == 16);
  }
#ifdef RDK_USE_BOOST_IOSTREAMS
  // make sure the boost::iostreams are working
  {
    io::filtering_istream strm;
    strm.push(io::file_source(fname));

    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm, line);
      if (!strm.eof()) {
        ++i;
      }
      if (i > 1000) {
        break;
      }
    }
    TEST_ASSERT(i == 998);
  }
  {
    gzstream strm(fname2);
    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm, line);
      if (!strm.eof()) {
        ++i;
      }
      if (i > 1000) {
        break;
      }
    }
    TEST_ASSERT(i == 997);
  }
  // looks good, now do a supplier:
  {
    gzstream strm(fname2);

    ForwardSDMolSupplier sdsup(&strm, false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        TEST_ASSERT(nmol->hasProp(common_properties::_Name));
        TEST_ASSERT(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    TEST_ASSERT(i == 16);
  }
#endif

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  // Now test that Maestro parsing of gz files works
  std::string maefname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.mae";
  std::string maefname2 =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.maegz";
  {
    io::filtering_istream strm;
    strm.push(io::file_source(maefname));

    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm, line);
      if (!strm.eof()) {
        ++i;
      }
      if (i > 1700) {
        break;
      }
    }
    TEST_ASSERT(i == 1663);
  }
#if RDK_USE_BOOST_IOSTREAMS  
  {
    gzstream strm(maefname2);

    unsigned int i = 0;
    while (!strm.eof()) {
      std::string line;
      std::getline(strm, line);
      if (!strm.eof()) {
        ++i;
      }
      if (i > 1700) {
        break;
      }
    }
    TEST_ASSERT(i == 1663);
  }
  // looks good, now do a supplier:
  {
    auto *strm = new gzstream(maefname2);

    MaeMolSupplier maesup(strm);
    unsigned int i = 0;
    std::shared_ptr<ROMol> nmol;
    while (!maesup.atEnd()) {
      nmol.reset(maesup.next());
      if (nmol != nullptr) {
        i++;
      }
    }
    TEST_ASSERT(i == 16);
  }
#endif
#endif  // RDK_BUILD_MAEPARSER_SUPPORT

  return 1;
}

void testMissingCRSDSupplier() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/missingCR.sdf";
  SDMolSupplier reader(infile);
  auto *mol = reader.next();
  delete mol;
  TEST_ASSERT(reader.atEnd());
}

void testIssue3482695() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3482695.sdf";
  SDMolSupplier reader(infile);
  ROMol *nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 0);
  TEST_ASSERT(nmol->hasProp("test"));
  delete nmol;
}

void testIssue3525673() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525673.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 37);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 58);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(!nmol);  // broken due to 'foo' in counts line!

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 58);
  delete nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  delete nmol;
}

void testBlankLinesInProps() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/BlankPropLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 19);
  TEST_ASSERT(nmol->hasProp("MultiLineProperty1"));
  nmol->getProp("MultiLineProperty1", pval);
  TEST_ASSERT(pval == "foo\nbar\n \nbaz");
  TEST_ASSERT(nmol->hasProp("MultiLineProperty2"));
  TEST_ASSERT(!(nmol->hasProp("fooprop")));
  nmol->getProp("MultiLineProperty2", pval);
  TEST_ASSERT(pval == "foo\n>  <fooprop>\nbaz\n ");
  delete nmol;
}

void testSkipLines() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/SkipLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 1);
  TEST_ASSERT(nmol->hasProp("prop1"));
  delete nmol;
}

void testGitHub23() {
  std::string rdbase = getenv("RDBASE");
  std::string ofile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/blah_molsupplier.sdf";
  auto *writer = new SDWriter(ofile);

  ROMol *mol = SmilesToMol("CCCC");
  INT_VECT iv;
  iv.push_back(1);
  iv.push_back(2);
  mol->setProp("pval", iv);
  writer->write(*mol);
  delete mol;

  writer->close();
  delete writer;
}

void testGitHub88() {
  std::string rdbase = getenv("RDBASE");
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/github88.v3k.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;

  nmol = reader.next();
  TEST_ASSERT(nmol);
  TEST_ASSERT(nmol->getNumAtoms() == 8);
  TEST_ASSERT(nmol->hasProp("prop1"));
  std::string pval;
  nmol->getProp("prop1", pval);
  TEST_ASSERT(pval == "4");
  delete nmol;
}

void testGitHub2285() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/github2285.sdf";

  std::vector<std::string> smiles;
  {
    SDMolSupplier sdsup(fname);
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      TEST_ASSERT(nmol);
      smiles.push_back(MolToSmiles(*nmol));
      delete nmol;
    }
  }
  {
    SDMolSupplier sdsup(fname, true, false);
    int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      TEST_ASSERT(nmol);
      ROMol *m = MolOps::removeHs(*nmol);
      TEST_ASSERT(MolToSmiles(*m) == smiles[i++]);
      delete nmol;
      delete m;
    }
    TEST_ASSERT(i > 0);
  }
}

void testGitHub2479() {
  std::string smiles1 = R"DATA(smiles id
c1ccccc duff
c1ccccc1 ok
C(C garbage
C1CC1 ok2
CC(C)(C)(C)C duff2
)DATA";
  {
    SmilesMolSupplier suppl;
    suppl.setData(smiles1);
    unsigned int cnt = 0;
    while (!suppl.atEnd()) {
      std::unique_ptr<ROMol> mol(suppl.next());
      if (cnt % 2) {
        TEST_ASSERT(mol);
      }
      ++cnt;
    }
    TEST_ASSERT(cnt == 5);
  }

  std::string sdf1 = R"SDF(
  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -13.3985    4.9850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.7066    5.4343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0654    4.9151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$

  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6408    5.3345    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0277    4.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
)SDF";
  {
    std::stringstream iss(sdf1);
    SDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(!mol2);
    TEST_ASSERT(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf1);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(!mol2);
    TEST_ASSERT(!suppl.atEnd());
    TEST_ASSERT(!suppl.getEOFHitOnRead());
    std::unique_ptr<ROMol> mol3(suppl.next());
    TEST_ASSERT(!mol3);
    TEST_ASSERT(suppl.atEnd());
    TEST_ASSERT(suppl.getEOFHitOnRead());
  }

  // truncated file1
  std::string sdf2 = R"SDF(
  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -13.3985    4.9850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.7066    5.4343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0654    4.9151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$

  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6408    5.3345    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0277    4.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$

  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6
)SDF";
  {
    std::stringstream iss(sdf2);
    SDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(!mol2);
    std::unique_ptr<ROMol> mol3(suppl.next());
    TEST_ASSERT(!mol3);
    TEST_ASSERT(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf2);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(!mol2);
    TEST_ASSERT(!suppl.atEnd());
    TEST_ASSERT(!suppl.getEOFHitOnRead());
    std::unique_ptr<ROMol> mol3(suppl.next());
    TEST_ASSERT(!mol3);
    TEST_ASSERT(suppl.atEnd());
    TEST_ASSERT(!suppl.getEOFHitOnRead());
  }
  // truncated file2
  std::string sdf3 = R"SDF(
  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -13.3985    4.9850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.7066    5.4343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0654    4.9151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <pval>  (1) 
[1,2,]

$$$$

  Mrv1810 06051911332D          

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6408    5.3345    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0277    4.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <pval>  (1) 
[1,2,]
)SDF";
  {
    std::stringstream iss(sdf3);
    SDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(mol2);
    TEST_ASSERT(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf3);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    TEST_ASSERT(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    TEST_ASSERT(mol2);
    TEST_ASSERT(suppl.atEnd());
  }
}

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
void testGitHub2881() {
  std::string data = R"DATA(f_m_ct { 
 s_m_title
 s_m_entry_id
 s_m_entry_name
 s_m_Source_Path
 s_m_Source_File
 i_m_Source_File_Index
 s_st_Chirality_1
 s_st_Chirality_2
 s_m_subgroup_title
 s_m_subgroupid
 b_m_subgroup_collapsed
 i_m_ct_format
 :::
 "Untitled Document-4" 
  17 
  newTemplates2.1 
  /Users/nicola/schrodinger/coordgen_standalone 
  templates.mae 
  17
  3_S_4_6_2 
  7_S_8_9_6_10 
  templates->templates->templates 
  templates->templates1->templates11 
  0
  2
 m_depend[2] { 
  # First column is dependency index #
  i_m_depend_dependency
  s_m_depend_property
  :::
  1 10 s_st_Chirality_1 
  2 10 s_st_Chirality_2 
  :::
 } 
 m_atom[15] { 
  # First column is atom index #
  i_m_mmod_type
  r_m_x_coord
  r_m_y_coord
  r_m_z_coord
  i_m_residue_number
  i_m_color
  i_m_atomic_number
  s_m_color_rgb
  s_m_atom_name
  :::
  1 5 1.186400 1.035900 0.000000 900 2 6 A0A0A0  C1 
  2 5 0.370300 1.157000 0.000000 900 2 6 A0A0A0  C2 
  3 4 -0.326500 0.715300 0.000000 900 2 6 A0A0A0  C3 
  4 5 0.085100 0.000400 0.000000 900 2 6 A0A0A0  C4 
  5 26 -0.328300 -0.713600 0.000000 900 43 7 5757FF  N5 
  6 5 -1.151500 0.716400 0.000000 900 2 6 A0A0A0  C6 
  7 5 -1.564900 0.002400 0.000000 900 2 6 A0A0A0  C7 
  8 5 -1.153300 -0.712600 0.000000 900 2 6 A0A0A0  C9 
  9 2 1.724800 0.410800 0.000000 900 2 6 A0A0A0  C12 
  10 2 1.723800 -0.414200 0.000000 900 2 6 A0A0A0  C13 
  11 5 1.183800 -1.037900 0.000000 900 2 6 A0A0A0  C14 
  12 5 0.367400 -1.157000 0.000000 900 2 6 A0A0A0  C15 
  13 7 2.508100 -0.670100 0.000000 900 2 6 A0A0A0  C16 
  14 7 2.993800 -0.003300 0.000000 900 2 6 A0A0A0  C17 
  15 29 2.509700 0.664800 0.000000 900 43 7 5757FF  N18 
  :::
 } 
 m_bond[17] { 
  # First column is bond index #
  i_m_from
  i_m_to
  i_m_order
  :::
  1 1 2 1
  2 1 9 1
  3 2 3 1
  4 3 4 1
  5 3 6 1
  6 4 5 1
  7 5 8 1
  8 5 12 1
  9 6 7 1
  10 7 8 1
  11 9 10 2
  12 9 15 1
  13 10 11 1
  14 10 13 1
  15 11 12 1
  16 13 14 2
  17 14 15 1
  :::
 } 
} 
)DATA";
  {
    auto *iss = new std::istringstream(data);
    bool sanitize = false;
    bool takeOwnership = true;
    MaeMolSupplier suppl(iss, takeOwnership, sanitize);
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (const Invar::Invariant &) {
    }
    TEST_ASSERT(!mol);
  }
}
#else
void testGitHub2881() {}
#endif

void testGitHub3517() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";

  SDMolSupplier sdsup(fname);
  TEST_ASSERT(!sdsup.atEnd());
  size_t l = sdsup.length();
  TEST_ASSERT(l > 0);
  TEST_ASSERT(!sdsup.atEnd());
}

int main() {
  RDLog::InitLogs();

#if 1
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testMolSup();
  BOOST_LOG(rdErrorLog) << "Finished: testMolSup()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testRandMolSup();
  BOOST_LOG(rdErrorLog) << "Finished: testRandMolSup()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesSup();
  BOOST_LOG(rdErrorLog) << "Finished: testSmilesSup()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesSupFromText();
  BOOST_LOG(rdErrorLog) << "Finished: testSmilesSupFromText()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSmilesWriter();
  BOOST_LOG(rdErrorLog) << "Finished: testSmilesWriter()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDWriter();
  BOOST_LOG(rdErrorLog) << "Finished: testSDWriter()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierEnding();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierEnding()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSuppliersEmptyFile();
  BOOST_LOG(rdErrorLog) << "Finished: testSuppliersEmptyFile()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testCisTrans();
  BOOST_LOG(rdErrorLog) << "Finished: testCisTrans()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testStereoRound();
  BOOST_LOG(rdErrorLog) << "Finished: testStereoRound()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue226();
  BOOST_LOG(rdErrorLog) << "Finished: testIssue226()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier1();
  BOOST_LOG(rdErrorLog) << "Finished: testTDTSupplier1()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier2();
  BOOST_LOG(rdErrorLog) << "Finished: testTDTSupplier2()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testTDTSupplier3();
  BOOST_LOG(rdErrorLog) << "Finished: testTDTSupplier3()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierFromText();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierFromText()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierStrLax1();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierStrLax1()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierStrLax2();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierStrLax2()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierFromTextStrLax1();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierFromTextStrLax1()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDSupplierFromTextStrLax2();
  BOOST_LOG(rdErrorLog) << "Finished: testSDSupplierFromTextStrLax2()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue265();
  BOOST_LOG(rdErrorLog) << "Finished: testIssue265()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSDErrorHandling();
  BOOST_LOG(rdErrorLog) << "Finished: testSDErrorHandling()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue381();
  BOOST_LOG(rdErrorLog) << "Finished: testIssue381()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSetStreamIndices();
  BOOST_LOG(rdErrorLog) << "Finished: testSetStreamIndices()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testMixIterAndRandom();
  BOOST_LOG(rdErrorLog) << "Finished: testMixIterAndRandom()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testRemoveHs();
  BOOST_LOG(rdErrorLog) << "Finished: testRemoveHs()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGetItemText();
  BOOST_LOG(rdErrorLog) << "Finished: testGetItemText()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testForwardSDSupplier();
  BOOST_LOG(rdErrorLog) << "Finished: testForwardSDSupplier()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testMissingCRSDSupplier();
  BOOST_LOG(rdErrorLog) << "Finished: testMissingCRSDSupplier()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue3482695();
  BOOST_LOG(rdErrorLog) << "Finished: testIssue3482695()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testIssue3525673();
  BOOST_LOG(rdErrorLog) << "Finished: testIssue3525673()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testBlankLinesInProps();
  BOOST_LOG(rdErrorLog) << "Finished: testBlankLinesInProps()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testSkipLines();
  BOOST_LOG(rdErrorLog) << "Finished: testSkipLines()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub23();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub23()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub88();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub88()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub2285();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub2285()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
#endif

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub2479();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub2479()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub2881();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub2881()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n";
  testGitHub3517();
  BOOST_LOG(rdErrorLog) << "Finished: testGitHub3517()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  return 0;
}
