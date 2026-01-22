//
//   Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <csignal>
#include <fstream>
#include <future>
#include <map>
#include <memory>
#include <string>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/RDKitBase.h>
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

namespace io = boost::iostreams;

static const std::string rdbase = getenv("RDBASE");

using namespace RDKit;

TEST_CASE("testMolSup") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";

  {
    SDMolSupplier sdsup(fname);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    REQUIRE(i == 16);
  }
  {
    SDMolSupplier sdsup(fname);
    for (unsigned int i = 0; i < 16; ++i) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
    }
    // test issue 3524949:
    REQUIRE(sdsup.atEnd());
    CHECK_THROWS_AS(sdsup.next(), FileParseException);
  }
  {
    std::ifstream strm(fname.c_str());
    SDMolSupplier sdsup(&strm, false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    REQUIRE(i == 16);
  }
  {
    auto *strm = new std::ifstream(fname.c_str());
    SDMolSupplier sdsup(strm, true);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
      }
      i++;
    }
    REQUIRE(i == 16);
  }
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  for (const bool useLegacy : {true, false}) {
    UseLegacyStereoPerceptionFixture fx(useLegacy);
    {  // Test reading properties
      fname = rdbase + "/Code/GraphMol/FileParsers/test_data/props_test.mae";

      MaeMolSupplier maesup(fname);
      std::unique_ptr<ROMol> nmol(maesup.next());
      REQUIRE(nmol);

      // Test mol properties
      REQUIRE(nmol->hasProp(common_properties::_Name));
      REQUIRE(nmol->hasProp("b_sd_chiral_flag"));
      REQUIRE(nmol->getProp<bool>("b_sd_chiral_flag") == false);
      REQUIRE(nmol->hasProp("i_sd_NSC"));
      REQUIRE(nmol->getProp<int>("i_sd_NSC") == 48);
      REQUIRE(nmol->hasProp("s_m_entry_name"));
      REQUIRE(nmol->getProp<std::string>("s_m_entry_name") == "NCI_aids_few.1");
      REQUIRE(nmol->hasProp("r_f3d_dummy"));
      REQUIRE(std::abs(nmol->getProp<double>("r_f3d_dummy") - 42.123) < 0.0001);

      // Test atom properties
      REQUIRE(nmol->getNumAtoms() == 19);
      for (int i = 0; i < 19; ++i) {
        const auto *atom = nmol->getAtomWithIdx(i);

        // The integer property is present for all atoms
        REQUIRE(atom->hasProp("i_m_minimize_atom_index"));
        REQUIRE(atom->getProp<int>("i_m_minimize_atom_index") == 1 + i);

        // The bool property is only defined for i < 10
        if (i < 10) {
          REQUIRE(atom->hasProp("b_m_dummy"));
          REQUIRE(atom->getProp<bool>("b_m_dummy") == static_cast<bool>(i % 2));
        } else {
          REQUIRE(!atom->hasProp("b_m_dummy"));
        }

        // The real property is only defined for i >= 10
        if (i >= 10) {
          REQUIRE(atom->hasProp("r_f3d_dummy"));
          REQUIRE(std::abs(atom->getProp<double>("r_f3d_dummy") - (19.1 - i)) <
                  0.0001);
        } else {
          REQUIRE(!atom->hasProp("r_f3d_dummy"));
        }

        // All atoms have the string prop
        REQUIRE(atom->hasProp("s_m_dummy"));
        REQUIRE(atom->getProp<std::string>("s_m_dummy") ==
                std::to_string(19 - i));
      }

      REQUIRE(maesup.atEnd());
    }

    {  // Test parsing stereo properties. Mol is 2D and has stereo labels.
      fname = rdbase + "/Code/GraphMol/FileParsers/test_data/stereochem.mae";
      MaeMolSupplier maesup(fname);

      {  // Stereo bonds. These get overwritten by the double bond detection.
        std::unique_ptr<ROMol> nmol(maesup.next());
        REQUIRE(nmol);
        {
          Bond *bnd = nmol->getBondWithIdx(1);
          REQUIRE(bnd);
          REQUIRE(bnd->getStereoAtoms() == INT_VECT({0, 3}));
          REQUIRE(bnd->getStereo() == Bond::STEREOTRANS);
        }
        {
          Bond *bnd = nmol->getBondWithIdx(3);
          REQUIRE(bnd);
          REQUIRE(bnd->getStereoAtoms() == INT_VECT({2, 5}));
          REQUIRE(bnd->getStereo() == Bond::STEREOCIS);
        }
      }
      {  // Chiralities (these get CIP codes)
        std::unique_ptr<ROMol> nmol(maesup.next());
        REQUIRE(nmol);
        {
          Atom *at = nmol->getAtomWithIdx(1);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
          if (useLegacy) {
            REQUIRE(at->getProp<std::string>(common_properties::_CIPCode) ==
                    "R");
          }
        }
        {
          Atom *at = nmol->getAtomWithIdx(3);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
          if (useLegacy) {
            REQUIRE(at->getProp<std::string>(common_properties::_CIPCode) ==
                    "S");
          }
        }
      }
      {  // Pseudochiralities (no CIP codes)
        std::unique_ptr<ROMol> nmol(maesup.next());
        REQUIRE(nmol);
        {
          Atom *at = nmol->getAtomWithIdx(2);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
          REQUIRE(!at->hasProp(common_properties::_CIPCode));
        }
        {
          Atom *at = nmol->getAtomWithIdx(5);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
          REQUIRE(!at->hasProp(common_properties::_CIPCode));
        }
      }
      {  // intentionally bad chirality label, intended to
        // make sure we can step over parse errors
        std::unique_ptr<ROMol> nmol(maesup.next());
        REQUIRE(nmol);
        REQUIRE(nmol->getNumAtoms() > 1);
      }
      {  // "Undefined" chirality label
        std::unique_ptr<ROMol> nmol(maesup.next());
        REQUIRE(nmol);
        {
          Atom *at = nmol->getAtomWithIdx(2);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_UNSPECIFIED);
          REQUIRE(!at->hasProp(common_properties::_CIPCode));
        }
        {
          Atom *at = nmol->getAtomWithIdx(5);
          REQUIRE(at);
          REQUIRE(at->getChiralTag() == Atom::CHI_UNSPECIFIED);
          REQUIRE(!at->hasProp(common_properties::_CIPCode));
        }
      }
      REQUIRE(maesup.atEnd());
    }
    {  // Test loop reading
      fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.mae";
      MaeMolSupplier maesup(fname);
      std::shared_ptr<ROMol> nmol;
      for (unsigned int i = 0; i < 16; ++i) {
        nmol.reset(maesup.next());
        if (nmol) {
          REQUIRE(nmol->hasProp(common_properties::_Name));
          REQUIRE(nmol->getNumAtoms() > 0);
          if (i == 0) {
            auto smiles = MolToSmiles(*nmol);
            REQUIRE(smiles ==
                    "CCC1=[O+][Cu@]2([O+]=C(CC)C1)[O+]=C(CC)CC(CC)=[O+]2");
          }
        }
      }
      REQUIRE(maesup.atEnd());
      CHECK_THROWS_AS(maesup.next(), FileParseException);
    }

    {
      fname = rdbase + "/Code/GraphMol/FileParsers/test_data/bad_ppty.mae";
      const std::string err_msg_substr = "Bad format for property";

      bool ok = false;
      std::unique_ptr<ROMol> mol;
      MaeMolSupplier maesup(fname);

      // This is in excess: there are only 3 mols in the file, and the second
      // one has an invalid property name, so it won't be read
      for (unsigned int i = 0; i < 5; ++i) {
        try {
          mol.reset(maesup.next());
        } catch (const FileParseException &e) {
          const std::string err_msg(e.what());
          REQUIRE(i == 1);
          REQUIRE(err_msg.find(err_msg_substr) != std::string::npos);
          ok = true;
          break;
        }
        REQUIRE(mol);
        REQUIRE(mol->hasProp(common_properties::_Name));
        REQUIRE(mol->getNumAtoms() == 1);
        REQUIRE(!maesup.atEnd());
      }
      REQUIRE(!maesup.atEnd());
      REQUIRE(ok);
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
      REQUIRE(info->getResidueName() == "ARG ");
      REQUIRE(info->getChainId() == "A");
      REQUIRE(info->getResidueNumber() == 5);
    }

#endif
  }
#endif  // RDK_BUILD_MAEPARSER_SUPPORT
}

TEST_CASE("testRandMolSup") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  SDMolSupplier sdsup(fname);

  ROMol *tmol = sdsup[7];
  delete tmol;

  REQUIRE(sdsup.length() == 16);

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
    REQUIRE(mname == names[i]);
    delete mol;
  }

  // get a random molecule
  ROMol *mol = sdsup[5];
  REQUIRE(mol);
  std::string mname;
  mol->getProp(common_properties::_Name, mname);
  delete mol;
  REQUIRE(mname == "170");

  // get the last molecule:
  mol = sdsup[15];
  REQUIRE(mol);
  delete mol;

  // and make sure we're at the end:
  REQUIRE(sdsup.atEnd());
  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  mol = sdsup[0];
  REQUIRE(mol);
  delete mol;

  // Issue 113: calling length before grabbing a molecule results in crashes:
  SDMolSupplier sdsup2(fname);
  REQUIRE(sdsup2.length() == 16);
}

TEST_CASE("testSmilesSup") {
  std::string mname;
  std::string fname;
  ROMol *mol;

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);
    REQUIRE(nSup2.length() == 10);
  }
  {
    SmilesMolSupplier nSup2(fname, ",", 1, 0, true);

    mol = nSup2[3];
    REQUIRE(!nSup2.atEnd());
    REQUIRE(nSup2.length() == 10);

    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "4");
    mol->getProp("TPSA", mname);
    REQUIRE(mname == "82.78");
    delete mol;

    mol = nSup2[9];
    REQUIRE(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    REQUIRE(mol);
    delete mol;
  }
  {
    std::ifstream strm(fname.c_str(), std::ios_base::binary);
    SmilesMolSupplier nSup2(&strm, false, ",", 1, 0, true);

    mol = nSup2[3];
    REQUIRE(nSup2.length() == 10);

    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "4");
    mol->getProp("TPSA", mname);
    REQUIRE(mname == "82.78");
    delete mol;

    mol = nSup2[9];
    REQUIRE(mol);
    delete mol;
    // now make sure we can grab earlier mols (was sf.net issue 1904170):
    mol = nSup2[0];
    REQUIRE(mol);
    delete mol;
  }

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  SmilesMolSupplier smiSup(fname, ",", 0, -1);

  mol = smiSup[16];

  mol->getProp("TPSA", mname);
  REQUIRE(mname == "46.25");
  delete mol;

  mol = smiSup[8];
  mol->getProp("TPSA", mname);
  REQUIRE(mname == "65.18");
  delete mol;

  int len = smiSup.length();
  REQUIRE(len == 200);

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

  REQUIRE(i == 200);

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);

  // check the length before we read anything out...
  //  this was a problem at one point (Issue 113)
  REQUIRE(nSup->length() == 10);
  mol = (*nSup)[3];

  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "82.78");

  delete nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  unsigned int nRead = 0;
  while (!nSup->atEnd()) {
    delete mol;
    mol = nSup->next();
    REQUIRE(mol);
    nRead++;
  }
  REQUIRE(nSup->length() == 10);
  REQUIRE(nRead == 10);

  delete nSup;
  delete mol;
}

TEST_CASE("testSmilesSupFromText") {
  std::string mname;
  std::string fname;
  ROMol *mol;

  SmilesMolSupplier nSup2;
  std::string text;
  int nAts;

  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  {
    nSup2.setData(text, " ", 0, -1, false, true);

    mol = nSup2.next();
    nAts = mol->getNumAtoms();
    delete mol;
    REQUIRE(nAts == 2);

    mol = nSup2[3];
    nAts = mol->getNumAtoms();
    delete mol;
    REQUIRE(nAts == 6);
    REQUIRE(nSup2.length() == 4);

    CHECK_THROWS_AS(nSup2[4], FileParseException);

    mol = nSup2[2];
    nAts = mol->getNumAtoms();
    REQUIRE(nAts == 4);
    REQUIRE(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "2");
    delete mol;
  }
  {
    nSup2.setData(text, " ", 0, -1, false, true);
    mol = nSup2[2];
    REQUIRE(mol);
    nAts = mol->getNumAtoms();
    REQUIRE(nAts == 4);
    REQUIRE(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "2");
    delete mol;

    mol = nSup2[3];
    REQUIRE(mol);
    nAts = mol->getNumAtoms();
    REQUIRE(nAts == 6);
    REQUIRE(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "3");
    delete mol;
  }
  {
    nSup2.setData(text, " ", 0, -1, false, true);
    mol = nSup2[3];
    REQUIRE(mol);
    nAts = mol->getNumAtoms();
    REQUIRE(nAts == 6);
    REQUIRE(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "3");

    delete mol;
    mol = nSup2[2];
    REQUIRE(mol);
    nAts = mol->getNumAtoms();
    REQUIRE(nAts == 4);
    REQUIRE(mol->hasProp(common_properties::_Name));
    mol->getProp(common_properties::_Name, mname);
    REQUIRE(mname == "2");
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
  nSup2.setData(text, " ", 1, 0, true, true);
  mol = nSup2[3];

  REQUIRE(nSup2.length() == 4);
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "16.0");
  delete mol;

  // ensure that we can call setData a second time:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-3 CCC 9.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);
  REQUIRE(nSup2.length() == 3);
  mol = nSup2[2];
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-3");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "9.0");
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

  REQUIRE(nSup2.length() == 4);
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "16.0");
  delete mol;

  // failures should give null molecules:
  mol = nSup2[2];
  REQUIRE(!mol);
  delete mol;

  // issue 114, no \n at EOF:
  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0\n";
  nSup2.setData(text, " ", 1, 0, true, true);

  REQUIRE(nSup2.length() == 3);
  mol = nSup2[2];
  REQUIRE(mol);
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "16.0");
  REQUIRE(nSup2.atEnd());
  delete mol;

  text =
      "Id SMILES Column_2\n"
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0";
  nSup2.setData(text, " ", 1, 0, true, true);

  REQUIRE(nSup2.length() == 3);
  mol = nSup2[2];
  REQUIRE(mol);
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "16.0");
  REQUIRE(nSup2.atEnd());
  delete mol;

  CHECK_THROWS_AS(nSup2[3], FileParseException);

  text =
      "mol-1 C 1.0\n"
      "mol-2 CC 4.0\n"
      "mol-4 CCCC 16.0";
  nSup2.setData(text, " ", 1, 0, false, true);

  REQUIRE(nSup2.length() == 3);
  mol = nSup2[2];
  REQUIRE(mol);
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-4");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "16.0");
  delete mol;

  text =
      "C\n"
      "CC\n"
      "CCCC";
  nSup2.setData(text, " ", 0, -1, false, true);

  REQUIRE(nSup2.length() == 3);
  mol = nSup2[2];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 4);
  delete mol;

  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  nSup2.setData(text, " ", 0, -1, false, true);

  mol = nSup2.next();
  delete mol;

  mol = nSup2[3];
  REQUIRE(nSup2.length() == 4);
  delete mol;

  CHECK_THROWS_AS(nSup2[4], FileParseException);

  text =
      "CC\n"
      "CCC\n"
      "CCOC\n"
      "CCCCOC";
  nSup2.setData(text, " ", 0, -1, false, true);

  CHECK_THROWS_AS(nSup2[4], FileParseException);

  // This may result in an infinite loop.  It should finish almost immediately:
  REQUIRE(nSup2.length() == 4);

  nSup2.reset();
  unsigned int nDone = 0;
  while (!nSup2.atEnd()) {
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  REQUIRE(nDone == nSup2.length());

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
  REQUIRE(mname == "mol-3");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "9.0");
  delete mol;

  mol = nSup2[1];
  mol->getProp(common_properties::_Name, mname);
  REQUIRE(mname == "mol-2");
  mol->getProp("Column_2", mname);
  REQUIRE(mname == "4.0");
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
  REQUIRE(nSup2.length() == 4);
  nSup2.reset();
  nDone = 0;
  while (!nSup2.atEnd()) {
    mol = nSup2.next();
    nDone++;
    delete mol;
  }
  REQUIRE(nDone == nSup2.length());
}

TEST_CASE("testSmilesWriter") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";

  SmilesMolSupplier *nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  std::string oname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/outSmiles_molsupplier.csv";

  STR_VECT propNames;
  propNames.push_back(std::string("Column_2"));
  SmilesWriter *writer = new SmilesWriter(oname, " ");
  writer->setProps(propNames);

  STR_VECT names;
  STR_VECT props;
  ROMol *mol = nSup->next();

  while (mol) {
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
    REQUIRE(mname == names[i]);
    REQUIRE(pval == props[i]);
    i++;
    delete mol;
    try {
      mol = nSup->next();
    } catch (FileParseException &) {
      break;
    }
  }
  REQUIRE(nSup->length() == writer->numMols());
  writer->close();
  delete writer;
  delete nSup;
}

TEST_CASE("testSDWriter") {
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
  REQUIRE(writer->numMols() == 16);
  writer->close();
  delete writer;

  // now read in the file we just finished writing

  SDMolSupplier reader(ofile);
  while (!reader.atEnd()) {
    ROMol *mol = reader.next();
    CHECK(mol->hasProp(common_properties::_Name));

    delete mol;
  }
}

TEST_CASE("testSDSupplierEnding") {
  // test the SD supplier to check if it properly handle the end of sd file
  // conditions should work fine if the sd file end with  a $$$$ followed by
  // blank line or no no blank lines
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
  REQUIRE(i == 6);
}

TEST_CASE("testSuppliersEmptyFile") {
  {  // contains no records
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty.sdf";
    SDMolSupplier reader(infile);
    REQUIRE(reader.atEnd());
  }
  {
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty.smi";
    SmilesMolSupplier smiSup(infile, ",", 0, -1);
    REQUIRE(smiSup.atEnd());
  }
  // tests for GitHub issue 19:
  {  // actually an empty file, throws an exception:
    std::string infile =
        rdbase + "/Code/GraphMol/FileParsers/test_data/empty2.sdf";
    CHECK_THROWS_AS(SDMolSupplier(infile), BadFileException);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    REQUIRE(reader.atEnd());
    CHECK_THROWS_AS(reader[0], FileParseException);
    CHECK(reader.length() == 0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    CHECK_THROWS_AS(reader[0], FileParseException);
    CHECK(reader.length() == 0);
  }
  {
    SDMolSupplier reader;
    reader.setData("");
    CHECK(reader.length() == 0);
  }
}

TEST_CASE("testCisTrans") {
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
    REQUIRE(mol);
    RDDepict::compute2DCoords(*mol);
    writer.write(*mol);
    delete mol;
  }
  writer.close();
  // do the round t;est
  // parse the sd file and write it out to smiles

  SDMolSupplier *reader;
  REQUIRE_NOTHROW(reader = new SDMolSupplier("cisTrans_molsupplier.sdf"));
  REQUIRE(reader);
  while (!reader->atEnd()) {
    ROMol *mol = reader->next();
    REQUIRE(mol->hasProp(common_properties::_Name));
    delete mol;
  }
  delete reader;
}

TEST_CASE("testStereoRound") {
  // - we will read ina bunch of cdk2 smiles with stereo on them
  // - generate the canonical smiles for each one
  // - generate 2D coordinates, write to an sdf file
  // - read the sdf file back in and compare the canonical smiles
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/cdk2_stereo.csv";
  SmilesMolSupplier *smiSup;
  REQUIRE_NOTHROW(smiSup =
                      new SmilesMolSupplier(infile, ",", 0, 1, false, true));
  REQUIRE(smiSup);
  std::map<std::string, std::string> nameSmi;
  std::string ofile =
      rdbase +
      "/Code/GraphMol/FileParsers/test_data/cdk2_stereo_molsupplier.sdf";
  auto *writer = new SDWriter(ofile);

  while (!smiSup->atEnd()) {
    ROMol *mol = smiSup->next();

    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    nameSmi[mname] = MolToSmiles(*mol, 1);

    ROMol *nmol = SmilesToMol(nameSmi[mname]);

    std::string nsmi = MolToSmiles(*nmol, 1);
    CHECK(nameSmi[mname] == nsmi);

    RDDepict::compute2DCoords(*mol);
    writer->write(*mol);
    delete mol;
    delete nmol;
  }
  writer->close();
  delete smiSup;
  delete writer;

  // now read the SD file back in check if the canonical smiles are the same
  SDMolSupplier *reader;
  REQUIRE_NOTHROW(reader = new SDMolSupplier(ofile));
  REQUIRE(reader);

  while (!reader->atEnd()) {
    ROMol *mol = reader->next();

    std::string smiles = MolToSmiles(*mol, 1);
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    CHECK(nameSmi[mname] == smiles);
    delete mol;
  }
  delete reader;
}

TEST_CASE("testIssue226") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue226.sdf";
  SDMolSupplier sdsup(fname);

  ROMol *mol;

  mol = sdsup.next();
  REQUIRE(mol);
  REQUIRE(mol->hasProp("E1"));
  REQUIRE(mol->hasProp("E2"));
  delete mol;

  mol = sdsup.next();
  REQUIRE(mol);
  REQUIRE(mol->hasProp("E1"));
  REQUIRE(mol->hasProp("E2"));
  delete mol;
}

TEST_CASE("testTDTSupplier1") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  {
    TDTMolSupplier suppl(fname, "PN");
    unsigned int i = 0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1, prop2;
        REQUIRE(nmol->getNumAtoms() > 0);
        REQUIRE(nmol->hasProp("PN"));
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("MFCD"));

        nmol->getProp("PN", prop1);
        nmol->getProp(common_properties::_Name, prop2);
        REQUIRE(prop1 == prop2);

        // we didn't ask for 2D conformers, so there should be a property 2D:
        REQUIRE(nmol->hasProp(common_properties::TWOD));
        // and no conformer:
        REQUIRE(!nmol->getNumConformers());

        delete nmol;
        i++;
      }
    }
    REQUIRE(i == 10);
  }
  {
    std::ifstream strm(fname.c_str(), std::ios_base::binary);
    TDTMolSupplier suppl(&strm, false, "PN");
    unsigned int i = 0;
    while (!suppl.atEnd()) {
      ROMol *nmol = suppl.next();
      if (nmol) {
        std::string prop1, prop2;
        REQUIRE(nmol->getNumAtoms() > 0);
        REQUIRE(nmol->hasProp("PN"));
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("MFCD"));

        nmol->getProp("PN", prop1);
        nmol->getProp(common_properties::_Name, prop2);
        REQUIRE(prop1 == prop2);

        // we didn't ask for 2D conformers, so there should be a property 2D:
        REQUIRE(nmol->hasProp(common_properties::TWOD));
        // and no conformer:
        REQUIRE(!nmol->getNumConformers());

        delete nmol;
        i++;
      }
    }
    REQUIRE(i == 10);
  }
}

TEST_CASE("testTDTSupplier2") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  int i;
  std::string prop1, prop2;

  TDTMolSupplier suppl(fname, "PN", 2);
  i = 0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      REQUIRE(nmol->getNumAtoms() > 0);
      REQUIRE(nmol->hasProp("PN"));
      REQUIRE(nmol->hasProp(common_properties::_Name));
      REQUIRE(nmol->hasProp("MFCD"));

      nmol->getProp("PN", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      REQUIRE(prop1 == prop2);

      // we asked for 2D conformers, so there should be no property 2D:
      REQUIRE(!nmol->hasProp(common_properties::TWOD));
      // and a conformer:
      REQUIRE(nmol->getNumConformers() == 1);
      // with id "2":
      REQUIRE(nmol->beginConformers()->get()->getId() == 2);

      delete nmol;
      i++;
    }
  }
  REQUIRE(i == 10);
}

TEST_CASE("testTDTSupplier3") {
  int i;
  std::string prop1, prop2;

  TDTMolSupplier suppl;

  constexpr const char *data =
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
      REQUIRE(nmol->getNumAtoms() > 0);
      REQUIRE(nmol->hasProp("CAS"));
      REQUIRE(nmol->hasProp(common_properties::_Name));

      nmol->getProp("CAS", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      REQUIRE(prop1 == prop2);

      // no conformers should have been read:
      REQUIRE(nmol->getNumConformers() == 0);

      delete nmol;
      i++;
    }
  }
  REQUIRE(i == 4);
  REQUIRE(suppl.length() == 4);

  // now make sure we can grab earlier mols (was sf.net issue 1904170):
  ROMol *mol = suppl[0];
  REQUIRE(mol);
  delete mol;

  // make sure we can reset the supplier and still process it properly;
  suppl.setData(data, "CAS");

  i = 0;
  while (!suppl.atEnd()) {
    ROMol *nmol = suppl.next();
    if (nmol) {
      REQUIRE(nmol->getNumAtoms() > 0);
      REQUIRE(nmol->hasProp("CAS"));
      REQUIRE(nmol->hasProp(common_properties::_Name));

      nmol->getProp("CAS", prop1);
      nmol->getProp(common_properties::_Name, prop2);
      REQUIRE(prop1 == prop2);

      // no conformers should have been read:
      REQUIRE(nmol->getNumConformers() == 0);

      delete nmol;
      i++;
    }
  }
  REQUIRE(i == 4);
}

TEST_CASE("testSDSupplierFromText") {
  int i = 0;
  SDMolSupplier reader;

  constexpr const char *text =
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
    REQUIRE(mol->hasProp(common_properties::_Name));
    REQUIRE(mol->hasProp("ID"));
    i++;
    delete mol;
  }
  REQUIRE(i == 2);
}

TEST_CASE("testSDSupplierFromTextStrLax1") {
  constexpr const char *text =
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
      REQUIRE(mol->hasProp(common_properties::_Name));
      if (i == 0) {
        REQUIRE(!mol->hasProp("ID"));
      }
      REQUIRE(!mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    REQUIRE(i == 2);
  }
  // lax
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    REQUIRE(i == 2);
  }
}

TEST_CASE("testSDSupplierFromTextStrLax2") {
  constexpr const char *text =
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
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      REQUIRE(s == "Lig1");
      mol->getProp("ANOTHER_PROPERTY", s);
      REQUIRE(s ==
              "No blank line before dollars\n"
              "$$$$\n"
              "Structure1\n"
              "csChFnd70/05230312262D");
      i++;
      delete mol;
    }
    REQUIRE(i == 1);
  }
  // lax
  {
    SDMolSupplier reader;

    reader.setData(text, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      REQUIRE(s == "Lig2");
      mol->getProp("ANOTHER_PROPERTY", s);
      REQUIRE(s == "Value2");
      i++;
      delete mol;
    }
    REQUIRE(i == 1);
  }
}

TEST_CASE("testSDSupplierStrLax1") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/strictLax1.sdf";
  // strict
  {
    SDMolSupplier reader(fname, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      if (i == 0) {
        REQUIRE(!mol->hasProp("ID"));
      }
      REQUIRE(!mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    REQUIRE(i == 2);
  }
  // lax
  {
    SDMolSupplier reader(fname, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      i++;
      delete mol;
    }
    REQUIRE(i == 2);
  }
}

TEST_CASE("testSDSupplierStrLax2") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/strictLax2.sdf";
  // strict
  {
    SDMolSupplier reader(fname, true, true, true);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      REQUIRE(s == "Lig1");
      mol->getProp("ANOTHER_PROPERTY", s);
      REQUIRE(s ==
              "No blank line before dollars\n"
              "$$$$\n"
              "Structure1\n"
              "csChFnd70/05230312262D");
      i++;
      delete mol;
    }
    REQUIRE(i == 1);
  }
  // lax
  {
    SDMolSupplier reader(fname, true, true, false);

    int i = 0;
    while (!reader.atEnd()) {
      ROMol *mol = reader.next();
      REQUIRE(mol->hasProp(common_properties::_Name));
      REQUIRE(mol->hasProp("ID"));
      REQUIRE(mol->hasProp("ANOTHER_PROPERTY"));
      std::string s;
      mol->getProp("ID", s);
      REQUIRE(s == "Lig2");
      mol->getProp("ANOTHER_PROPERTY", s);
      REQUIRE(s == "Value2");
      i++;
      delete mol;
    }
    REQUIRE(i == 1);
  }
}

TEST_CASE("testIssue265") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NotThere.sdf";
  CHECK_THROWS_AS(SDMolSupplier(fname), BadFileException);

  CHECK_THROWS_AS(SmilesMolSupplier(fname), BadFileException);

  CHECK_THROWS_AS(TDTMolSupplier(fname), BadFileException);
}

TEST_CASE("testSDErrorHandling") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors1.sdf";
  SDMolSupplier *sdsup;
  ROMol *nmol = nullptr;

  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  REQUIRE(!sdsup->atEnd());
  nmol = sdsup->next();
  REQUIRE(nmol);
  REQUIRE(!nmol->hasProp("ID"));
  delete sdsup;
  delete nmol;

  // case 2: can't be sanitized
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors2.sdf";
  sdsup = new SDMolSupplier(fname);
  REQUIRE(!sdsup->atEnd());
  nmol = sdsup->next();
  REQUIRE(!nmol);
  REQUIRE(sdsup->atEnd());
  delete sdsup;
  delete nmol;

  // entry 3: bad number of atoms
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors3.sdf";
  sdsup = new SDMolSupplier(fname);
  REQUIRE(!sdsup->atEnd());
  nmol = sdsup->next();
  REQUIRE(!nmol);
  REQUIRE(sdsup->atEnd());
  delete sdsup;
  delete nmol;

  // entry 4: bad number of bonds
  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/sdErrors4.sdf";
  sdsup = new SDMolSupplier(fname);
  REQUIRE(!sdsup->atEnd());
  nmol = sdsup->next();
  REQUIRE(!nmol);
  REQUIRE(sdsup->atEnd());
  delete sdsup;
  delete nmol;
}

TEST_CASE("testIssue381") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue381.sdf";
  SDMolSupplier *sdsup;

  ROMol *nmol = nullptr;
  int count;

  // entry 1: bad properties
  sdsup = new SDMolSupplier(fname);
  REQUIRE(!sdsup->atEnd());
  count = 0;
  while (!sdsup->atEnd()) {
    nmol = sdsup->next();
    if (nmol) {
      delete nmol;
    }
    count++;
  }
  REQUIRE(sdsup->atEnd());
  REQUIRE(count == 9);

  REQUIRE(sdsup->length() == 9);

  delete sdsup;
}

TEST_CASE("testSetStreamIndices") {
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
  REQUIRE(!sdsup->atEnd());
  REQUIRE(sdsup->length() == 16);

  count = 0;
  while (!sdsup->atEnd()) {
    nmol = sdsup->next();
    if (nmol) {
      delete nmol;
    }
    count++;
  }
  REQUIRE(sdsup->atEnd());
  REQUIRE(count == 16);

  REQUIRE(sdsup->length() == 16);

  delete sdsup;
}

TEST_CASE("testMixIterAndRandom") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/esters.sdf";

  SDMolSupplier *sdsup;
  ROMol *mol;
  std::string name;

  sdsup = new SDMolSupplier(fname);
  REQUIRE(sdsup);
  unsigned int i = 0;
  while (!sdsup->atEnd()) {
    mol = sdsup->next();
    if (mol) {
      REQUIRE(mol->hasProp("ID"));
      delete mol;
    }
    i++;
  }
  REQUIRE(i == 6);
  REQUIRE(sdsup->length() == 6);

  delete sdsup;
  sdsup = new SDMolSupplier(fname);
  REQUIRE(sdsup);
  REQUIRE(sdsup->length() == 6);

  mol = sdsup->next();
  REQUIRE(mol);
  REQUIRE(mol->hasProp("ID"));
  mol->getProp("ID", name);
  REQUIRE(name == "Lig1");
  delete mol;

  mol = (*sdsup)[0];
  REQUIRE(mol);
  REQUIRE(mol->hasProp("ID"));
  mol->getProp("ID", name);
  REQUIRE(name == "Lig1");
  delete mol;

  sdsup->reset();
  mol = sdsup->next();
  REQUIRE(mol);
  REQUIRE(mol->hasProp("ID"));
  mol->getProp("ID", name);
  REQUIRE(name == "Lig1");
  delete mol;
  mol = sdsup->next();
  REQUIRE(mol);
  REQUIRE(mol->hasProp("ID"));
  mol->getProp("ID", name);
  REQUIRE(name == "Lig2");
  delete mol;
  delete sdsup;

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  SmilesMolSupplier *nSup;
  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  REQUIRE(nSup);
  REQUIRE(nSup->length() == 10);
  mol = (*nSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  REQUIRE(nSup);
  mol = (*nSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(nSup->length() == 10);
  delete mol;
  delete nSup;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  REQUIRE(nSup);
  mol = nSup->next();
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(nSup->length() == 10);
  delete mol;
  mol = (*nSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(nSup->length() == 10);
  delete mol;
  mol = nSup->next();
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 20);
  delete nSup;
  delete mol;

  nSup = new SmilesMolSupplier(fname, ",", 1, 0, false);
  REQUIRE(nSup);
  CHECK_THROWS_AS((*nSup)[20], FileParseException);
  delete nSup;

  fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
  TDTMolSupplier *tSup;
  tSup = new TDTMolSupplier(fname);
  REQUIRE(tSup);
  REQUIRE(tSup->length() == 10);
  mol = (*tSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  REQUIRE(tSup);
  mol = (*tSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(tSup->length() == 10);
  delete mol;
  delete tSup;

  tSup = new TDTMolSupplier(fname);
  REQUIRE(tSup);
  mol = tSup->next();
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(tSup->length() == 10);
  delete mol;

  mol = (*tSup)[0];
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 9);
  REQUIRE(tSup->length() == 10);
  delete mol;

  mol = tSup->next();
  REQUIRE(mol);
  delete mol;

  mol = tSup->next();
  REQUIRE(mol);
  delete mol;

  mol = tSup->next();
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 10);
  delete tSup;
  delete mol;

  tSup = new TDTMolSupplier(fname);
  REQUIRE(tSup);
  CHECK_THROWS_AS((*tSup)[20], FileParseException);
  delete tSup;
}

TEST_CASE("testRemoveHs") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/withHs.sdf";

  SDMolSupplier sdsup(fname);
  ROMol *nmol;

  nmol = sdsup.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 23);
  delete nmol;
  nmol = sdsup.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 28);
  delete nmol;

  std::cerr << "build:" << std::endl;
  SDMolSupplier sdsup2(fname, true, false);
  nmol = sdsup2.next();
  REQUIRE(nmol);

  REQUIRE(nmol->getNumAtoms() == 39);
  delete nmol;
  nmol = sdsup2.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 30);
  delete nmol;
}

TEST_CASE("testGetItemText") {
  std::string fname;

  ROMol *mol1, *mol2;
  std::string molB, smiles;

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);
    REQUIRE(sdsup.length() == 16);

    molB = sdsup.getItemText(0);
    mol1 = sdsup[0];
    REQUIRE(mol1);
    mol2 = MolBlockToMol(molB);
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    // make sure getItemText() doesn't screw up the current position:
    molB = sdsup.getItemText(10);
    mol1 = sdsup.next();
    molB = sdsup.getItemText(1);
    REQUIRE(mol1);
    mol2 = MolBlockToMol(molB);
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = sdsup.getItemText(15);
    mol1 = sdsup[15];
    mol2 = MolBlockToMol(molB);
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    CHECK_THROWS_AS(sdsup.getItemText(16), FileParseException);

    CHECK_THROWS_AS(sdsup.getItemText(20), FileParseException);
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
    SDMolSupplier sdsup(fname);

    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = sdsup.getItemText(0);
    mol2 = MolBlockToMol(molB);
    REQUIRE(mol2);
    mol1 = sdsup[0];
    REQUIRE(mol1);
    REQUIRE(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;

    molB = sdsup.getItemText(5);
    mol2 = MolBlockToMol(molB);
    REQUIRE(mol2);
    REQUIRE(mol2->getNumAtoms() == 16);
    mol1 = sdsup[5];
    REQUIRE(mol1);
    REQUIRE(mol2->getNumAtoms() == mol1->getNumAtoms());
    delete mol1;
    delete mol2;
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);
    REQUIRE(smisup.length() == 10);

    molB = smisup.getItemText(0);
    REQUIRE(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");
    mol1 = smisup[0];
    REQUIRE(mol1);
    delete mol1;

    molB = smisup.getItemText(5);
    REQUIRE(
        molB ==
        "6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");
    mol1 = smisup.next();
    REQUIRE(mol1);
    REQUIRE(mol1->getNumAtoms() == 20);
    delete mol1;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = smisup.getItemText(8);
    REQUIRE(molB == "9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    REQUIRE(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");

    mol1 = smisup[0];
    REQUIRE(mol1);
    smiles = MolToSmiles(*mol1, 1);
    REQUIRE(smiles == "CC1=CC(=O)C=CC1=O");
    REQUIRE(mol1->getNumAtoms() == 9);
    delete mol1;
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);

    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = smisup.getItemText(0);
    REQUIRE(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");

    molB = smisup.getItemText(5);
    REQUIRE(
        molB ==
        "6, OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br, 87.74");

    molB = smisup.getItemText(8);
    REQUIRE(molB == "9, CC(=NO)C(C)=NO, 65.18");
    molB = smisup.getItemText(9);
    REQUIRE(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
    SmilesMolSupplier smisup(fname, ",", 1, 0, false);

    // make sure getItemText() flags EOF
    // (this was sf.net issue 3299878)
    molB = smisup.getItemText(0);
    REQUIRE(molB == "1, CC1=CC(=O)C=CC1=O, 34.14");

    ROMol *m = smisup[9];
    REQUIRE(m);
    delete m;
    REQUIRE(smisup.atEnd());
    molB = smisup.getItemText(9);
    REQUIRE(molB == "10, C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3, 0.00");
    REQUIRE(smisup.atEnd());
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    // make sure getItemText() works if we haven't read at all from the
    // supplier:
    // (this was sf.net issue 2632960)
    molB = tdtsup.getItemText(0);
    REQUIRE(molB != "");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    REQUIRE(tdtsup.length() == 10);

    molB = tdtsup.getItemText(0);
    REQUIRE(molB != "");

    mol1 = tdtsup[0];
    REQUIRE(mol1);
    smiles = MolToSmiles(*mol1, 1);
    REQUIRE(smiles == "Cc1nnc(N)nc1C");
    REQUIRE(mol1->getNumAtoms() == 9);
    delete mol1;

    // make sure getItemText doesn't screw up next()
    molB = tdtsup.getItemText(5);
    mol1 = tdtsup.next();
    REQUIRE(mol1);
    REQUIRE(mol1->getNumAtoms() == 9);
    smiles = MolToSmiles(*mol1, 1);
    REQUIRE(smiles == "Cc1n[nH]c(=O)nc1N");
    delete mol1;

    // make sure getItemText() works on the last molecule
    // (this was sf.net issue 1874882
    molB = tdtsup.getItemText(9);
    REQUIRE(molB != "");
    REQUIRE(molB.substr(0, 12) == "$SMI<Cc1n[nH");
  }

  {
    fname = rdbase + "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    TDTMolSupplier tdtsup(fname);
    REQUIRE(tdtsup.length() == 10);

    ROMol *mol = tdtsup[9];
    REQUIRE(mol);
    delete mol;
    REQUIRE(tdtsup.atEnd());

    // (this was sf.net issue 3299878
    molB = tdtsup.getItemText(9);
    REQUIRE(molB != "");
    REQUIRE(molB.substr(0, 12) == "$SMI<Cc1n[nH");
    REQUIRE(tdtsup.atEnd());
  }
}

TEST_CASE("testForwardSDSupplier") {
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
      REQUIRE((nmol || sdsup.atEnd()));
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    REQUIRE(i == 16);
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
    REQUIRE(i == 998);
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
    REQUIRE(i == 997);
  }
  // looks good, now do a supplier:
  {
    gzstream strm(fname2);

    ForwardSDMolSupplier sdsup(&strm, false);
    unsigned int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      if (nmol) {
        REQUIRE(nmol->hasProp(common_properties::_Name));
        REQUIRE(nmol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
        delete nmol;
        i++;
      }
    }
    REQUIRE(i == 16);
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
    REQUIRE(i == 1663);
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
    REQUIRE(i == 1663);
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
    REQUIRE(i == 16);
  }
#endif
#endif  // RDK_BUILD_MAEPARSER_SUPPORT
}

TEST_CASE("testMissingCRSDSupplier") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/missingCR.sdf";
  SDMolSupplier reader(infile);
  auto *mol = reader.next();
  delete mol;
  REQUIRE(reader.atEnd());
}

TEST_CASE("testIssue3482695") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3482695.sdf";
  SDMolSupplier reader(infile);
  ROMol *nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 0);
  REQUIRE(nmol->hasProp("test"));
  delete nmol;
}

TEST_CASE("testIssue3525673") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525673.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 37);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 58);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  delete nmol;

  nmol = reader.next();
  REQUIRE(!nmol);  // broken due to 'foo' in counts line!

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 58);
  delete nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  delete nmol;
}

TEST_CASE("testBlankLinesInProps") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/BlankPropLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 19);
  REQUIRE(nmol->hasProp("MultiLineProperty1"));
  nmol->getProp("MultiLineProperty1", pval);
  REQUIRE(pval == "foo\nbar\n \nbaz");
  REQUIRE(nmol->hasProp("MultiLineProperty2"));
  REQUIRE(!(nmol->hasProp("fooprop")));
  nmol->getProp("MultiLineProperty2", pval);
  REQUIRE(pval == "foo\n>  <fooprop>\nbaz\n ");
  delete nmol;
}

TEST_CASE("testSkipLines") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/SkipLines.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;
  std::string pval;

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 1);
  REQUIRE(nmol->hasProp("prop1"));
  delete nmol;
}

TEST_CASE("testGitHub23") {
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

TEST_CASE("testGitHub88") {
  std::string infile =
      rdbase + "/Code/GraphMol/FileParsers/test_data/github88.v3k.sdf";
  std::ifstream ins(infile.c_str());
  ForwardSDMolSupplier reader(&ins, false);
  ROMol *nmol;

  nmol = reader.next();
  REQUIRE(nmol);
  REQUIRE(nmol->getNumAtoms() == 8);
  REQUIRE(nmol->hasProp("prop1"));
  std::string pval;
  nmol->getProp("prop1", pval);
  REQUIRE(pval == "4");
  delete nmol;
}

TEST_CASE("testGitHub2285") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/github2285.sdf";

  std::vector<std::string> smiles;
  {
    SDMolSupplier sdsup(fname);
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      REQUIRE(nmol);
      smiles.push_back(MolToSmiles(*nmol));
      delete nmol;
    }
  }
  {
    SDMolSupplier sdsup(fname, true, false);
    int i = 0;
    while (!sdsup.atEnd()) {
      ROMol *nmol = sdsup.next();
      REQUIRE(nmol);
      ROMol *m = MolOps::removeHs(*nmol);
      REQUIRE(MolToSmiles(*m) == smiles[i++]);
      delete nmol;
      delete m;
    }
    REQUIRE(i > 0);
  }
}

TEST_CASE("testGitHub2479") {
  constexpr const char *smiles1 = R"DATA(smiles id
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
        REQUIRE(mol);
      }
      ++cnt;
    }
    REQUIRE(cnt == 5);
  }

  constexpr const char *sdf1 = R"SDF(
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
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(!mol2);
    REQUIRE(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf1);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(!mol2);
    REQUIRE(!suppl.atEnd());
    REQUIRE(!suppl.getEOFHitOnRead());
    std::unique_ptr<ROMol> mol3(suppl.next());
    REQUIRE(!mol3);
    REQUIRE(suppl.atEnd());
    REQUIRE(suppl.getEOFHitOnRead());
  }

  // truncated file1
  constexpr const char *sdf2 = R"SDF(
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
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(!mol2);
    std::unique_ptr<ROMol> mol3(suppl.next());
    REQUIRE(!mol3);
    REQUIRE(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf2);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(!mol2);
    REQUIRE(!suppl.atEnd());
    REQUIRE(!suppl.getEOFHitOnRead());
    std::unique_ptr<ROMol> mol3(suppl.next());
    REQUIRE(!mol3);
    REQUIRE(suppl.atEnd());
    REQUIRE(!suppl.getEOFHitOnRead());
  }
  // truncated file2
  constexpr const char *sdf3 = R"SDF(
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
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(mol2);
    REQUIRE(suppl.atEnd());
  }
  {
    std::stringstream iss(sdf3);
    ForwardSDMolSupplier suppl(&iss, false);
    std::unique_ptr<ROMol> mol1(suppl.next());
    REQUIRE(mol1);
    std::unique_ptr<ROMol> mol2(suppl.next());
    REQUIRE(mol2);
    REQUIRE(suppl.atEnd());
  }
}

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
TEST_CASE("testGitHub2881") {
  constexpr const char *data = R"DATA(f_m_ct {
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
    std::unique_ptr<ROMol> mol(suppl.next());
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 15);
    REQUIRE(mol->getNumBonds() == 17);
  }
}
#endif

TEST_CASE("testGitHub3517") {
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";

  SDMolSupplier sdsup(fname);
  REQUIRE(!sdsup.atEnd());
  size_t l = sdsup.length();
  REQUIRE(l > 0);
  REQUIRE(!sdsup.atEnd());
}

TEST_CASE(
    "GitHub Issue #9014: SDMolSupplier enters an infinite loop if number of SGroups is incorrect") {
  constexpr const char *molblock = R"CTAB(
                    2D

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.299038 0.750000 0.000000 0
M  V30 2 C 1.299038 2.250000 0.000000 0
M  V30 3 O 2.598076 3.000000 0.000000 0
M  V30 4 C 2.598076 -0.000000 0.000000 0
M  V30 5 C 0.000000 0.000000 0.000000 0
M  V30 6 H 3.897114 0.750000 0.000000 0
M  V30 7 H 2.598076 4.500000 0.000000 0
M  V30 8 H 2.598076 -1.500000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 5 1 4 6
M  V30 6 1 3 7
M  V30 7 1 4 8
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(8 1 2 3 4 5 6 7 8)
M  V30 END SGROUP
M  V30 END CTAB
M  END
$$$$)CTAB";

  v2::FileParsers::SDMolSupplier sdsup;
  v2::FileParsers::MolFileParserParams p{.strictParsing = false};

  sdsup.setData(molblock, p);

  std::unique_ptr<ROMol> mol;
  auto parser_thread = std::async(std::launch::async, [&sdsup, &mol]() {
    // This might run forever if the code is bugged...
    mol = sdsup.next();
  });

  // wait a few seconds (allow extra time for debugging),
  // and send sigint if the parser thread hasn't finished
#ifndef DNDEBUG
  auto timeout = std::chrono::seconds(900);
#else
  auto timeout = std::chrono::seconds(10);
#endif
  auto status = parser_thread.wait_for(timeout);
  if (status != std::future_status::ready) {
    // There's no safe way to terminate a rogue thread in C++,
    //  so just send Ctrl+C and abort the whole test run
    raise(SIGINT);
  } else {
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 5);  // only heavy atoms are kept
    CHECK(getSubstanceGroups(*mol).size() == 1);
  }
}

TEST_CASE("Read SD properties till last '>'") {
  auto m = v2::SmilesParse::MolFromSmiles("C");
  REQUIRE(m);

  constexpr const char *prop_name = "654 > 321";
  m->setProp(prop_name, "this is not important");

  auto molblock = SDWriter::getText(*m);
  REQUIRE_THAT(molblock, Catch::Matchers::ContainsSubstring(prop_name));

  SDMolSupplier supplier;
  supplier.setData(molblock);

  auto m2 = supplier.next();
  REQUIRE(m2);
  CHECK(m2->hasProp(prop_name));
}