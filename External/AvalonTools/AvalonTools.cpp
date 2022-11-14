//
//  Copyright (C) 2008-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Created by Greg Landrum, July 2008
//

#include <DataStructs/ExplicitBitVect.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Geometry/point.h>
#include "AvalonTools.h"

extern "C" {
#include "local.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "ssmatch.h"
#include "smi2mol.h"
#include "canonizer.h"
#include "layout.h"
#include "struchk.h"

extern int RunStruchk(struct reaccs_molecule_t **mpp,
                      struct data_line_t *data_list);
extern void ClearParameters();
}

// already defined in struchk.c
// FILE *log_file=NULL;

namespace AvalonTools {
using namespace RDKit;
namespace {
int *getCountFp(struct reaccs_molecule_t *molPtr, unsigned int bitFlags,
                bool isQuery, unsigned int nBytes) {
  PRECONDITION(molPtr, "bad molecule");
  int *res = TypeAlloc(nBytes * sizeof(int), int);
  memset(res, 0, nBytes * sizeof(int));
  SetFingerprintCountsWithFocus(molPtr, res, static_cast<int>(nBytes),
                                static_cast<int>(bitFlags),
                                static_cast<int>(isQuery), 0, 0);
  if (!isQuery) {
    SetFingerprintCountsWithFocus(
        molPtr, res, static_cast<int>(nBytes), static_cast<int>(bitFlags),
        static_cast<int>(1), ACCUMULATE_BITS | USE_DY_AROMATICITY, 0);
  }
  return res;
}
char *getFp(struct reaccs_molecule_t *molPtr, unsigned int bitFlags,
            bool isQuery, unsigned int nBytes) {
  PRECONDITION(molPtr, "bad molecule");
  while (nBytes % 4) {
    ++nBytes;
  }
  char *fingerprint = TypeAlloc(nBytes, char);
  SetFingerprintBits(molPtr, fingerprint, static_cast<int>(nBytes),
                     static_cast<int>(bitFlags), static_cast<int>(isQuery), 0);
  if (!isQuery) {
    SetFingerprintBits(molPtr, fingerprint, static_cast<int>(nBytes),
                       static_cast<int>(bitFlags), static_cast<int>(0),
                       ACCUMULATE_BITS | USE_DY_AROMATICITY);
  }
  return fingerprint;
}
void reaccsToFingerprint(struct reaccs_molecule_t *molPtr,
                         std::vector<boost::uint32_t> &res,
                         unsigned int bitFlags = 32767U, bool isQuery = false,
                         bool resetVect = true, unsigned int nBytes = 64) {
  if (resetVect) {
    res.clear();
  }
  char *fingerprint = getFp(molPtr, bitFlags, isQuery, nBytes);
  for (unsigned int i = 0; i < nBytes; i += 4) {
    boost::uint32_t word;
    word = fingerprint[i] | (fingerprint[i + 1] << 8) |
           (fingerprint[i + 2] << 16) | (fingerprint[i + 3] << 24);
    res.push_back(word);
  }

  MyFree(fingerprint);
};

void reaccsToCounts(struct reaccs_molecule_t *molPtr,
                    SparseIntVect<boost::uint32_t> &res,
                    unsigned int bitFlags = 32767U, bool isQuery = false,
                    unsigned int nBytes = 64) {
  PRECONDITION(molPtr, "bad molecule");
  PRECONDITION(res.getLength() >= nBytes, "res too small");

  int *fingerprint = getCountFp(molPtr, bitFlags, isQuery, nBytes);

  for (unsigned int i = 0; i < nBytes; ++i) {
    res.setVal(i, fingerprint[i]);
  }
  MyFree((char *)fingerprint);
};

void reaccsToFingerprint(struct reaccs_molecule_t *molPtr, ExplicitBitVect &res,
                         unsigned int bitFlags = 32767U, bool isQuery = false,
                         bool resetVect = true, unsigned int nBytes = 64) {
  PRECONDITION(molPtr, "bad molecule");
  PRECONDITION(res.getNumBits() >= nBytes * 8U, "res too small");
  if (resetVect) {
    res.clearBits();
  }

  char *fingerprint = getFp(molPtr, bitFlags, isQuery, nBytes);

  for (unsigned int i = 0; i < nBytes; ++i) {
    char byte = fingerprint[i];
    if (byte) {
      char mask = 1;
      for (int j = 0; j < 8; ++j) {
        if (byte & mask) {
          res.setBit(i * 8 + j);
        }
        mask = mask << 1;
      }
    }
  }
  MyFree(fingerprint);
};

struct reaccs_molecule_t *reaccsGetCoords(struct reaccs_molecule_t *molPtr) {
  PRECONDITION(molPtr, "bad molecule");

  RecolorMolecule(molPtr);
  struct reaccs_molecule_t *res = LayoutMolecule(molPtr);
  POSTCONDITION(res, "could not layout molecule");
  return res;
};

struct reaccs_molecule_t *molToReaccs(const ROMol &mol) {
  std::string molB = MolToMolBlock(mol, true);
  Utils::LocaleSwitcher ls;
  struct reaccs_molecule_t *res = MolStr2Mol((char *)molB.c_str());
  POSTCONDITION(res, "could not build a molecule");
  return res;
}

struct reaccs_molecule_t *stringToReaccs(const std::string &data,
                                         bool isSmiles) {
  struct reaccs_molecule_t *res;
  if (isSmiles) {
    res = SMIToMOL(data.c_str(), DY_AROMATICITY);
  } else {
    Utils::LocaleSwitcher ls;
    res = MolStr2Mol((char *)data.c_str());
  }
  if (!res) {
    if (isSmiles) {
      BOOST_LOG(rdErrorLog)
          << "ERROR could not build molecule from smiles: " << data
          << std::endl;
    } else {
      BOOST_LOG(rdErrorLog)
          << "ERROR could not build molecule from molblock: \n"
          << data << std::endl;
    }
  }
  return res;
}

}  // end of anonymous namespace

std::string getCanonSmiles(ROMol &mol, int flags) {
  if (flags == -1) {
    flags = DB_STEREO | CENTER_STEREO;
  }
  std::string res;
  if (!mol.getNumConformers()) {
    std::string rdSmi = MolToSmiles(mol, true);
    res = getCanonSmiles(rdSmi, true, flags);
  } else {
    std::string rdMB = MolToMolBlock(mol);
    res = getCanonSmiles(rdMB, false, flags);
  }
  return res;
}

void getAvalonCountFP(const ROMol &mol, SparseIntVect<boost::uint32_t> &res,
                      unsigned int nBits, bool isQuery, bool resetVect,
                      unsigned int bitFlags) {
  (void)resetVect;
  struct reaccs_molecule_t *mp = molToReaccs(mol);
  reaccsToCounts(mp, res, bitFlags, isQuery, nBits);
  FreeMolecule(mp);
}

void getAvalonFP(const ROMol &mol, ExplicitBitVect &res, unsigned int nBits,
                 bool isQuery, bool resetVect, unsigned int bitFlags) {
  if (nBits % 8) {
    BOOST_LOG(rdWarningLog)
        << "Warning: number of bits (" << nBits
        << ") is not evenly divisible by 8. Rounding to the nearest byte."
        << std::endl;
  }
  unsigned int nBytes = nBits / 8;
  struct reaccs_molecule_t *mp = molToReaccs(mol);
  reaccsToFingerprint(mp, res, bitFlags, isQuery, resetVect, nBytes);
  FreeMolecule(mp);
}
void getAvalonFP(const ROMol &mol, std::vector<boost::uint32_t> &res,
                 unsigned int nBits, bool isQuery, bool resetVect,
                 unsigned int bitFlags) {
  if (nBits % 8) {
    BOOST_LOG(rdWarningLog)
        << "Warning: number of bits (" << nBits
        << ") is not evenly divisible by 8. Rounding to the nearest byte."
        << std::endl;
  }
  unsigned int nBytes = nBits / 8;
  struct reaccs_molecule_t *mp = molToReaccs(mol);
  reaccsToFingerprint(mp, res, bitFlags, isQuery, resetVect, nBytes);
  FreeMolecule(mp);
}

unsigned int set2DCoords(ROMol &mol, bool clearConfs) {
  bool origMolHasHs = false;
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() == 1 && !atom->getIsotope()) {
      origMolHasHs = true;
      break;
    }
  }
  auto smiles = MolToSmiles(mol);
  struct reaccs_molecule_t *mp = stringToReaccs(smiles, true);
  struct reaccs_molecule_t *mp2 = reaccsGetCoords(mp);
  TEST_ASSERT(mp2->n_atoms >= mol.getNumAtoms());

  auto *conf = new RDKit::Conformer(mol.getNumAtoms());
  conf->set3D(false);

  std::vector<int> matchedIndices(mol.getNumAtoms());
  // the toolkit may add chiral Hs without putting them at the end of the list
  //    we need to be able to ignore them:
  if (origMolHasHs || mp2->n_atoms > mol.getNumAtoms()) {
    // the toolkit may have rearranged the atoms, we need to do a substructure
    // match to figure out what's what
    char *mb = MolToMolStr(mp2);
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<RWMol> avmol(MolBlockToMol(mb, sanitize, removeHs));
    MyFree(mb);
    CHECK_INVARIANT(avmol, "could not parse mol block from avalon toolkit");
    TEST_ASSERT(avmol);
    auto match = SubstructMatch(*avmol, mol);
    CHECK_INVARIANT(
        !match.empty(),
        "no substructure match found between avalon mol and input mol");
    for (const auto &pr : match[0]) {
      matchedIndices[pr.first] = pr.second;
    }
  } else {
    // Atoms in the intermediate smiles representation may be ordered
    // differently compared to the original input molecule.
    // Make sure that output coordinates are assigned in the correct order.
    std::vector<unsigned int> atomOrdering;
    mol.getProp(common_properties::_smilesAtomOutputOrder, atomOrdering);
    for (size_t ai = 0; ai < mol.getNumAtoms(); ++ai) {
      matchedIndices[ai] = atomOrdering[ai];
    }
  }
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    auto x = mp2->atom_array[matchedIndices[i]].x;
    auto y = mp2->atom_array[matchedIndices[i]].y;
    RDGeom::Point3D loc(x, y, 0.);
    conf->setAtomPos(i, loc);
  }

  unsigned int res;
  if (clearConfs) {
    mol.clearConformers();
    conf->setId(0);
    mol.addConformer(conf);
    res = 0;
  } else {
    res = mol.addConformer(conf, true);
  }

  FreeMolecule(mp);
  FreeMolecule(mp2);

  return res;
}
std::string set2DCoords(const std::string &data, bool isSmiles) {
  struct reaccs_molecule_t *mp = stringToReaccs(data, isSmiles);
  std::string res = "";
  if (mp) {
    struct reaccs_molecule_t *mp2 = reaccsGetCoords(mp);
    Utils::LocaleSwitcher ls;
    char *molB = MolToMolStr(mp2);
    res = molB;
    FreeMolecule(mp);
    FreeMolecule(mp2);
    MyFree(molB);
  }
  return res;
}

std::string getCanonSmiles(const std::string &data, bool isSmiles, int flags) {
  if (flags == -1) {
    flags = DB_STEREO | CENTER_STEREO;
  }
  char *smiles = nullptr, *canSmiles = nullptr;
  if (!isSmiles) {
    struct reaccs_molecule_t *mp = stringToReaccs(data, isSmiles);
    if (mp) {
      smiles = MOLToSMI(mp, ISOMERIC_SMILES);
      FreeMolecule(mp);
      if (smiles) {
        canSmiles = CanSmiles(smiles, flags);
        MyFree(smiles);
      }
    }
  } else {
    canSmiles = CanSmiles((char *)data.c_str(), flags);
  }
  std::string res = "";
  if (canSmiles) {
    res = canSmiles;
    MyFree(canSmiles);
  } else {
    BOOST_LOG(rdErrorLog) << "ERROR: no smiles generated for molecule."
                          << std::endl;
  }
  return res;
}

void getAvalonCountFP(const std::string &data, bool isSmiles,
                      SparseIntVect<boost::uint32_t> &res, unsigned int nBits,
                      bool isQuery, unsigned int bitFlags) {
  struct reaccs_molecule_t *mp = stringToReaccs(data, isSmiles);
  if (mp) {
    reaccsToCounts(mp, res, bitFlags, isQuery, nBits);
    FreeMolecule(mp);
  } else {
    BOOST_LOG(rdErrorLog) << "ERROR: no fingeprint generated for molecule."
                          << std::endl;
  }
}
void getAvalonFP(const std::string &data, bool isSmiles, ExplicitBitVect &res,
                 unsigned int nBits, bool isQuery, bool resetVect,
                 unsigned int bitFlags) {
  if (nBits % 8) {
    BOOST_LOG(rdWarningLog)
        << "Warning: number of bits (" << nBits
        << ") is not evenly divisible by 8. Rounding to the nearest byte."
        << std::endl;
  }
  unsigned int nBytes = nBits / 8;
  struct reaccs_molecule_t *mp = stringToReaccs(data, isSmiles);
  if (mp) {
    reaccsToFingerprint(mp, res, bitFlags, isQuery, resetVect, nBytes);
    FreeMolecule(mp);
  } else {
    BOOST_LOG(rdErrorLog) << "ERROR: no fingeprint generated for molecule."
                          << std::endl;
  }
}
void getAvalonFP(const std::string &data, bool isSmiles,
                 std::vector<boost::uint32_t> &res, unsigned int nBits,
                 bool isQuery, bool resetVect, unsigned int bitFlags) {
  if (nBits % 8) {
    BOOST_LOG(rdWarningLog)
        << "Warning: number of bits (" << nBits
        << ") is not evenly divisible by 8. Rounding to the nearest byte."
        << std::endl;
  }
  unsigned int nBytes = nBits / 8;
  struct reaccs_molecule_t *mp = stringToReaccs(data, isSmiles);
  if (mp) {
    reaccsToFingerprint(mp, res, bitFlags, isQuery, resetVect, nBytes);
    FreeMolecule(mp);
  } else {
    BOOST_LOG(rdErrorLog) << "ERROR: no fingeprint generated for molecule."
                          << std::endl;
  }
}

int _checkMolWrapper(struct reaccs_molecule_t **mpp) {
  if (!*mpp) {
    return BAD_MOLECULE;
  }
  int res;
  struct reaccs_molecule_t *tmp = *mpp;
  res = RunStruchk(mpp, nullptr);
  if (*mpp != tmp) {
    FreeMolecule(tmp);
  }
  return res;
}

/**
 * Wrapper around struchk.CheckMol
 * The molecule to check is passed in as a string. isSmiles
 * should be set to TRUE if the molecule is encoded as SMILES,
 * to FALSE if the molecule is encoded as sn MDL CTAB.
 * mp is an output parameter - it will point to the checked
 * molecule upon successful checking. In case of errors, mp may be 0.
 **/
int checkMolString(const std::string &data, const bool isSmiles,
                   struct reaccs_molecule_t **mp) {
  // clean msg list from previous call (if no previous call, freemsglist does
  // nothing)
  FreeMsgList();

  int errs = 0;
  if (isSmiles) {
    *mp = SMIToMOL(data.c_str(), DY_AROMATICITY);
  } else {
    Utils::LocaleSwitcher ls;
    *mp = MolStr2Mol((char *)data.c_str());
  }
  if (*mp) {
    errs = _checkMolWrapper(mp);
  } else {
    errs = BAD_MOLECULE;
  }
  return errs;
}

int initCheckMol(const std::string &optString) {
  // n.b. always add a cr to the end for safety
  auto *optBuffer = new char[optString.size() + 2];
  optString.copy(optBuffer, optString.size());
  optBuffer[optString.size() - 1] = '\n';
  optBuffer[optString.size()] = '\0';
  int res = InitCheckMol(optBuffer);
  delete[] optBuffer;
  return res;
}

std::string getCheckMolLog() {
  char *buf = GetMsgList();
  std::string res = buf;
  MyFree(buf);

  return res;
}

RDKit::ROMOL_SPTR checkMol(int &errs, RDKit::ROMol &inMol) {
  // clean msg list from previous call (if no previous call, freemsglist does
  // nothing)
  FreeMsgList();

  struct reaccs_molecule_t *mp;
  RDKit::ROMol *rMol = nullptr;
  mp = molToReaccs(inMol);
  errs = _checkMolWrapper(&mp);
  if (mp) {
    Utils::LocaleSwitcher ls;
    char *molStr = MolToMolStr(mp);
    FreeMolecule(mp);
    if (molStr) {
      rMol = MolBlockToMol(molStr);
      MyFree(molStr);
    }
  }
  return RDKit::ROMOL_SPTR(rMol);
}

RDKit::ROMOL_SPTR checkMol(int &errs, const std::string &data,
                           const bool isSmiles) {
  struct reaccs_molecule_t *mp;
  errs = checkMolString(data, isSmiles, &mp);
  if (mp) {
    Utils::LocaleSwitcher ls;
    char *molStr = MolToMolStr(mp);
    RDKit::ROMol *rMol = MolBlockToMol(molStr);
    FreeMolecule(mp);
    MyFree(molStr);
    return RDKit::ROMOL_SPTR(rMol);
  } else {
    return RDKit::ROMOL_SPTR();
  }
}

std::pair<std::string, int> checkMolString(const std::string &data,
                                           bool isSmiles) {
  struct reaccs_molecule_t *mp;
  int errs = checkMolString(data, isSmiles, &mp);
  std::string molStr;
  if (mp) {
    Utils::LocaleSwitcher ls;
    char *tmp = MolToMolStr(mp);
    molStr = std::string(tmp);
    FreeMolecule(mp);
    MyFree(tmp);
  } else {
    molStr = "";
  }
  return std::make_pair(molStr, errs);
}

void closeCheckMolFiles() {
  ClearParameters();
  CloseOpenFiles();
}

}  // namespace AvalonTools
