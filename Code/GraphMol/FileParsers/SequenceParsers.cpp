//
//  Copyright (C) 2015,2016 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>

#include <string>
#include <vector>
#include <map>

#include "SequenceParsers.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>

namespace RDKit {

static Atom *CreateAAAtom(RWMol *mol, const char *name,
                          AtomPDBResidueInfo &info) {
  Atom *atom = (Atom *)nullptr;

  if (name[0] == ' ' && name[1] == 'C') {
    atom = new Atom(6);
  } else if (name[0] == ' ' && name[1] == 'N') {
    atom = new Atom(7);
  } else if (name[0] == ' ' && name[1] == 'O') {
    atom = new Atom(8);
  } else if (name[0] == ' ' && name[1] == 'P') {
    atom = new Atom(15);
  } else if (name[0] == ' ' && name[1] == 'S') {
    atom = new Atom(16);
  } else if (name[0] == 'S' && name[1] == 'E') {
    atom = new Atom(34);
  } else {
    atom = new Atom(0);
  }
  mol->addAtom(atom, true, true);
  auto *copy = (AtomPDBResidueInfo *)info.copy();
  copy->setName(name);
  atom->setMonomerInfo(copy);

  unsigned int serno = info.getSerialNumber();
  info.setSerialNumber(serno + 1);
  return atom;
}

static void CreateAABond(RWMol *mol, Atom *beg, Atom *end, unsigned int order) {
  Bond *bond;
  if (order == 2) {
    bond = new Bond(Bond::DOUBLE);
  } else {
    bond = new Bond(Bond::SINGLE);
  }
  bond->setOwningMol(mol);
  bond->setBeginAtom(beg);
  bond->setEndAtom(end);
  mol->addBond(bond, true);
}

static void CreateAABackbone(RWMol *mol, Atom *&r1, Atom *&r2, Atom *&cb,
                             AtomPDBResidueInfo &info, int ldstereo) {
  r1 = CreateAAAtom(mol, " N  ", info);
  Atom *ca = CreateAAAtom(mol, " CA ", info);
  r2 = CreateAAAtom(mol, " C  ", info);
  Atom *o = CreateAAAtom(mol, " O  ", info);
  cb = CreateAAAtom(mol, " CB ", info);
  CreateAABond(mol, r1, ca, 1);
  CreateAABond(mol, ca, r2, 1);
  CreateAABond(mol, r2, o, 2);
  CreateAABond(mol, ca, cb, 1);

  if (ldstereo > 0) {  // L-stereo
    ca->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  } else if (ldstereo < 0) {  // D-stereo
    ca->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  }
}

// aa is a three letter PDB residue code
static void CreateAminoAcid(RWMol *mol, const char *aa, Atom *&r1, Atom *&r2,
                            Atom *&r3, AtomPDBResidueInfo &info) {
  Atom *atom[10];

  r1 = (Atom *)nullptr;
  r2 = (Atom *)nullptr;
  r3 = (Atom *)nullptr;

  int resno = info.getResidueNumber();
  info.setResidueNumber(resno + 1);
  info.setIsHeteroAtom(false);
  info.setResidueName(aa);

  // Standard amino acids before non-standard, in PDB code alphabetical order
  switch (aa[0]) {
    case 'A':
      if (!strcmp(aa, "ALA")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
      } else if (!strcmp(aa, "ARG")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " NE ", info);
        atom[4] = CreateAAAtom(mol, " CZ ", info);
        atom[5] = CreateAAAtom(mol, " NH1", info);
        atom[6] = CreateAAAtom(mol, " NH2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        CreateAABond(mol, atom[3], atom[4], 1);
        CreateAABond(mol, atom[4], atom[5], 2);
        CreateAABond(mol, atom[4], atom[6], 1);
      } else if (!strcmp(aa, "ASP")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " OD1", info);
        atom[3] = CreateAAAtom(mol, " OD2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "ASN")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " OD1", info);
        atom[3] = CreateAAAtom(mol, " ND2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "ABA")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
      } else if (!strcmp(aa, "ACE")) {
        info.setIsHeteroAtom(true);
        r2 = CreateAAAtom(mol, " C  ", info);
        atom[0] = CreateAAAtom(mol, " O  ", info);
        atom[1] = CreateAAAtom(mol, " CH3", info);
        CreateAABond(mol, r2, atom[1], 1);
        CreateAABond(mol, r2, atom[0], 2);
      }
      break;

    case 'C':
      if (!strcmp(aa, "CYS")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        r3 = CreateAAAtom(mol, " SG ", info);
        CreateAABond(mol, atom[0], r3, 1);
      }
      break;

    case 'D':
      if (!strcmp(aa, "DAL")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
      } else if (!strcmp(aa, "DAR")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " NE ", info);
        atom[4] = CreateAAAtom(mol, " CZ ", info);
        atom[5] = CreateAAAtom(mol, " NH1", info);
        atom[6] = CreateAAAtom(mol, " NH2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        CreateAABond(mol, atom[3], atom[4], 1);
        CreateAABond(mol, atom[4], atom[5], 2);
        CreateAABond(mol, atom[4], atom[6], 1);
      } else if (!strcmp(aa, "DAS")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " OD1", info);
        atom[3] = CreateAAAtom(mol, " OD2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "DBB")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
      } else if (!strcmp(aa, "DCY")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        r3 = CreateAAAtom(mol, " SG ", info);
        CreateAABond(mol, atom[0], r3, 1);
      } else if (!strcmp(aa, "DGL")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " OE1", info);
        atom[4] = CreateAAAtom(mol, " OE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
      } else if (!strcmp(aa, "DGN")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " OE1", info);
        atom[4] = CreateAAAtom(mol, " NE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
      } else if (!strcmp(aa, "DHI")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " ND1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " NE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[4], 2);
        CreateAABond(mol, atom[4], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 1);
        CreateAABond(mol, atom[3], atom[1], 2);
      } else if (!strcmp(aa, "DIL")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        atom[3] = CreateAAAtom(mol, " CD1", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        atom[0]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      } else if (!strcmp(aa, "DLE")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "DLY")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        atom[4] = CreateAAAtom(mol, " NZ ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        CreateAABond(mol, atom[3], atom[4], 1);
      } else if (!strcmp(aa, "DPN")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CZ ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[4], atom[6], 2);
        CreateAABond(mol, atom[6], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 2);
        CreateAABond(mol, atom[3], atom[1], 1);
      } else if (!strcmp(aa, "DPR")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], r1, 1);
      } else if (!strcmp(aa, "DSG")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " OD1", info);
        atom[3] = CreateAAAtom(mol, " ND2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "DSN")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " OG ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
      } else if (!strcmp(aa, "DTH")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " OG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
        atom[0]->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      } else if (!strcmp(aa, "DTR")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " NE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CE3", info);
        atom[7] = CreateAAAtom(mol, " CZ2", info);
        atom[8] = CreateAAAtom(mol, " CZ3", info);
        atom[9] = CreateAAAtom(mol, " CH2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[3], atom[5], 2);
        CreateAABond(mol, atom[3], atom[6], 1);
        CreateAABond(mol, atom[4], atom[5], 1);
        CreateAABond(mol, atom[5], atom[7], 1);
        CreateAABond(mol, atom[6], atom[8], 2);
        CreateAABond(mol, atom[7], atom[9], 2);
        CreateAABond(mol, atom[8], atom[9], 1);
      } else if (!strcmp(aa, "DTY")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CZ ", info);
        atom[7] = CreateAAAtom(mol, " OH ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[4], atom[6], 2);
        CreateAABond(mol, atom[6], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 2);
        CreateAABond(mol, atom[3], atom[1], 1);
        CreateAABond(mol, atom[6], atom[7], 1);
      } else if (!strcmp(aa, "DVA")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
      }
      break;

    case 'G':
      if (!strcmp(aa, "GLN")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " OE1", info);
        atom[4] = CreateAAAtom(mol, " NE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
      } else if (!strcmp(aa, "GLU")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " OE1", info);
        atom[4] = CreateAAAtom(mol, " OE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
      } else if (!strcmp(aa, "GLY")) {
        r1 = CreateAAAtom(mol, " N  ", info);
        atom[0] = CreateAAAtom(mol, " CA ", info);
        r2 = CreateAAAtom(mol, " C  ", info);
        atom[1] = CreateAAAtom(mol, " O  ", info);
        CreateAABond(mol, r1, atom[0], 1);
        CreateAABond(mol, atom[0], r2, 1);
        CreateAABond(mol, r2, atom[1], 2);
      }
      break;

    case 'H':
      if (!strcmp(aa, "HIS")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " ND1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " NE2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[4], 2);
        CreateAABond(mol, atom[4], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 1);
        CreateAABond(mol, atom[3], atom[1], 2);
      }
      break;

    case 'I':
      if (!strcmp(aa, "ILE")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        atom[3] = CreateAAAtom(mol, " CD1", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        atom[0]->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      }
      break;

    case 'L':
      if (!strcmp(aa, "LEU")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[1], atom[3], 1);
      } else if (!strcmp(aa, "LYS")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        atom[4] = CreateAAAtom(mol, " NZ ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
        CreateAABond(mol, atom[3], atom[4], 1);
      }
      break;

    case 'M':
      if (!strcmp(aa, "MET")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " SD ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
      } else if (!strcmp(aa, "MED")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, -1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " SD ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
      } else if (!strcmp(aa, "MSE")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, "SE  ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
      }
      break;

    case 'N':
      if (!strcmp(aa, "NLE")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " CE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
      } else if (!strcmp(aa, "NVA")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
      }
      break;

    case 'O':
      if (!strcmp(aa, "ORN")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " NE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], atom[3], 1);
      }
      break;

    case 'P':
      if (!strcmp(aa, "PHE")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CZ ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[4], atom[6], 2);
        CreateAABond(mol, atom[6], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 2);
        CreateAABond(mol, atom[3], atom[1], 1);
      } else if (!strcmp(aa, "PRO")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], r1, 1);
      } else if (!strcmp(aa, "PCA")) {
        info.setIsHeteroAtom(true);
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD ", info);
        atom[3] = CreateAAAtom(mol, " OE ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 1);
        CreateAABond(mol, atom[2], r1, 1);
        CreateAABond(mol, atom[2], atom[3], 2);
      }
      break;

    case 'S':
      if (!strcmp(aa, "SER")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " OG ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
      } else if (!strcmp(aa, "SAR")) {
        info.setIsHeteroAtom(true);
        r1 = CreateAAAtom(mol, " N  ", info);
        atom[0] = CreateAAAtom(mol, " CA ", info);
        r2 = CreateAAAtom(mol, " C  ", info);
        atom[1] = CreateAAAtom(mol, " O  ", info);
        atom[2] = CreateAAAtom(mol, " CN ", info);
        CreateAABond(mol, r1, atom[0], 1);
        CreateAABond(mol, atom[0], r2, 1);
        CreateAABond(mol, r2, atom[1], 2);
        CreateAABond(mol, r1, atom[2], 1);
      }
      break;

    case 'T':
      if (!strcmp(aa, "THR")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " OG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
        atom[0]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      } else if (!strcmp(aa, "TRP")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " NE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CE3", info);
        atom[7] = CreateAAAtom(mol, " CZ2", info);
        atom[8] = CreateAAAtom(mol, " CZ3", info);
        atom[9] = CreateAAAtom(mol, " CH2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[1], atom[3], 1);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[3], atom[5], 2);
        CreateAABond(mol, atom[3], atom[6], 1);
        CreateAABond(mol, atom[4], atom[5], 1);
        CreateAABond(mol, atom[5], atom[7], 1);
        CreateAABond(mol, atom[6], atom[8], 2);
        CreateAABond(mol, atom[7], atom[9], 2);
        CreateAABond(mol, atom[8], atom[9], 1);
      } else if (!strcmp(aa, "TYR")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG ", info);
        atom[2] = CreateAAAtom(mol, " CD1", info);
        atom[3] = CreateAAAtom(mol, " CD2", info);
        atom[4] = CreateAAAtom(mol, " CE1", info);
        atom[5] = CreateAAAtom(mol, " CE2", info);
        atom[6] = CreateAAAtom(mol, " CZ ", info);
        atom[7] = CreateAAAtom(mol, " OH ", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[1], atom[2], 2);
        CreateAABond(mol, atom[2], atom[4], 1);
        CreateAABond(mol, atom[4], atom[6], 2);
        CreateAABond(mol, atom[6], atom[5], 1);
        CreateAABond(mol, atom[5], atom[3], 2);
        CreateAABond(mol, atom[3], atom[1], 1);
        CreateAABond(mol, atom[6], atom[7], 1);
      }
      break;

    case 'V':
      if (!strcmp(aa, "VAL")) {
        CreateAABackbone(mol, r1, r2, atom[0], info, 1);
        atom[1] = CreateAAAtom(mol, " CG1", info);
        atom[2] = CreateAAAtom(mol, " CG2", info);
        CreateAABond(mol, atom[0], atom[1], 1);
        CreateAABond(mol, atom[0], atom[2], 1);
      }
      break;
  }
}

static RWMol *AASequenceToMol(const char *seq, bool lowerD) {
  AtomPDBResidueInfo info;
  Atom *prev = (Atom *)nullptr;
  char chain[2];

  chain[0] = 'A';
  chain[1] = '\0';

  info.setSerialNumber(1);
  info.setAltLoc(" ");
  info.setResidueNumber(0);
  info.setInsertionCode(" ");
  info.setChainId(chain);
  auto *mol = new RWMol();

  while (*seq) {
    Atom *r1 = (Atom *)nullptr;
    Atom *r2 = (Atom *)nullptr;
    Atom *r3 = (Atom *)nullptr;

    switch (*seq) {
      case '\n':
      case '\r':
      case '-':
        seq++;
        continue;

      case ' ':
      case '\t':
      case '\0':
        break;

      case '.':
        if (prev) {
          Atom *oxt = CreateAAAtom(mol, " OXT", info);
          CreateAABond(mol, prev, oxt, 1);
          if (chain[0] < 'Z') {
            chain[0]++;
            info.setChainId(chain);
          }
          info.setResidueNumber(0);
          prev = (Atom *)nullptr;
        }
        seq++;
        continue;

      default:
        delete mol;
        return (RWMol *)nullptr;

      case 'A':
        CreateAminoAcid(mol, "ALA", r1, r2, r3, info);
        break;
      case 'C':
        CreateAminoAcid(mol, "CYS", r1, r2, r3, info);
        break;
      case 'D':
        CreateAminoAcid(mol, "ASP", r1, r2, r3, info);
        break;
      case 'E':
        CreateAminoAcid(mol, "GLU", r1, r2, r3, info);
        break;
      case 'F':
        CreateAminoAcid(mol, "PHE", r1, r2, r3, info);
        break;
      case 'G':
      case 'g':
        CreateAminoAcid(mol, "GLY", r1, r2, r3, info);
        break;
      case 'H':
        CreateAminoAcid(mol, "HIS", r1, r2, r3, info);
        break;
      case 'I':
        CreateAminoAcid(mol, "ILE", r1, r2, r3, info);
        break;
      case 'K':
        CreateAminoAcid(mol, "LYS", r1, r2, r3, info);
        break;
      case 'L':
        CreateAminoAcid(mol, "LEU", r1, r2, r3, info);
        break;
      case 'M':
        CreateAminoAcid(mol, "MET", r1, r2, r3, info);
        break;
      case 'N':
        CreateAminoAcid(mol, "ASN", r1, r2, r3, info);
        break;
      case 'P':
        CreateAminoAcid(mol, "PRO", r1, r2, r3, info);
        break;
      case 'Q':
        CreateAminoAcid(mol, "GLN", r1, r2, r3, info);
        break;
      case 'R':
        CreateAminoAcid(mol, "ARG", r1, r2, r3, info);
        break;
      case 'S':
        CreateAminoAcid(mol, "SER", r1, r2, r3, info);
        break;
      case 'T':
        CreateAminoAcid(mol, "THR", r1, r2, r3, info);
        break;
      case 'V':
        CreateAminoAcid(mol, "VAL", r1, r2, r3, info);
        break;
      case 'W':
        CreateAminoAcid(mol, "TRP", r1, r2, r3, info);
        break;
      case 'Y':
        CreateAminoAcid(mol, "TYR", r1, r2, r3, info);
        break;

      case 'a':
        CreateAminoAcid(mol, lowerD ? "DAL" : "ALA", r1, r2, r3, info);
        break;
      case 'c':
        CreateAminoAcid(mol, lowerD ? "DCY" : "CYS", r1, r2, r3, info);
        break;
      case 'd':
        CreateAminoAcid(mol, lowerD ? "DAS" : "ASP", r1, r2, r3, info);
        break;
      case 'e':
        CreateAminoAcid(mol, lowerD ? "DGL" : "GLU", r1, r2, r3, info);
        break;
      case 'f':
        CreateAminoAcid(mol, lowerD ? "DPN" : "PHE", r1, r2, r3, info);
        break;
      case 'h':
        CreateAminoAcid(mol, lowerD ? "DHI" : "HIS", r1, r2, r3, info);
        break;
      case 'i':
        CreateAminoAcid(mol, lowerD ? "DIL" : "ILE", r1, r2, r3, info);
        break;
      case 'k':
        CreateAminoAcid(mol, lowerD ? "DLY" : "LYS", r1, r2, r3, info);
        break;
      case 'l':
        CreateAminoAcid(mol, lowerD ? "DLE" : "LEU", r1, r2, r3, info);
        break;
      case 'm':
        CreateAminoAcid(mol, lowerD ? "MED" : "MET", r1, r2, r3, info);
        break;
      case 'n':
        CreateAminoAcid(mol, lowerD ? "DSG" : "ASN", r1, r2, r3, info);
        break;
      case 'p':
        CreateAminoAcid(mol, lowerD ? "DPR" : "PRO", r1, r2, r3, info);
        break;
      case 'q':
        CreateAminoAcid(mol, lowerD ? "DGN" : "GLN", r1, r2, r3, info);
        break;
      case 'r':
        CreateAminoAcid(mol, lowerD ? "DAR" : "ARG", r1, r2, r3, info);
        break;
      case 's':
        CreateAminoAcid(mol, lowerD ? "DSN" : "SER", r1, r2, r3, info);
        break;
      case 't':
        CreateAminoAcid(mol, lowerD ? "DTH" : "THR", r1, r2, r3, info);
        break;
      case 'v':
        CreateAminoAcid(mol, lowerD ? "DVA" : "VAL", r1, r2, r3, info);
        break;
      case 'w':
        CreateAminoAcid(mol, lowerD ? "DTR" : "TRP", r1, r2, r3, info);
        break;
      case 'y':
        CreateAminoAcid(mol, lowerD ? "DTY" : "TYR", r1, r2, r3, info);
        break;
    }
    if (prev && r1) {
      CreateAABond(mol, prev, r1, 1);
    }
    prev = r2;
    seq++;
  }

  if (prev) {
    Atom *oxt = CreateAAAtom(mol, " OXT", info);
    CreateAABond(mol, prev, oxt, 1);
  }
  return mol;
}

static void CreateNucleicAcid(RWMol *mol, const char *na, Atom *&r1, Atom *&r2,
                              AtomPDBResidueInfo &info, bool PCap5) {
  Atom *atom[32];

  r1 = (Atom *)nullptr;
  r2 = (Atom *)nullptr;

  int resno = info.getResidueNumber();
  info.setResidueNumber(resno + 1);
  info.setIsHeteroAtom(false);
  info.setResidueName(na);

  if (resno) {
    atom[1] = CreateAAAtom(mol, " P  ", info);
    atom[2] = CreateAAAtom(mol, " OP1", info);
    atom[3] = CreateAAAtom(mol, " OP2", info);
    atom[4] = CreateAAAtom(mol, " O5'", info);
    CreateAABond(mol, atom[1], atom[2], 2);
    CreateAABond(mol, atom[1], atom[3], 1);
    CreateAABond(mol, atom[1], atom[4], 1);
    r1 = atom[1];
  } else if (PCap5) {
    atom[0] = CreateAAAtom(mol, " OP3", info);
    atom[1] = CreateAAAtom(mol, " P  ", info);
    atom[2] = CreateAAAtom(mol, " OP1", info);
    atom[3] = CreateAAAtom(mol, " OP2", info);
    atom[4] = CreateAAAtom(mol, " O5'", info);
    CreateAABond(mol, atom[0], atom[1], 1);
    CreateAABond(mol, atom[1], atom[2], 2);
    CreateAABond(mol, atom[1], atom[3], 1);
    CreateAABond(mol, atom[1], atom[4], 1);
    r1 = atom[0];
  } else {
    atom[4] = CreateAAAtom(mol, " O5'", info);
    r1 = atom[4];
  }

  atom[5] = CreateAAAtom(mol, " C5'", info);
  atom[6] = CreateAAAtom(mol, " C4'", info);
  atom[7] = CreateAAAtom(mol, " O4'", info);
  atom[8] = CreateAAAtom(mol, " C3'", info);
  atom[9] = CreateAAAtom(mol, " O3'", info);
  atom[10] = CreateAAAtom(mol, " C2'", info);
  CreateAABond(mol, atom[4], atom[5], 1);
  CreateAABond(mol, atom[5], atom[6], 1);
  CreateAABond(mol, atom[6], atom[7], 1);
  CreateAABond(mol, atom[6], atom[8], 1);
  CreateAABond(mol, atom[8], atom[9], 1);
  CreateAABond(mol, atom[8], atom[10], 1);
  atom[6]->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  atom[8]->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  if (na[1] == ' ') {
    atom[11] = CreateAAAtom(mol, " O2'", info);
    atom[12] = CreateAAAtom(mol, " C1'", info);
    CreateAABond(mol, atom[10], atom[11], 1);
    CreateAABond(mol, atom[7], atom[12], 1);
    CreateAABond(mol, atom[10], atom[12], 1);
    atom[10]->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  } else {
    atom[12] = CreateAAAtom(mol, " C1'", info);
    CreateAABond(mol, atom[7], atom[12], 1);
    CreateAABond(mol, atom[10], atom[12], 1);
  }
  r2 = atom[9];

  switch (na[2]) {
    case 'A':
      atom[13] = CreateAAAtom(mol, " N9 ", info);
      atom[14] = CreateAAAtom(mol, " C8 ", info);
      atom[15] = CreateAAAtom(mol, " N7 ", info);
      atom[16] = CreateAAAtom(mol, " C5 ", info);
      atom[17] = CreateAAAtom(mol, " C6 ", info);
      atom[18] = CreateAAAtom(mol, " N6 ", info);
      atom[19] = CreateAAAtom(mol, " N1 ", info);
      atom[20] = CreateAAAtom(mol, " C2 ", info);
      atom[21] = CreateAAAtom(mol, " N3 ", info);
      atom[22] = CreateAAAtom(mol, " C4 ", info);
      CreateAABond(mol, atom[12], atom[13], 1);
      CreateAABond(mol, atom[13], atom[14], 1);
      CreateAABond(mol, atom[13], atom[22], 1);
      CreateAABond(mol, atom[14], atom[15], 2);
      CreateAABond(mol, atom[15], atom[16], 1);
      CreateAABond(mol, atom[16], atom[17], 1);
      CreateAABond(mol, atom[16], atom[22], 2);
      CreateAABond(mol, atom[17], atom[18], 1);
      CreateAABond(mol, atom[17], atom[19], 2);
      CreateAABond(mol, atom[19], atom[20], 1);
      CreateAABond(mol, atom[20], atom[21], 2);
      CreateAABond(mol, atom[21], atom[22], 1);
      atom[12]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
    case 'C':
      atom[13] = CreateAAAtom(mol, " N1 ", info);
      atom[14] = CreateAAAtom(mol, " C2 ", info);
      atom[15] = CreateAAAtom(mol, " O2 ", info);
      atom[16] = CreateAAAtom(mol, " N3 ", info);
      atom[17] = CreateAAAtom(mol, " C4 ", info);
      atom[18] = CreateAAAtom(mol, " N4 ", info);
      atom[19] = CreateAAAtom(mol, " C5 ", info);
      atom[20] = CreateAAAtom(mol, " C6 ", info);
      CreateAABond(mol, atom[12], atom[13], 1);
      CreateAABond(mol, atom[13], atom[14], 1);
      CreateAABond(mol, atom[13], atom[20], 1);
      CreateAABond(mol, atom[14], atom[15], 2);
      CreateAABond(mol, atom[14], atom[16], 1);
      CreateAABond(mol, atom[16], atom[17], 2);
      CreateAABond(mol, atom[17], atom[18], 1);
      CreateAABond(mol, atom[17], atom[19], 1);
      CreateAABond(mol, atom[19], atom[20], 2);
      atom[12]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
    case 'G':
      atom[13] = CreateAAAtom(mol, " N9 ", info);
      atom[14] = CreateAAAtom(mol, " C8 ", info);
      atom[15] = CreateAAAtom(mol, " N7 ", info);
      atom[16] = CreateAAAtom(mol, " C5 ", info);
      atom[17] = CreateAAAtom(mol, " C6 ", info);
      atom[18] = CreateAAAtom(mol, " O6 ", info);
      atom[19] = CreateAAAtom(mol, " N1 ", info);
      atom[20] = CreateAAAtom(mol, " C2 ", info);
      atom[21] = CreateAAAtom(mol, " N2 ", info);
      atom[22] = CreateAAAtom(mol, " N3 ", info);
      atom[23] = CreateAAAtom(mol, " C4 ", info);
      CreateAABond(mol, atom[12], atom[13], 1);
      CreateAABond(mol, atom[13], atom[14], 1);
      CreateAABond(mol, atom[13], atom[23], 1);
      CreateAABond(mol, atom[14], atom[15], 2);
      CreateAABond(mol, atom[15], atom[16], 1);
      CreateAABond(mol, atom[16], atom[17], 1);
      CreateAABond(mol, atom[16], atom[23], 2);
      CreateAABond(mol, atom[17], atom[18], 2);
      CreateAABond(mol, atom[17], atom[19], 1);
      CreateAABond(mol, atom[19], atom[20], 1);
      CreateAABond(mol, atom[20], atom[21], 1);
      CreateAABond(mol, atom[20], atom[22], 2);
      CreateAABond(mol, atom[22], atom[23], 1);
      atom[12]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
    case 'T':
      atom[13] = CreateAAAtom(mol, " N1 ", info);
      atom[14] = CreateAAAtom(mol, " C2 ", info);
      atom[15] = CreateAAAtom(mol, " O2 ", info);
      atom[16] = CreateAAAtom(mol, " N3 ", info);
      atom[17] = CreateAAAtom(mol, " C4 ", info);
      atom[18] = CreateAAAtom(mol, " O4 ", info);
      atom[19] = CreateAAAtom(mol, " C5 ", info);
      atom[20] = CreateAAAtom(mol, " C7 ", info);
      atom[21] = CreateAAAtom(mol, " C6 ", info);
      CreateAABond(mol, atom[12], atom[13], 1);
      CreateAABond(mol, atom[13], atom[14], 1);
      CreateAABond(mol, atom[13], atom[21], 1);
      CreateAABond(mol, atom[14], atom[15], 2);
      CreateAABond(mol, atom[14], atom[16], 1);
      CreateAABond(mol, atom[16], atom[17], 1);
      CreateAABond(mol, atom[17], atom[18], 2);
      CreateAABond(mol, atom[17], atom[19], 1);
      CreateAABond(mol, atom[19], atom[20], 1);
      CreateAABond(mol, atom[19], atom[21], 2);
      atom[12]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
    case 'U':
      atom[13] = CreateAAAtom(mol, " N1 ", info);
      atom[14] = CreateAAAtom(mol, " C2 ", info);
      atom[15] = CreateAAAtom(mol, " O2 ", info);
      atom[16] = CreateAAAtom(mol, " N3 ", info);
      atom[17] = CreateAAAtom(mol, " C4 ", info);
      atom[18] = CreateAAAtom(mol, " O4 ", info);
      atom[19] = CreateAAAtom(mol, " C5 ", info);
      atom[20] = CreateAAAtom(mol, " C6 ", info);
      CreateAABond(mol, atom[12], atom[13], 1);
      CreateAABond(mol, atom[13], atom[14], 1);
      CreateAABond(mol, atom[13], atom[20], 1);
      CreateAABond(mol, atom[14], atom[15], 2);
      CreateAABond(mol, atom[14], atom[16], 1);
      CreateAABond(mol, atom[16], atom[17], 1);
      CreateAABond(mol, atom[17], atom[18], 2);
      CreateAABond(mol, atom[17], atom[19], 1);
      CreateAABond(mol, atom[19], atom[20], 2);
      atom[12]->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
  }
}

static void CreatePCap3(RWMol *mol, Atom *prev, AtomPDBResidueInfo &info) {
  Atom *atom[4];

  int resno = info.getResidueNumber();
  info.setResidueNumber(resno + 1);
  info.setIsHeteroAtom(true);
  info.setResidueName("2PO");

  atom[0] = CreateAAAtom(mol, " P  ", info);
  atom[1] = CreateAAAtom(mol, " OP1", info);
  atom[2] = CreateAAAtom(mol, " OP2", info);
  atom[3] = CreateAAAtom(mol, " OP3", info);
  CreateAABond(mol, prev, atom[0], 1);
  CreateAABond(mol, atom[0], atom[1], 2);
  CreateAABond(mol, atom[0], atom[2], 1);
  CreateAABond(mol, atom[0], atom[3], 1);
}

static RWMol *NASequenceToMol(const char *seq, bool Dna, bool PCap5,
                              bool PCap3) {
  char chain[2];
  chain[0] = 'A';
  chain[1] = '\0';

  Atom *prev = (Atom *)nullptr;
  AtomPDBResidueInfo info;
  info.setSerialNumber(1);
  info.setAltLoc(" ");
  info.setResidueNumber(0);
  info.setInsertionCode(" ");
  info.setChainId(chain);
  auto *mol = new RWMol();

  while (*seq) {
    Atom *r1 = (Atom *)nullptr;
    Atom *r2 = (Atom *)nullptr;

    switch (*seq) {
      case '\n':
      case '\r':
      case '-':
        seq++;
        continue;

      case ' ':
      case '\t':
      case '\0':
        break;

      case '.':
        if (prev) {
          if (PCap3) {
            CreatePCap3(mol, prev, info);
          }
          if (chain[0] < 'Z') {
            chain[0]++;
            info.setChainId(chain);
          }
          info.setResidueNumber(0);
        }
        prev = (Atom *)nullptr;
        seq++;
        continue;

      default:
        delete mol;
        return (RWMol *)nullptr;

      case 'A':
      case 'a':
        CreateNucleicAcid(mol, Dna ? " DA" : "  A", r1, r2, info, PCap5);
        break;
      case 'C':
      case 'c':
        CreateNucleicAcid(mol, Dna ? " DC" : "  C", r1, r2, info, PCap5);
        break;
      case 'G':
      case 'g':
        CreateNucleicAcid(mol, Dna ? " DG" : "  G", r1, r2, info, PCap5);
        break;
      case 'T':
      case 't':
        CreateNucleicAcid(mol, Dna ? " DT" : "  T", r1, r2, info, PCap5);
        break;
      case 'U':
      case 'u':
        CreateNucleicAcid(mol, Dna ? " DU" : "  U", r1, r2, info, PCap5);
        break;
    }
    if (prev && r1) {
      CreateAABond(mol, prev, r1, 1);
    }
    prev = r2;
    seq++;
  }
  if (prev && PCap3) {
    CreatePCap3(mol, prev, info);
  }
  return mol;
}

RWMol *SequenceToMol(const char *seq, bool sanitize, int flavor) {
  if (!seq) {
    return (RWMol *)nullptr;
  }
  RWMol *mol;

  switch (flavor) {
    /* Protein */
    case 0:
      mol = AASequenceToMol(seq, false);
      break;
    case 1:
      mol = AASequenceToMol(seq, true);
      break;

    /* RNA */
    case 2:
      mol = NASequenceToMol(seq, false, false, false);
      break;
    case 3:
      mol = NASequenceToMol(seq, false, true, false);
      break;
    case 4:
      mol = NASequenceToMol(seq, false, false, true);
      break;
    case 5:
      mol = NASequenceToMol(seq, false, true, true);
      break;

    /* DNA */
    case 6:
      mol = NASequenceToMol(seq, true, false, false);
      break;
    case 7:
      mol = NASequenceToMol(seq, true, true, false);
      break;
    case 8:
      mol = NASequenceToMol(seq, true, false, true);
      break;
    case 9:
      mol = NASequenceToMol(seq, true, true, true);
      break;

    default:
      return (RWMol *)nullptr;
  }
  if (sanitize && mol) {
    MolOps::sanitizeMol(*mol);
  }
  return mol;
}

RWMol *SequenceToMol(const std::string &seq, bool sanitize, int flavor) {
  return SequenceToMol(seq.c_str(), sanitize, flavor);
}

RWMol *SequenceToMol(const char *seq, bool sanitize, bool lowerD) {
  return SequenceToMol(seq, sanitize, lowerD ? 1 : 0);
}

RWMol *SequenceToMol(const std::string &seq, bool sanitize, bool lowerD) {
  return SequenceToMol(seq.c_str(), sanitize, lowerD ? 1 : 0);
}

RWMol *FASTAToMol(const char *seq, bool sanitize, int flavor) {
  if (!seq) {
    return (RWMol *)nullptr;
  }

  std::string title;
  if (seq[0] == '>') {
    seq++;
    while (*seq && *seq != '\n' && *seq != '\r') {
      title += *seq++;
    }
  }
  RWMol *mol = SequenceToMol(seq, sanitize, flavor);
  if (!title.empty()) {
    mol->setProp(common_properties::_Name, title);
  }
  return mol;
}

RWMol *FASTAToMol(const std::string &seq, bool sanitize, int flavor) {
  return FASTAToMol(seq.c_str(), sanitize, flavor);
}

RWMol *FASTAToMol(const char *seq, bool sanitize, bool lowerD) {
  return FASTAToMol(seq, sanitize, lowerD ? 1 : 0);
}

RWMol *FASTAToMol(const std::string &seq, bool sanitize, bool lowerD) {
  return FASTAToMol(seq.c_str(), sanitize, lowerD ? 1 : 0);
}

struct HELMMonomer {
  Atom *r1{nullptr};
  Atom *r2{nullptr};
  Atom *r3{nullptr};
  Atom *oxt{nullptr};

  HELMMonomer() {}
  HELMMonomer(Atom *x, Atom *y, Atom *z) : r1(x), r2(y), r3(z), oxt(nullptr) {}
};

static const char *GetHELMOneLetterCode(char ch) {
  switch (ch) {
    case 'A':
      return "ALA";
    case 'C':
      return "CYS";
    case 'D':
      return "ASP";
    case 'E':
      return "GLU";
    case 'F':
      return "PHE";
    case 'G':
      return "GLY";
    case 'H':
      return "HIS";
    case 'I':
      return "ILE";
    case 'K':
      return "LYS";
    case 'L':
      return "LEU";
    case 'M':
      return "MET";
    case 'N':
      return "ASN";
    case 'P':
      return "PRO";
    case 'Q':
      return "GLN";
    case 'R':
      return "ARG";
    case 'S':
      return "SER";
    case 'T':
      return "THR";
    case 'V':
      return "VAL";
    case 'W':
      return "TRP";
    case 'Y':
      return "TYR";
  }
  return (char *)nullptr;
}

static bool IsHELMMonomerIDChar(char ch) {
  if (ch >= 'A' && ch <= 'Z') {
    return true;
  }
  if (ch >= 'a' && ch <= 'z') {
    return true;
  }
  if (ch >= '0' && ch <= '9') {
    return true;
  }
  return false;
}

static const char *LookupHELMPeptideMonomer(const char *ptr) {
  switch (ptr[0]) {
    case 'A':
      if (ptr[1] == '\0') {
        return "ALA";
      }
      if (ptr[1] == 'b' && ptr[2] == 'u' && ptr[3] == '\0') {
        return "ABA";
      }
      break;
    case 'C':
      if (ptr[1] == '\0') {
        return "CYS";
      }
      break;
    case 'D':
      if (ptr[1] == '\0') {
        return "ASP";
      }
      break;
    case 'E':
      if (ptr[1] == '\0') {
        return "GLU";
      }
      break;
    case 'F':
      if (ptr[1] == '\0') {
        return "PHE";
      }
      break;
    case 'G':
      if (ptr[1] == '\0') {
        return "GLY";
      }
      if (ptr[1] == 'l' && ptr[2] == 'p' && ptr[3] == '\0') {
        return "PCA";
      }
      break;
    case 'H':
      if (ptr[1] == '\0') {
        return "HIS";
      }
      break;
    case 'I':
      if (ptr[1] == '\0') {
        return "ILE";
      }
      break;
    case 'K':
      if (ptr[1] == '\0') {
        return "LYS";
      }
      break;
    case 'L':
      if (ptr[1] == '\0') {
        return "LEU";
      }
      break;
    case 'M':
      if (ptr[1] == '\0') {
        return "MET";
      }
      break;
    case 'N':
      if (ptr[1] == '\0') {
        return "ASN";
      }
      if (ptr[1] == 'a' && ptr[2] == 'l' && ptr[3] == '\0') {
        return "NAL";
      }
      if (ptr[1] == 'l' && ptr[2] == 'e' && ptr[3] == '\0') {
        return "NLE";
      }
      if (ptr[1] == 'v' && ptr[2] == 'a' && ptr[3] == '\0') {
        return "NVA";
      }
      break;
    case 'O':
      if (ptr[1] == 'r' && ptr[2] == 'n' && ptr[3] == '\0') {
        return "ORN";
      }
      break;
    case 'P':
      if (ptr[1] == '\0') {
        return "PRO";
      }
      break;
    case 'Q':
      if (ptr[1] == '\0') {
        return "GLN";
      }
      break;
    case 'R':
      if (ptr[1] == '\0') {
        return "ARG";
      }
      break;
    case 'S':
      if (ptr[1] == '\0') {
        return "SER";
      }
      if (ptr[1] == 'a' && ptr[2] == 'r' && ptr[3] == '\0') {
        return "SAR";
      }
      break;
    case 'T':
      if (ptr[1] == '\0') {
        return "THR";
      }
      break;
    case 'V':
      if (ptr[1] == '\0') {
        return "VAL";
      }
      break;
    case 'W':
      if (ptr[1] == '\0') {
        return "TRP";
      }
      break;
    case 'Y':
      if (ptr[1] == '\0') {
        return "TYR";
      }
      break;
    case 'd':
      switch (ptr[1]) {
        case 'A':
          if (ptr[2] == '\0') {
            return "DAL";
          }
          break;
        case 'C':
          if (ptr[2] == '\0') {
            return "DCY";
          }
          break;
        case 'D':
          if (ptr[2] == '\0') {
            return "DAS";
          }
          break;
        case 'E':
          if (ptr[2] == '\0') {
            return "DGL";
          }
          break;
        case 'F':
          if (ptr[2] == '\0') {
            return "DPN";
          }
          break;
        case 'H':
          if (ptr[2] == '\0') {
            return "DHI";
          }
          break;
        case 'I':
          if (ptr[2] == '\0') {
            return "DIL";
          }
          break;
        case 'K':
          if (ptr[2] == '\0') {
            return "DLY";
          }
          break;
        case 'L':
          if (ptr[2] == '\0') {
            return "DLE";
          }
          break;
        case 'M':
          if (ptr[2] == '\0') {
            return "MED";
          }
          break;
        case 'N':
          if (ptr[2] == '\0') {
            return "DSG";
          }
          break;
        case 'P':
          if (ptr[2] == '\0') {
            return "DPR";
          }
          break;
        case 'Q':
          if (ptr[2] == '\0') {
            return "DGN";
          }
          break;
        case 'R':
          if (ptr[2] == '\0') {
            return "DAR";
          }
          break;
        case 'S':
          if (ptr[2] == '\0') {
            return "DSN";
          }
          break;
        case 'T':
          if (ptr[2] == '\0') {
            return "DTH";
          }
          break;
        case 'V':
          if (ptr[2] == '\0') {
            return "DVA";
          }
          break;
        case 'W':
          if (ptr[2] == '\0') {
            return "DTR";
          }
          break;
        case 'Y':
          if (ptr[2] == '\0') {
            return "DTY";
          }
          break;
      }
      break;
    case 's':
      if (ptr[1] == 'e' && ptr[2] == 'C' && ptr[3] == '\0') {
        return "MSE";
      }
      break;
  }
  return (const char *)nullptr;
}

static const char *ParseHELMPeptide(RWMol *mol, const char *ptr,
                                    const char *chain,
                                    std::vector<HELMMonomer> &vseq) {
  unsigned int len = 0;
  HELMMonomer curr;

  vseq.clear();
  if (ptr[0] == '}') {
    return ptr;
  }

  AtomPDBResidueInfo info;
  info.setSerialNumber(1);
  info.setAltLoc(" ");
  info.setResidueNumber(0);
  info.setInsertionCode(" ");
  info.setChainId(chain);

  if (ptr[0] == '[' && ptr[1] == 'a' && ptr[2] == 'c' && ptr[3] == ']') {
    if (ptr[4] != '.') {
      return (const char *)nullptr;
    }
    info.setResidueNumber(-2);
    CreateAminoAcid(mol, "ACE", curr.r1, curr.r2, curr.r3, info);
    vseq.push_back(curr);
    info.setResidueNumber(0);
    ptr += 5;
    len = 1;
  }

  for (;;) {
    const char *name = nullptr;
    if (*ptr == '[') {
      std::string tmp;
      ptr++;
      while (IsHELMMonomerIDChar(*ptr)) {
        tmp += *ptr++;
      }
      if (*ptr != ']') {
        return (char *)nullptr;
      }
      name = LookupHELMPeptideMonomer(tmp.c_str());
    } else {
      name = GetHELMOneLetterCode(*ptr);
    }
    if (!name) {
      return (const char *)nullptr;
    }
    ptr++;

    CreateAminoAcid(mol, name, curr.r1, curr.r2, curr.r3, info);
    if (len && vseq[len - 1].r2 && curr.r1) {
      CreateAABond(mol, vseq[len - 1].r2, curr.r1, 1);
      vseq[len - 1].r2 = nullptr;
    }
    vseq.push_back(curr);
    len++;

    if (*ptr == '.') {
      if (ptr[1] == '[' && ptr[2] == 'a' && ptr[3] == 'm' && ptr[4] == ']' &&
          ptr[5] == '}') {
        if (!vseq[len - 1].r2) {
          return (const char *)nullptr;
        }
        int resno = info.getResidueNumber();
        info.setResidueNumber(resno + 1);
        info.setIsHeteroAtom(true);
        info.setResidueName("NH2");
        Atom *n = CreateAAAtom(mol, " N  ", info);
        CreateAABond(mol, vseq[len - 1].r2, n, 1);
        vseq[len - 1].r2 = (Atom *)nullptr;
        vseq.emplace_back();
        len++;
        return ptr + 5;
      }
      ptr++;
    } else if (*ptr == '}') {
      if (!vseq[len - 1].r2) {
        return (const char *)nullptr;
      }
      Atom *oxt = CreateAAAtom(mol, " OXT", info);
      CreateAABond(mol, vseq[len - 1].r2, oxt, 1);
      vseq[len - 1].oxt = oxt;
      return ptr;
    } else {
      return (const char *)nullptr;
    }
  }
}

static const char *ParseHELMNucleic(RWMol *mol, const char *ptr,
                                    const char *chain) {
  if (ptr[0] == '}') {
    return ptr;
  }

  bool PCap5 = false;
  Atom *prev = nullptr;
  Atom *r1 = nullptr;
  Atom *r2 = nullptr;

  AtomPDBResidueInfo info;
  info.setSerialNumber(1);
  info.setAltLoc(" ");
  info.setResidueNumber(0);
  info.setInsertionCode(" ");
  info.setChainId(chain);

  if (*ptr == 'P') {
    PCap5 = true;
    ptr++;
    if (*ptr == '.') {
      ptr++;
    }
  }

  for (;;) {
    const char *name = nullptr;
    if (*ptr == 'R' && ptr[1] == '(') {
      if (ptr[2] == 'A') {
        if (ptr[3] == ')') {
          name = "  A";
          ptr += 4;
        }
      } else if (ptr[2] == 'C') {
        if (ptr[3] == ')') {
          name = "  C";
          ptr += 4;
        }
      } else if (ptr[2] == 'G') {
        if (ptr[3] == ')') {
          name = "  G";
          ptr += 4;
        }
      } else if (ptr[2] == 'T') {
        if (ptr[3] == ')') {
          name = "  T";
          ptr += 4;
        }
      } else if (ptr[2] == 'U') {
        if (ptr[3] == ')') {
          name = "  U";
          ptr += 4;
        }
      }
    } else if (*ptr == '[' && ptr[1] == 'd' && ptr[2] == 'R' && ptr[3] == ']' &&
               ptr[4] == '(') {
      if (ptr[5] == 'A') {
        if (ptr[6] == ')') {
          name = " DA";
          ptr += 7;
        }
      } else if (ptr[5] == 'C') {
        if (ptr[6] == ')') {
          name = " DC";
          ptr += 7;
        }
      } else if (ptr[5] == 'G') {
        if (ptr[6] == ')') {
          name = " DG";
          ptr += 7;
        }
      } else if (ptr[5] == 'T') {
        if (ptr[6] == ')') {
          name = " DT";
          ptr += 7;
        }
      } else if (ptr[5] == 'U') {
        if (ptr[6] == ')') {
          name = " DU";
          ptr += 7;
        }
      }
    }
    if (!name) {
      return (const char *)nullptr;
    }

    CreateNucleicAcid(mol, name, r1, r2, info, PCap5);
    if (prev && r1) {
      CreateAABond(mol, prev, r1, 1);
    }
    prev = r2;

    if (*ptr == '}') {
      return ptr;
    }
    if (*ptr == '.') {
      ptr++;
    }
    if (*ptr != 'P') {
      return (const char *)nullptr;
    }
    ptr++;
    if (*ptr == '}') {
      CreatePCap3(mol, prev, info);
      return ptr;
    }
    PCap5 = true;
    if (*ptr == '.') {
      ptr++;
    }
  }
}

static bool ParseHELM(RWMol *mol, const char *ptr) {
  std::map<std::string, std::vector<HELMMonomer>> seqs;
  const char *orig;
  char chain[2];
  chain[0] = 'A';
  chain[1] = '\0';

  for (;;) {
    orig = ptr;
    if (ptr[0] == 'P' && ptr[1] == 'E' && ptr[2] == 'P' && ptr[3] == 'T' &&
        ptr[4] == 'I' && ptr[5] == 'D' && ptr[6] == 'E' && ptr[7] >= '1' &&
        ptr[7] <= '9') {
      ptr += 8;
      while (*ptr >= '0' && *ptr <= '9') {
        ptr++;
      }
      if (*ptr != '{') {
        return false;
      }
      std::string id(orig, ptr - orig);
      chain[0] = 'A' + (orig[7] - '1');
      ptr = ParseHELMPeptide(mol, ptr + 1, chain, seqs[id]);
      if (!ptr || *ptr != '}') {
        return false;
      }
      ptr++;
    } else if (ptr[0] == 'R' && ptr[1] == 'N' && ptr[2] == 'A' &&
               ptr[3] >= '1' && ptr[3] <= '9') {
      ptr += 4;
      while (*ptr >= '0' && *ptr <= '9') {
        ptr++;
      }
      if (*ptr != '{') {
        return false;
      }
      chain[0] = 'A' + (orig[3] - '1');
      ptr = ParseHELMNucleic(mol, ptr + 1, chain);
      if (!ptr || *ptr != '}') {
        return false;
      }
      ptr++;
    } else {
      return false;
    }

    if (*ptr == '$') {
      break;
    }
    if (*ptr == '\0') {
      return true;
    }
    if (*ptr != '|') {
      return false;
    }
    ptr++;
  }
  ptr++;

  if (ptr[0] == '$' && ptr[1] == '$' && ptr[2] == '$') {
    return true;
  }

  for (;;) {
    orig = ptr;
    if (ptr[0] == 'P' && ptr[1] == 'E' && ptr[2] == 'P' && ptr[3] == 'T' &&
        ptr[4] == 'I' && ptr[5] == 'D' && ptr[6] == 'E' && ptr[7] >= '1' &&
        ptr[7] <= '9') {
      ptr += 8;
    } else {
      return false;
    }
    while (*ptr >= '0' && *ptr <= '9') {
      ptr++;
    }
    if (*ptr != ',') {
      return false;
    }

    std::string id1(orig, ptr - orig);
    ptr++;

    orig = ptr;
    if (ptr[0] == 'P' && ptr[1] == 'E' && ptr[2] == 'P' && ptr[3] == 'T' &&
        ptr[4] == 'I' && ptr[5] == 'D' && ptr[6] == 'E' && ptr[7] >= '1' &&
        ptr[7] <= '9') {
      ptr += 8;
    } else {
      return false;
    }
    while (*ptr >= '0' && *ptr <= '9') {
      ptr++;
    }
    if (*ptr != ',') {
      return false;
    }

    std::string id2(orig, ptr - orig);
    ptr++;

    unsigned int res1;
    unsigned int res2;
    unsigned int res1r;
    unsigned int res2r;

    if (*ptr >= '1' && *ptr <= '9') {
      res1 = (*ptr++) - '0';
      while (*ptr >= '0' && *ptr <= '9') {
        res1 = 10 * res1 + ((*ptr++) - '0');
      }
    } else {
      return false;
    }
    if (ptr[0] == ':' && ptr[1] == 'R' && ptr[2] >= '1' && ptr[2] <= '9') {
      res1r = ptr[2] - '0';
      ptr += 3;
    } else {
      return false;
    }
    if (*ptr != '-') {
      return false;
    }
    ptr++;

    if (*ptr >= '1' && *ptr <= '9') {
      res2 = (*ptr++) - '0';
      while (*ptr >= '0' && *ptr <= '9') {
        res2 = 10 * res2 + ((*ptr++) - '0');
      }
    } else {
      return false;
    }
    if (ptr[0] == ':' && ptr[1] == 'R' && ptr[2] >= '1' && ptr[2] <= '9') {
      res2r = ptr[2] - '0';
      ptr += 3;
    } else {
      return false;
    }

    // printf("%s:%u:R%u - %s:%u:R%u\n",id1.c_str(),res1,res1r,
    //                                  id2.c_str(),res2,res2r);

    if (res1 < 1 || res2 < 1) {
      return false;
    }
    if (seqs.find(id1) == seqs.end() || seqs.find(id2) == seqs.end()) {
      return false;
    }
    std::vector<HELMMonomer> *vseq1 = &seqs[id1];
    if (res1 > (unsigned int)vseq1->size()) {
      return false;
    }
    std::vector<HELMMonomer> *vseq2 = &seqs[id2];
    if (res2 > (unsigned int)vseq2->size()) {
      return false;
    }

    if (res1r == 3 && res2r == 3) {
      Atom *src = (*vseq1)[res1 - 1].r3;
      Atom *dst = (*vseq2)[res2 - 1].r3;
      if (src && dst && src != dst) {
        CreateAABond(mol, src, dst, 1);
        (*vseq1)[res1 - 1].r3 = (Atom *)nullptr;
        (*vseq2)[res2 - 1].r3 = (Atom *)nullptr;
      } else {
        return false;
      }
    } else if (res1r == 1 && res2r == 2) {
      Atom *src = (*vseq1)[res1 - 1].r1;
      Atom *dst = (*vseq2)[res2 - 1].r2;
      Atom *oxt = (*vseq2)[res2 - 1].oxt;
      if (src && dst && oxt && src != dst) {
        mol->removeAtom(oxt);
        CreateAABond(mol, src, dst, 1);
        (*vseq1)[res1 - 1].r1 = (Atom *)nullptr;
        (*vseq2)[res2 - 1].r2 = (Atom *)nullptr;
      } else {
        return false;
      }
    } else if (res1r == 2 && res2r == 1) {
      Atom *src = (*vseq2)[res2 - 1].r1;
      Atom *dst = (*vseq1)[res1 - 1].r2;
      Atom *oxt = (*vseq1)[res1 - 1].oxt;
      if (src && dst && oxt && src != dst) {
        mol->removeAtom(oxt);
        CreateAABond(mol, dst, src, 1);
        (*vseq1)[res1 - 1].r2 = (Atom *)nullptr;
        (*vseq2)[res2 - 1].r1 = (Atom *)nullptr;
      } else {
        return false;
      }
    } else {
      return false;
    }

    if (*ptr == '$') {
      break;
    }
    if (*ptr != '|') {
      return false;
    }
    ptr++;
  }
  ptr++;
  return ptr[0] == '$' && ptr[1] == '$';
}

RWMol *HELMToMol(const char *helm, bool sanitize) {
  auto *mol = new RWMol();

  const char *ptr = helm;
  if (ptr[0] == '$' && ptr[1] == '$' && ptr[2] == '$' && ptr[3] == '$') {
    return mol;
  }

  if (ParseHELM(mol, ptr)) {
    if (sanitize) {
      MolOps::sanitizeMol(*mol);
    }
    return mol;
  }
  delete mol;
  return (RWMol *)nullptr;
}

RWMol *HELMToMol(const std::string &helm, bool sanitize) {
  return HELMToMol(helm.c_str(), sanitize);
}

}  // namespace RDKit
