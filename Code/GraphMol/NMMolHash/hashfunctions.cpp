/*==============================================*/
/* Copyright (C)  2011-2019  NextMove Software  */
/* All rights reserved.                         */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#define _CRT_SECURE_NO_WARNINGS

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <string>

#include "toolkit.h"
#include "molhash.h"
#include "mf.h"


static unsigned int NMDetermineComponents(NMS_pMOL mol, unsigned int *parts,
  unsigned int acount)
{
  memset(parts, 0, acount * sizeof(unsigned int));
  std::vector<NMS_pATOM> todo;

  unsigned int result = 0;
  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    unsigned int idx = NMS_ATOM_GET_IDX(aptr);
    if (parts[idx] == 0) {
      parts[idx] = ++result;
      todo.push_back(aptr);

      while (!todo.empty()) {
        aptr = todo.back();
        todo.pop_back();
        NMS_FOR_NBR_OF_ATOM(nbor, aptr) {
          NMS_pATOM nptr = NMS_ITER_ATOM_NBR(nbor, aptr);
          idx = NMS_ATOM_GET_IDX(nptr);
          if (parts[idx] == 0) {
            parts[idx] = result;
            todo.push_back(nptr);
          }
        }
      }
    }
  }
  return result;
}


static std::string NMMolecularFormula(NMS_pMOL mol,
                                      const unsigned int *parts,
                                      unsigned int part)
{
  unsigned int hist[256];
  int charge = 0;

  memset(hist, 0, sizeof(hist));

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    unsigned int idx = NMS_ATOM_GET_IDX(aptr);
    if (part == 0 || parts[idx] == part) {
      unsigned int elem = NMS_ATOM_GET_ATOMICNUM(aptr);
      if (elem < 256)
        hist[elem]++;
      else hist[0]++;
      hist[1] += NMS_ATOM_GET_IMPLICITHCOUNT(aptr);
      charge += NMS_ATOM_GET_FORMALCHARGE(aptr);
    }
  }

  char buffer[16];
  std::string result;
  const unsigned char *perm = hist[6] ? OrganicHillOrder : InorganicHillOrder;
  for (unsigned int i = 0; i < 119; i++) {
    unsigned int elem = perm[i];
    if (hist[elem] > 0) {
      result += symbol[elem];
      if (hist[elem] > 1) {
        sprintf(buffer, "%u", hist[elem]);
        result += buffer;
      }
    }
  }
  if (charge != 0) {
    if (charge > 0) {
      result += '+';
      if (charge > 1) {
        sprintf(buffer, "%d", charge);
        result += buffer;
      }
    }
    else { // charge < 0
      result += '-';
      if (charge < -1) {
        sprintf(buffer, "%d", -charge);
        result += buffer;
      }
    }
  }
  return result;
}


static std::string NMMolecularFormula(NMS_pMOL mol, bool sep=false)
{
  if (!sep)
    return NMMolecularFormula(mol, 0, 0);

  unsigned int acount = NMS_MOL_GET_MAXATOMIDX(mol);
  if (acount == 0)
    return "";

  unsigned int size = (unsigned int)(acount * sizeof(int));
  unsigned int *parts = (unsigned int*)malloc(size);
  unsigned int pcount = NMDetermineComponents(mol, parts, acount);

  std::string result;
  if (pcount > 1) {
    std::vector<std::string> vmf;
    for (unsigned int i = 1; i <= pcount; i++)
      vmf.push_back(NMMolecularFormula(mol, parts, i));

    // sort
    result = vmf[0];
    for (unsigned int i = 1; i < pcount; i++) {
      result += ".";
      result += vmf[i];
    }
  }
  else // pcount == 1
    result = NMMolecularFormula(mol, parts, 1);
  free(parts);
  return result;
}

static void NormalizeHCount(NMS_pATOM aptr)
{
  unsigned int hcount;

  switch (NMS_ATOM_GET_ATOMICNUM(aptr)) {
  case 9:  // Fluorine
  case 17: // Chlorine
  case 35: // Bromine
  case 53: // Iodine
    hcount = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
    hcount = hcount < 1 ? 1 - hcount : 0;
    break;
  case 8:  // Oxygen
  case 16: // Sulfur
    hcount = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
    hcount = hcount < 2 ? 2 - hcount : 0;
    break;
  case 5:  // Boron
    hcount = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
    hcount = hcount < 3 ? 3 - hcount : 0;
    break;
  case 7:  // Nitogen
  case 15: // Phosphorus
    hcount = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
    if (hcount < 3)
      hcount = 3 - hcount;
    else if (hcount == 4)
      hcount = 1;
    else
      hcount = 0;
    break;
  case 6:  // Carbon
    hcount = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
    hcount = hcount < 4 ? 4 - hcount : 0;
    break;
  default:
    hcount = 0;
  }
  NMS_ATOM_SET_IMPLICITHCOUNT(aptr, hcount);
}

static std::string AnonymousGraph(NMS_pMOL mol, bool elem)
{
  std::string result;
  int charge = 0;

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    charge += NMS_ATOM_GET_FORMALCHARGE(aptr);
    NMS_ATOM_SET_IS_AROMATIC(aptr, false);
    NMS_ATOM_SET_FORMALCHARGE(aptr, 0);
    if (!elem) {
      NMS_ATOM_SET_IMPLICITHCOUNT(aptr, 0);
      NMS_ATOM_SET_ATOMICNUM(aptr, 0);
    } else {
      NormalizeHCount(aptr);
    }
  }

  NMS_FOR_BOND_IN_MOL(bond, mol) {
    NMS_pBOND bptr = NMS_ITER_MOL_BOND(bond, mol);
    NMS_BOND_SET_ORDER(bptr, 1);
  }

  NMS_MOL_ASSIGN_RADICALS(mol);
  NMS_GENERATE_SMILES(mol, result);
  return result;
}

static std::string MesomerHash(NMS_pMOL mol, bool netq)
{
  std::string result;
  char buffer[32];
  int charge = 0;

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    charge += NMS_ATOM_GET_FORMALCHARGE(aptr);
    NMS_ATOM_SET_IS_AROMATIC(aptr, false);
    NMS_ATOM_SET_FORMALCHARGE(aptr, 0);
  }

  NMS_FOR_BOND_IN_MOL(bond, mol) {
    NMS_pBOND bptr = NMS_ITER_MOL_BOND(bond, mol);
    NMS_BOND_SET_ORDER(bptr, 1);
  }

  NMS_MOL_ASSIGN_RADICALS(mol);
  NMS_GENERATE_SMILES(mol, result);
  if (netq) {
    sprintf(buffer,"_%d",charge);
    result += buffer;
  }
  return result;
}

static std::string TautomerHash(NMS_pMOL mol, bool proto)
{
  std::string result;
  char buffer[32];
  int hcount = 0;
  int charge = 0;

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    charge += NMS_ATOM_GET_FORMALCHARGE(aptr);
    NMS_ATOM_SET_IS_AROMATIC(aptr, false);
    NMS_ATOM_SET_FORMALCHARGE(aptr, 0);
    if (NMS_ATOM_GET_ATOMICNUM(aptr) != 6) {
      hcount += NMS_ATOM_GET_IMPLICITHCOUNT(aptr);
      NMS_ATOM_SET_IMPLICITHCOUNT(aptr, 0);
    }
  }

  NMS_FOR_BOND_IN_MOL(bond, mol) {
    NMS_pBOND bptr = NMS_ITER_MOL_BOND(bond, mol);
    NMS_BOND_SET_ORDER(bptr, 1);
  }

  NMS_MOL_ASSIGN_RADICALS(mol);
  NMS_GENERATE_SMILES(mol, result);
  if (!proto) {
    sprintf(buffer, "_%d_%d", hcount, charge);
  } else sprintf(buffer, "_%d", hcount - charge);
  result += buffer;
  return result;
}


static bool TraverseForRing(NMS_pATOM atom,
                            unsigned char *visit)
{
  visit[NMS_ATOM_GET_IDX(atom)] = 1;

  NMS_FOR_NBR_OF_ATOM(nbr, atom) {
    NMS_pATOM nptr = NMS_ITER_ATOM_NBR(nbr, atom);
    if (visit[NMS_ATOM_GET_IDX(nptr)] == 0) {
      if (NMS_ATOM_IS_IN_RING(nptr))
        return true;

      if (TraverseForRing(nptr, visit))
        return true;
    }
  }
  return false;
}


static bool DepthFirstSearchForRing(NMS_pATOM root,
                             NMS_pATOM nbor,
                             unsigned int maxatomidx)
{
  unsigned int natoms = maxatomidx;
  unsigned char *visit = (unsigned char*)alloca(natoms);
  memset(visit, 0, natoms);

  visit[NMS_ATOM_GET_IDX(root)] = true;
  return TraverseForRing(nbor, visit);
}

bool IsInScaffold(NMS_pATOM atom, unsigned int maxatomidx)
{
  if (NMS_ATOM_IS_IN_RING(atom))
    return true;

  unsigned int count = 0;
  NMS_FOR_NBR_OF_ATOM(nbr, atom) {
    NMS_pATOM nptr = NMS_ITER_ATOM_NBR(nbr, atom);
    if (DepthFirstSearchForRing(atom, nptr, maxatomidx))
      ++count;
  }
  return count > 1;
}


static bool HasNbrInScaffold(NMS_pATOM aptr, unsigned char* is_in_scaffold)
{
  NMS_FOR_NBR_OF_ATOM(nbr, aptr) {
    NMS_pATOM nptr = NMS_ITER_ATOM_NBR(nbr, aptr);
    if (is_in_scaffold[NMS_ATOM_GET_IDX(nptr)])
      return true;
  }
  return false;
}

static std::string ExtendedMurckoScaffold(NMS_pMOL mol)
{
  NMS_MOL_CALCULATE_RINGINFO(mol);

  unsigned int maxatomidx = NMS_MOL_GET_MAXATOMIDX(mol);
  unsigned char* is_in_scaffold = (unsigned char*)alloca(maxatomidx);
  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    is_in_scaffold[NMS_ATOM_GET_IDX(aptr)] = IsInScaffold(aptr, maxatomidx);
  }

  std::vector<NMS_pATOM> for_deletion;
  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    unsigned int aidx = NMS_ATOM_GET_IDX(aptr);
    if (is_in_scaffold[aidx]) continue;
    if (HasNbrInScaffold(aptr, is_in_scaffold)) {
      NMS_ATOM_SET_ATOMICNUM(aptr, 0);
      NMS_ATOM_SET_FORMALCHARGE(aptr, 0);
      NMS_ATOM_SET_IMPLICITHCOUNT(aptr, 0);
    } else {
      for_deletion.push_back(aptr);
    }
  }

  for(unsigned int i=0; i<for_deletion.size(); ++i)
    NMS_MOL_DELETE_ATOM(mol, for_deletion[i]);

  NMS_MOL_ASSIGN_RADICALS(mol);
  std::string result;
  NMS_GENERATE_SMILES(mol, result);
  return result;
}

static std::string MurckoScaffoldHash(NMS_pMOL mol)
{
  std::vector<NMS_pATOM> for_deletion;
  do {
    for_deletion.clear();
    NMS_FOR_ATOM_IN_MOL(atom, mol) {
      NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
      unsigned int deg = NMS_ATOM_GET_EXPLICITDEGREE(aptr);
      if (deg < 2) {
        if (deg == 1) { // i.e. not 0 and the last atom in the molecule
          NMS_FOR_BOND_OF_ATOM(bond, aptr) {
            NMS_pBOND bptr = NMS_ITER_ATOM_BOND(bond, aptr);
            NMS_pATOM nbr = NMS_BOND_GET_NBR(bptr, aptr);
            unsigned int hcount = NMS_ATOM_GET_IMPLICITHCOUNT(nbr);
            NMS_ATOM_SET_IMPLICITHCOUNT(nbr, hcount + NMS_BOND_GET_ORDER(bptr));
          }
        }
        for_deletion.push_back(aptr);
      }
    }
    for (unsigned int i = 0; i < for_deletion.size(); ++i) {
      NMS_MOL_DELETE_ATOM(mol, for_deletion[i]);
    }
  } while (!for_deletion.empty());

  NMS_MOL_ASSIGN_RADICALS(mol);
  std::string result;
  NMS_GENERATE_SMILES(mol, result);
  return result;
}

static std::string NetChargeHash(NMS_pMOL mol)
{
  int totalq = 0;

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    totalq += NMS_ATOM_GET_FORMALCHARGE(aptr);
  }

  char buffer[16];
  sprintf(buffer, "%d", totalq);
  return buffer;
}


static std::string SmallWorldHash(NMS_pMOL mol, bool brl)
{
  char buffer[64];

  unsigned int acount = NMS_MOL_GET_NUMATOMS(mol);
  unsigned int bcount = NMS_MOL_GET_NUMBONDS(mol);
  unsigned int rcount = (bcount + 1) - acount;

  if (brl) {
    unsigned int lcount = 0;
    NMS_FOR_ATOM_IN_MOL(atom, mol) {
      NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
      if (NMS_ATOM_GET_EXPLICITDEGREE(aptr) == 2)
        lcount++;
    }
    sprintf(buffer, "B%uR%uL%u", bcount, rcount, lcount);
  }
  else {
    sprintf(buffer, "B%uR%u", bcount, rcount);
  }
  return buffer;
}

static void DegreeVector(NMS_pMOL mol, unsigned int *v)
{
  memset(v, 0, 4 * sizeof(unsigned int));
  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    switch (NMS_ATOM_GET_EXPLICITDEGREE(aptr)) {
    case 4:  v[0]++;  break;
    case 3:  v[1]++;  break;
    case 2:  v[2]++;  break;
    case 1:  v[3]++;  break;
    }
  }
}


static bool HasDoubleBond(NMS_pATOM atom)
{
  NMS_FOR_BOND_OF_ATOM(bond, atom) {
    NMS_pBOND bptr = NMS_ITER_ATOM_BOND(bond, atom);
    if (NMS_BOND_GET_ORDER(bptr) == 2)
      return true;
  }
  return false;
}

// Determine whether/how to fragment bond
// -1 means don't fragment bond
// 0 means break, with hydrogens on both beg and end
// 1 means break, with asterisk on beg and hydrogen on end
// 2 means break, with hydrogen on beg and asterisk on end
// 3 means break, with asterisks on both beg and end

static int RegioisomerBond(NMS_pBOND bnd)
{
  if (NMS_BOND_GET_ORDER(bnd) != 1)
    return -1;
  if (NMS_BOND_IS_IN_RING(bnd))
    return -1;

  NMS_pATOM beg = NMS_BOND_GET_BEG(bnd);
  NMS_pATOM end = NMS_BOND_GET_END(bnd);
  unsigned int beg_elem = NMS_ATOM_GET_ATOMICNUM(beg);
  unsigned int end_elem = NMS_ATOM_GET_ATOMICNUM(end);

  if (beg_elem == 0 || end_elem == 0)
    return -1;

  if (NMS_ATOM_IS_IN_RING(beg)) {
    if (NMS_ATOM_IS_IN_RING(end))
      return 0;
    return 2;
  }
  if (NMS_ATOM_IS_IN_RING(end))
    return 1;

  if (beg_elem != 6 &&
      end_elem == 6 &&
      !HasDoubleBond(end))
    return 1;
  if (beg_elem == 6 &&
      end_elem != 6 &&
      !HasDoubleBond(beg))
    return 2;

  return -1;
}


static void ClearEZStereo(NMS_pATOM atm)
{
  NMS_FOR_BOND_OF_ATOM(bond, atm) {
    NMS_pBOND bptr = NMS_ITER_ATOM_BOND(bond, atm);
    if (NMS_BOND_HAS_STEREO(bptr))
      NMS_BOND_REMOVE_STEREO(bptr);
  }
}


static std::string RegioisomerHash(NMS_pMOL mol)
{
  NMS_MOL_CALCULATE_RINGINFO(mol);

  NMS_FOR_BOND_IN_MOL(bond, mol) {
    NMS_pBOND bptr = NMS_ITER_MOL_BOND(bond, mol);
    int split = RegioisomerBond(bptr);
    if (split >= 0) {
      NMS_pATOM beg = NMS_BOND_GET_BEG(bptr);
      NMS_pATOM end = NMS_BOND_GET_END(bptr);
      NMS_MOL_DELETE_BOND(mol, bptr);

      ClearEZStereo(beg);
      ClearEZStereo(end);

      if (split & 1) {
        NMS_pATOM star = NMS_MOL_NEW_ATOM(mol, 0);
        NMS_MOL_NEW_BOND(mol, beg, star, 1, false);
      } else {
        unsigned int hcount = NMS_ATOM_GET_IMPLICITHCOUNT(beg);
        NMS_ATOM_SET_IMPLICITHCOUNT(beg, hcount + 1);
      }
      if (split & 2) {
        NMS_pATOM star = NMS_MOL_NEW_ATOM(mol, 0);
        NMS_MOL_NEW_BOND(mol, end, star, 1, false);
      } else {
        unsigned int hcount = NMS_ATOM_GET_IMPLICITHCOUNT(end);
        NMS_ATOM_SET_IMPLICITHCOUNT(end, hcount + 1);
      }
    }
  }

  std::string result;
  NMS_GENERATE_SMILES(mol, result);
  return result;
}


static std::string ArthorSubOrderHash(NMS_pMOL mol)
{
  char buffer[256];

  unsigned int acount = NMS_MOL_GET_NUMATOMS(mol);
  unsigned int bcount = NMS_MOL_GET_NUMBONDS(mol);

  unsigned int pcount = 1;
  unsigned int size = 4*NMS_MOL_GET_MAXATOMIDX(mol)+4;
  unsigned int *parts = (unsigned int*)malloc(size);
  if (parts) {
    memset(parts,0,size);
    pcount = NMDetermineComponents(mol, parts, acount);
    free(parts);
  }

  unsigned int ccount = 0;
  unsigned int ocount = 0;
  unsigned int zcount = 0;
  unsigned int icount = 0;
  unsigned int qcount = 0;
  unsigned int rcount = 0;

  NMS_FOR_ATOM_IN_MOL(atom, mol) {
    NMS_pATOM aptr = NMS_ITER_MOL_ATOM(atom, mol);
    unsigned int elem = NMS_ATOM_GET_ATOMICNUM(aptr);
    int charge = NMS_ATOM_GET_FORMALCHARGE(aptr);
    switch (elem) {
    case 6:  // Carbon
      ccount++;
      if (charge==0 && NMS_ATOM_GET_VALENCE(aptr)!=4)
        rcount++;
      break;
    case 7:  // Nitrogen
    case 15: // Phosphorus
      ocount++;
      if (charge==0) {
        unsigned int valence = NMS_ATOM_GET_VALENCE(aptr);
        if (valence!=3 && valence!=5)
          rcount++;
      }
      break;
    case 8:  // Oxygen
      ocount++;
      if (charge && NMS_ATOM_GET_VALENCE(aptr) !=2)
        rcount++;
      break;
    case 9:  // Fluorine
      ocount++;
      if (charge && NMS_ATOM_GET_VALENCE(aptr) !=1)
        rcount++;
      break;
    case 17: // Chlorine
    case 35: // Bromine
    case 53: // Iodine
      ocount++;
      if (charge==0) {
        unsigned int valence = NMS_ATOM_GET_VALENCE(aptr);
        if (valence!=1 && valence!=3 && valence!=5 && valence!=7)
          rcount++;
      }
      break;
    case 16: // Sulfur
      ocount++;
      if (charge==0) {
        unsigned int valence = NMS_ATOM_GET_VALENCE(aptr);
        if (valence!=2 && valence!=4 && valence!=6)
          rcount++;
      }
      break;
    }
    zcount += elem;
    if (NMS_ATOM_GET_ISOTOPE(aptr) != 0)
      icount++;
    if (charge != 0)
      qcount++;
  }

  if (acount > 0xffff)   acount = 0xffff;
  if (bcount > 0xffff)   bcount = 0xffff;
  if (pcount > 0xff)     pcount = 0xff;
  if (ccount > 0xffff)   ccount = 0xffff;
  if (ocount > 0xffff)   ocount = 0xffff;
  if (zcount > 0xffffff) zcount = 0xffffff;
  if (rcount > 0xff)     rcount = 0xff;
  if (qcount > 0xff)     qcount = 0xff;
  if (icount > 0xff)     icount = 0xff;

  sprintf(buffer,"%04x%04x%02x%04x%04x%06x%02x%02x%02x",
          acount,bcount,pcount,ccount,ocount,zcount,rcount,qcount,icount);
  return buffer;
}


std::string MolHash(NMS_pMOL mol, unsigned int func)
{
  std::string result;
  char buffer[32];
  NMS_SANITIZE_HYDROGENS(mol);

  switch (func) {
  default:
  case HashFunction::AnonymousGraph:
    result = AnonymousGraph(mol, false);
    break;
  case HashFunction::ElementGraph:
    result = AnonymousGraph(mol, true);
    break;
  case HashFunction::CanonicalSmiles:
    NMS_GENERATE_SMILES(mol, result);
    break;
  case HashFunction::MurckoScaffold:
    result = MurckoScaffoldHash(mol);
    break;
  case HashFunction::ExtendedMurcko:
    result = ExtendedMurckoScaffold(mol);
    break;
  case HashFunction::Mesomer:
    result = MesomerHash(mol, true);
    break;
  case HashFunction::RedoxPair:
    result = MesomerHash(mol, false);
    break;
  case HashFunction::HetAtomTautomer:
    result = TautomerHash(mol, false);
    break;
  case HashFunction::HetAtomProtomer:
    result = TautomerHash(mol, true);
    break;
  case HashFunction::MolFormula:
    result = NMMolecularFormula(mol);
    break;
  case HashFunction::AtomBondCounts:
    sprintf(buffer, "%u,%u",
            NMS_MOL_GET_NUMATOMS(mol),
            NMS_MOL_GET_NUMBONDS(mol));
    result = buffer;
    break;
  case HashFunction::NetCharge:
    result = NetChargeHash(mol);
    break;
  case HashFunction::SmallWorldIndexBR:
    result = SmallWorldHash(mol, false);
    break;
  case HashFunction::SmallWorldIndexBRL:
    result = SmallWorldHash(mol, true);
    break;
  case HashFunction::DegreeVector:
    {
      unsigned int dv[4];
      DegreeVector(mol, dv);
      sprintf(buffer, "%u,%u,%u,%u", dv[0], dv[1], dv[2], dv[3]);
      result = buffer;
    }
    break;
  case HashFunction::ArthorSubstructureOrder:
    result = ArthorSubOrderHash(mol);
    break;
  case HashFunction::Regioisomer:
    result = RegioisomerHash(mol);
    break;
  }
  return result;
}
