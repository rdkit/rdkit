/*==============================================*/
/* Copyright (c)  2012-2019  NextMove Software  */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#ifndef NMS_RDKTOOLKIT_H
#define NMS_RDKTOOLKIT_H
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#if NMS_RDKIT_VERSION < 20180301
  #define NMS_RDKIT_ITER_TO_BOND(M,I)  (M)[I].get()
  #define NMS_RDKIT_ITER_TO_ATOM(M,I)  (M)[I].get()
#else
  #define NMS_RDKIT_ITER_TO_BOND(M,I)  (M)[I]
  #define NMS_RDKIT_ITER_TO_ATOM(M,I)  (M)[I]
#endif

#define NMS_MOL			RDKit::RWMol
#define NMS_pMOL		RDKit::RWMol*
#define NMS_MOL_TO_pMOL(X)	(&(X))

#define NMS_pATOM		RDKit::Atom*
#define NMS_pBOND		RDKit::Bond*

#define NMS_MOL_SET_DIMENSION(M,D)	NMRDKitMolSetDimension(M,D)
#define NMS_MOL_SET_TITLE(M,X)		(M)->setProp("_Name",(X))
#define NMS_MOL_SET_TITLEP(M,X)		(M)->setProp("_Name",std::string(X))
#define NMS_MOL_SET_ISRXN(M,X)		(M)->setProp("isRxn",(X)?1:0)
#define NMS_MOL_SET_ATOM_COORDS(M,A,C)   NMRDKitMolSetAtomCoords(M,A,C)
#define NMS_MOL_SET_ATOM_COORDSF(M,A,C)  NMRDKitMolSetAtomCoordsF(M,A,C)

#define NMS_MOL_GET_ATOM_COORDS(M,A,C)	NMRDKitMolGetAtomCoords(M,A,C)
#define NMS_MOL_GET_ATOM_COORDSF(M,A,C)	NMRDKitMolGetAtomCoordsF(M,A,C)
#define NMS_MOL_GET_DIMENSION(M)	NMRDKitMolGetDimension(M)
#define NMS_MOL_GET_ISRXN(M)		NMRDKitMolGetIsRxn(M)
#define NMS_MOL_GET_MAXATOMIDX(M)	(M)->getNumAtoms()
#define NMS_MOL_GET_MAXBONDIDX(M)	(M)->getNumBonds()
#define NMS_MOL_GET_NUMATOMS(M)		(M)->getNumAtoms()
#define NMS_MOL_GET_NUMBONDS(M)		(M)->getNumBonds()
#define NMS_MOL_GET_TITLE(M)		NMRDKitMolGetTitle(M)

#define NMS_MOL_ADD_MOL(D,S)		NMRDKitMolInsert(D,S)
#define NMS_MOL_CLEAR(M)		    NMRDKitMolClear(M)
#define NMS_MOL_COPY_TITLE(D,S)		NMRDKitMolCopyTitle(D,S)
#define NMS_MOL_DELETE_ATOM(M,A)	(M)->removeAtom(A)
#define NMS_MOL_DELETE_BOND(M,B)	\
	(M)->removeBond((B)->getBeginAtomIdx(),(B)->getEndAtomIdx())
#define NMS_MOL_NEW_ATOM(X,Y)		NMRDKitMolNewAtom((X),(Y))
#define NMS_MOL_NEW_BOND(M,B,E,O,A)	NMRDKitMolNewBond((M),(B),(E),(O),(A))
#define NMS_MOL_NONEMPTY(M)		((M)->getNumAtoms()>0)
#define NMS_MOL_ASSIGN_RADICALS(M) RDKit::MolOps::assignRadicals(*(M));
#define NMS_MOL_SUPPRESS_HYDROGENS(M) RDKit::MolOps::removeHs(*(M))
#define NMS_MOL_SPLIT_INTO_FRAGMENTS(M, N)   NMRDKitMolSplitFragments(M, N)
#define NMS_MOL_CALCULATE_RINGINFO(M) NMRDKitMolCalculateRingInfo(M);

#define NMS_ATOM_SET_ATOMICNUM(A,X)		(A)->setAtomicNum(X)
#define NMS_ATOM_SET_FORMALCHARGE(A,X)		(A)->setFormalCharge(X)
#define NMS_ATOM_SET_IMPLICITHCOUNT(A,X)	NMRDKitAtomSetImplicitHCount(A,X)
#define NMS_ATOM_SET_RADICALCOUNT(A,X)    (A)->setNumRadicalElectrons(X)
#define NMS_ATOM_SET_ISOTOPE(A,X)		(A)->setIsotope(X)
#define NMS_ATOM_SET_ALIAS(A,X)    (A)->setProp("molFileAlias",(X))
#define NMS_ATOM_SET_ALIASP(A,X)   (A)->setProp("molFileAlias",\
                                                std::string(X))
#define NMS_ATOM_SET_MAPIDX(A,X)    NMRDKitAtomSetMapIdx(A,X)
#define NMS_ATOM_SET_RXNROLE(A,X)   (A)->setProp("molRxnRole",(int)(X))
#define NMS_ATOM_SET_RXNGROUP(A,X)  (A)->setProp("molRxnComponent",(int)(X))
#define NMS_ATOM_SET_IS_AROMATIC(A,X)	(A)->setIsAromatic(X)
#define NMS_ATOM_IS_IN_RING(X)      RDKit::queryIsAtomInRing(X)

#define NMS_ATOM_GET_ATOMICNUM(A)       (A)->getAtomicNum()
#define NMS_ATOM_GET_BOND(A,B)  \
	(A)->getOwningMol().getBondBetweenAtoms((A)->getIdx(),(B)->getIdx())
#define NMS_ATOM_GET_EXPLICITDEGREE(A)  (A)->getDegree()
#define NMS_ATOM_GET_EXPLICITVALENCE(A)  NMRDKitAtomGetExplicitValence(A)
#define NMS_ATOM_GET_VALENCE(A)         (A)->getTotalValence()
#define NMS_ATOM_GET_FORMALCHARGE(A)    (A)->getFormalCharge()
#define NMS_ATOM_GET_IDX(A)             (A)->getIdx()
#define NMS_ATOM_GET_IMPLICITHCOUNT(A)  (A)->getTotalNumHs(false)
#define NMS_ATOM_GET_RADICALCOUNT(A)    (A)->getNumRadicalElectrons()
#define NMS_ATOM_GET_ISOTOPE(A)         (A)->getIsotope()
#define NMS_ATOM_GET_MAPIDX(A)          NMRDKitAtomGetMapIdx(A)
#define NMS_ATOM_GET_TOTALHCOUNT(A)     (A)->getTotalNumHs(true)
#define NMS_ATOM_IS_CONNECTED(A, B)     NMRDKitAtomIsConnected(A, B)
#define NMS_ATOM_IS_AROMATIC(A)         (A)->getIsAromatic()
#define NMS_ATOM_HAS_STEREO(A)         ((A)->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED)
#define NMS_ATOM_REMOVE_STEREO(A)       (A)->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED)

#define NMS_BOND_SET_ORDER(B,O)  NMRDKitBondSetOrder(B,O)
#define NMS_BOND_SET_WEDGE(B)    (B)->setBondDir(RDKit::Bond::BEGINWEDGE)
#define NMS_BOND_SET_HASH(B)     (B)->setBondDir(RDKit::Bond::BEGINDASH)

#define NMS_BOND_GET_BEG(B)      (B)->getBeginAtom()
#define NMS_BOND_GET_BEGIDX(B)   (B)->getBeginAtomIdx()
#define NMS_BOND_GET_END(B)      (B)->getEndAtom()
#define NMS_BOND_GET_ENDIDX(B)   (B)->getEndAtomIdx()
#define NMS_BOND_GET_IDX(B)      (B)->getIdx()
#define NMS_BOND_GET_NBR(B,A)    (B)->getOtherAtom(A)
#define NMS_BOND_GET_ORDER(B)    NMRDKitBondGetOrder(B)
#define NMS_BOND_IS_AROMATIC(B)  \
	((B)->getBondType() == RDKit::Bond::AROMATIC)
#define NMS_BOND_IS_DOUBLE(B) \
	((B)->getBondType() == RDKit::Bond::DOUBLE)
#define NMS_BOND_IS_IN_RING(X)      RDKit::queryIsBondInRing(X)
#define NMS_BOND_HAS_STEREO(B)   ((B)->getStereo() > RDKit::Bond::STEREOANY)
#define NMS_BOND_REMOVE_STEREO(B) (B)->setStereo(RDKit::Bond::STEREOANY)


#define NMS_FOR_ATOM_IN_MOL(X,Y) \
  for (RDKit::RWMol::AtomIterator X=(Y)->beginAtoms();X!=(Y)->endAtoms();++X)
#define NMS_ITER_MOL_ATOM(X,Y) (*(X))
#define NMS_FOR_BOND_IN_MOL(X,Y) \
  for (RDKit::RWMol::BondIterator X=(Y)->beginBonds();X!=(Y)->endBonds();++X)
#define NMS_ITER_MOL_BOND(X,Y) (*(X))

#define NMS_FOR_BOND_OF_ATOM(X,Y) \
  RDKit::ROMol::OEDGE_ITER X##_beg, X##_end; \
  for (boost::tie(X##_beg,X##_end) = Y->getOwningMol().getAtomBonds(Y); \
       X##_beg != X##_end; ++X##_beg)
#define NMS_ITER_ATOM_BOND(X,Y) \
	NMS_RDKIT_ITER_TO_BOND(Y->getOwningMol(),*X##_beg)

#define NMS_FOR_NBR_OF_ATOM(X,Y) \
  RDKit::ROMol::ADJ_ITER X##_beg, X##_end; \
  for (boost::tie(X##_beg,X##_end) = Y->getOwningMol().getAtomNeighbors(Y); \
       X##_beg != X##_end; ++X##_beg)
#define NMS_ITER_ATOM_NBR(X,Y) \
	NMS_RDKIT_ITER_TO_ATOM(Y->getOwningMol(),*X##_beg)


#define NMS_GENERATE_SMILES(M,S)	S = RDKit::MolToSmiles(*(M));
#define NMS_SMILES_TO_MOL(X)      NMRDKitSmilesToMol(X)
#define NMS_SANITIZE_HYDROGENS(M) NMRDKitSanitizeHydrogens(M)


/* Implemented in rdktoolkit.cpp */
RDKit::RWMol *NMRDKitSmilesToMol(const char *str);
RDKit::Atom *NMRDKitMolNewAtom(RDKit::RWMol *mol, unsigned int elem);
RDKit::Bond *NMRDKitMolNewBond(RDKit::RWMol *mol,
                               RDKit::Atom *src, RDKit::Atom *dst,
                               unsigned int order, bool arom);
unsigned int NMRDKitAtomGetExplicitValence(RDKit::Atom *atm);
unsigned int NMRDKitBondGetOrder(const RDKit::Bond *bnd);
void NMRDKitAtomSetImplicitHCount(RDKit::Atom *atm, unsigned int hcount);
void NMRDKitBondSetOrder(RDKit::Bond *bnd, unsigned int order);
void NMRDKitSanitizeHydrogens(RDKit::RWMol *mol);
void NMRDKitAtomSetMapIdx(RDKit::Atom *atm, unsigned int idx);
void NMRDKitMolSplitFragments(NMS_pMOL mol, std::vector<NMS_MOL*> &fragments);
void NMRDKitMolCalculateRingInfo(NMS_pMOL mol);
std::string NMRDKitMolGetTitle(RDKit::RWMol *mol);

#endif // NMS_RDKTOOLKIT_H

