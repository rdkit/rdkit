// $Id$
//
//  Copyright (C) 2008 Greg Landrum
//   All Rights Reserved
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>

using namespace RDKit;

void BuildSimpleMolecule(){
  // build the molecule: C/C=C\C
  RWMol *mol=new RWMol();

  // add atoms and bonds:
  mol->addAtom(new Atom(6)); // atom 0
  mol->addAtom(new Atom(6)); // atom 1
  mol->addAtom(new Atom(6)); // atom 2
  mol->addAtom(new Atom(6)); // atom 3
  mol->addBond(0,1,Bond::SINGLE); // bond 0
  mol->addBond(1,2,Bond::DOUBLE); // bond 1
  mol->addBond(2,3,Bond::SINGLE); // bond 2
  // setup the stereochem:
  mol->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
  mol->getBondWithIdx(2)->setBondDir(Bond::ENDDOWNRIGHT);
  
  // do the chemistry perception:
  RDKit::MolOps::sanitizeMol(*mol);

  // Get the canonical SMILES, include stereochemistry:
  std::string smiles;
  smiles = MolToSmiles(*(static_cast<ROMol *>(mol)),true); 
  BOOST_LOG(rdInfoLog)<<" sample 1 SMILES: " <<smiles<<std::endl;
}

void WorkWithRingInfo(){
  // use a more complicated molecule to demonstrate querying about
  // ring information
  ROMol *mol=SmilesToMol("OC1CCC2C1CCCC2");
  // the molecule from SmilesToMol is already sanitized, so we don't
  // need to worry about that.

  // work with ring information
  RingInfo *ringInfo = mol->getRingInfo();

  TEST_ASSERT(ringInfo->numRings()==2);

  // can ask how many rings an atom is in:
  TEST_ASSERT(ringInfo->numAtomRings(0)==0);
  TEST_ASSERT(ringInfo->numAtomRings(1)==1);
  TEST_ASSERT(ringInfo->numAtomRings(4)==2);
  // same with bonds:
  TEST_ASSERT(ringInfo->numBondRings(0)==0);
  TEST_ASSERT(ringInfo->numBondRings(1)==1);

  // can check if an atom is in a ring of a particular size:
  TEST_ASSERT(!ringInfo->isAtomInRingOfSize(0,5));
  TEST_ASSERT(ringInfo->isAtomInRingOfSize(1,5));
  TEST_ASSERT(ringInfo->isAtomInRingOfSize(4,5));
  TEST_ASSERT(ringInfo->isAtomInRingOfSize(4,6));
  // same with bonds:
  TEST_ASSERT(!ringInfo->isBondInRingOfSize(0,5));
  TEST_ASSERT(ringInfo->isBondInRingOfSize(1,5));

  // can also get the full list of rings as atom indices:
  VECT_INT_VECT atomRings; // VECT_INT_VECT is vector< vector<int> >
  atomRings=ringInfo->atomRings();
  TEST_ASSERT(atomRings.size()==2);
  TEST_ASSERT(atomRings[0].size()==5);
  TEST_ASSERT(atomRings[1].size()==6);
  // this sort is just here for test/demo purposes:
  std::sort(atomRings[0].begin(),atomRings[0].end());
  TEST_ASSERT(atomRings[0][0]==1);
  TEST_ASSERT(atomRings[0][1]==2);
  TEST_ASSERT(atomRings[0][2]==3);
  TEST_ASSERT(atomRings[0][3]==4);
  TEST_ASSERT(atomRings[0][4]==5);
  // same with bonds:
  VECT_INT_VECT bondRings; // VECT_INT_VECT is vector< vector<int> >
  bondRings=ringInfo->bondRings();
  TEST_ASSERT(bondRings.size()==2);
  TEST_ASSERT(bondRings[0].size()==5);
  TEST_ASSERT(bondRings[1].size()==6);
  // the same trick played above with the contents of each ring
  // can be played, but we won't
}

void WorkWithSmarts(){
  // demonstrate the use of substructure searching
  ROMol *mol=SmilesToMol("ClCC=CCC");
  // a simple SMARTS pattern for rotatable bonds:
  ROMol *pattern=SmartsToMol("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");

  std::vector<MatchVectType> matches;
  unsigned int nMatches;
  nMatches=SubstructMatch(*mol,*pattern,matches);
  TEST_ASSERT(nMatches==2);
  TEST_ASSERT(matches.size()==2); // <- there are two rotatable bonds

  // a MatchVect is a vector of std::pairs with (patternIdx, molIdx):
  TEST_ASSERT(matches[0].size()==2);
  TEST_ASSERT(matches[0][0].first==0);
  TEST_ASSERT(matches[0][0].second==1);
  TEST_ASSERT(matches[0][1].first==1);
  TEST_ASSERT(matches[0][1].second==2);
  
}

void DepictDemo(){
  // demonstrate the use of the depiction-generation code2D coordinates:
  ROMol *mol=SmilesToMol("ClCC=CCC");

  // generate the 2D coordinates:
  RDDepict::compute2DCoords(*mol);

  // generate a mol block (could also go to a file):
  std::string molBlock=MolToMolBlock(*mol);
  BOOST_LOG(rdInfoLog)<<molBlock;
}

int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  BuildSimpleMolecule();
  WorkWithRingInfo();
  WorkWithSmarts();
  DepictDemo();
}
