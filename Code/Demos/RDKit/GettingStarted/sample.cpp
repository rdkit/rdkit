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

  // count the number of rings of size 5:
  unsigned int nRingsSize5=0;
  for(VECT_INT_VECT_CI ringIt=atomRings.begin();
      ringIt!=atomRings.end();++ringIt){
    if(ringIt->size()==5) nRingsSize5++;
  }
  TEST_ASSERT(nRingsSize5==1);
  delete mol;

  // count the number of atoms in 5-rings where all the atoms
  // are aromatic:
  mol=SmilesToMol("C1CC2=C(C1)C1=C(NC3=C1C=CC=C3)C=C2");
  ringInfo = mol->getRingInfo();
  atomRings=ringInfo->atomRings();

  unsigned int nMatchingAtoms=0;
  for(VECT_INT_VECT_CI ringIt=atomRings.begin();
      ringIt!=atomRings.end();++ringIt){
    if(ringIt->size()!=5){
      continue;
    }
    bool isAromatic=true;
    for(INT_VECT_CI atomIt=ringIt->begin();
        atomIt!=ringIt->end();++atomIt){
      if(!mol->getAtomWithIdx(*atomIt)->getIsAromatic()){
        isAromatic=false;
        break;
      }
    }
    if(isAromatic){
      nMatchingAtoms+=5;
    }
  }
  TEST_ASSERT(nMatchingAtoms==5);
  delete mol;

  // count the number of rings where all the bonds
  // are aromatic.
  mol=SmilesToMol("c1cccc2c1CCCC2");
  ringInfo = mol->getRingInfo();
  bondRings=ringInfo->bondRings();
  
  unsigned int nAromaticRings=0;
  for(VECT_INT_VECT_CI ringIt=bondRings.begin();
      ringIt!=bondRings.end();++ringIt){
    bool isAromatic=true;
    for(INT_VECT_CI bondIt=ringIt->begin();
        bondIt!=ringIt->end();++bondIt){
      if(!mol->getBondWithIdx(*bondIt)->getIsAromatic()){
        isAromatic=false;
        break;
      }
    }
    if(isAromatic) nAromaticRings++;
  }
  TEST_ASSERT(nAromaticRings==1);
  delete mol;
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

  delete pattern;
  delete mol;
}

void DepictDemo(){
  // demonstrate the use of the depiction-generation code2D coordinates:
  ROMol *mol=SmilesToMol("ClCC=CCC");

  // generate the 2D coordinates:
  RDDepict::compute2DCoords(*mol);

  // generate a mol block (could also go to a file):
  std::string molBlock=MolToMolBlock(*mol);
  BOOST_LOG(rdInfoLog)<<molBlock;

  delete mol;
}


void CleanupMolecule(){
  // build: C1CC1C(:O):O
  RWMol *mol=new RWMol();

  // add atoms and bonds:
  mol->addAtom(new Atom(6)); // atom 0
  mol->addAtom(new Atom(6)); // atom 1
  mol->addAtom(new Atom(6)); // atom 2
  mol->addAtom(new Atom(6)); // atom 3
  mol->addAtom(new Atom(8)); // atom 4
  mol->addAtom(new Atom(8)); // atom 5
  mol->addBond(3,4,Bond::AROMATIC); // bond 0
  mol->addBond(3,5,Bond::AROMATIC); // bond 1
  mol->addBond(3,2,Bond::SINGLE); // bond 2
  mol->addBond(2,1,Bond::SINGLE); // bond 3
  mol->addBond(1,0,Bond::SINGLE); // bond 4
  mol->addBond(0,2,Bond::SINGLE); // bond 5
  
  // instead of calling sanitize mol, which would generate an error,
  // we'll perceive the rings, then take care of aromatic bonds
  // that aren't in a ring, then sanitize:
  MolOps::findSSSR(*mol);
  for(ROMol::BondIterator bondIt=mol->beginBonds();
      bondIt!=mol->endBonds();++bondIt){
    if( ((*bondIt)->getIsAromatic() ||
         (*bondIt)->getBondType()==Bond::AROMATIC)
        && !mol->getRingInfo()->numBondRings((*bondIt)->getIdx()) ){
      (*bondIt)->setIsAromatic(false);
      // NOTE: this isn't really reasonable:
      (*bondIt)->setBondType(Bond::SINGLE);      
    }
  }
    
  // now it's safe to sanitize:
  RDKit::MolOps::sanitizeMol(*mol);

  // Get the canonical SMILES, include stereochemistry:
  std::string smiles;
  smiles = MolToSmiles(*(static_cast<ROMol *>(mol)),true); 
  BOOST_LOG(rdInfoLog)<<" fixed SMILES: " <<smiles<<std::endl;
}



int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  CleanupMolecule();
  BuildSimpleMolecule();
  WorkWithRingInfo();
  WorkWithSmarts();
  DepictDemo();
}
