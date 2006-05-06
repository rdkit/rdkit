// $Id$
//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>


#include <iostream>
using namespace std;
using namespace RDKit;

void test1(){
  Atom *a1 = new Atom(6);
  Atom *a2 = new Atom(*a1);
  Atom *a3 = a1->copy();
  delete a1;
  delete a2;
  delete a3;

}

void test2(){
  RWMol *mol = new RWMol();

  mol->addAtom(new Atom(6));
  mol->addAtom(new Atom(6));
  mol->addAtom(new Atom(8));
  mol->addBond(0,1,Bond::SINGLE);
  mol->addBond(1,2,Bond::SINGLE);
  mol->setAtomBookmark(mol->getAtomWithIdx(1),1);
  mol->setBondBookmark(mol->getBondWithIdx(0),2);
  CHECK_INVARIANT(mol->hasAtomBookmark(1),"");
  CHECK_INVARIANT(mol->getAtomWithBookmark(1)->getIdx()==1,"");
  CHECK_INVARIANT(mol->hasBondBookmark(2),"");
  CHECK_INVARIANT(mol->getBondWithBookmark(2)->getIdx()==0,"");

  RWMol *mol2 = new RWMol(*mol);
  CHECK_INVARIANT(mol2->hasAtomBookmark(1),"");
  CHECK_INVARIANT(mol2->getAtomWithBookmark(1)->getIdx()==1,"");
  CHECK_INVARIANT(mol2->hasBondBookmark(2),"");
  CHECK_INVARIANT(mol2->getBondWithBookmark(2)->getIdx()==0,"");

}

// -------------------------------------------------------------------
int main()
{
  test1();
  test2();
  return 0;
}
