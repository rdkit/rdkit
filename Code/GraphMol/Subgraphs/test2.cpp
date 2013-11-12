// $Id$
//
//  Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <boost/foreach.hpp>


#include <iostream>
using namespace std;
using namespace RDKit;

void test1()
{
  std::cout << "-----------------------\n Test1: pathToSubmol" << std::endl;
  {
    std::string smiles="CC1CC1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_LIST sgs;
    sgs = findAllSubgraphsOfLengthN(*mol,3,false,0);
    TEST_ASSERT(sgs.size()==3);
    BOOST_FOREACH(PATH_TYPE tmp,sgs){
      TEST_ASSERT(tmp[0]==0);
      TEST_ASSERT(tmp.size()==3);
      ROMol *frag=Subgraphs::pathToSubmol(*mol,tmp,false);
      smiles = MolToSmiles(*frag,true,false,0,false);
      if(tmp[1]==1){
        if(tmp[2]==2){
          TEST_ASSERT(smiles=="CCCC");
        } else if(tmp[2]==3) {
          TEST_ASSERT(smiles=="CC(C)C");
        } else {
          TEST_ASSERT(0);
        }
      } else if(tmp[1]==3){
        if(tmp[2]==2){
          TEST_ASSERT(smiles=="CCCC");
        } else if(tmp[2]==1) {
          TEST_ASSERT(smiles=="CC(C)C");
        } else {
          TEST_ASSERT(0);
        }
      } else {
        TEST_ASSERT(0);
      }
      delete frag;
    } 
    delete mol;
  }
  std::cout << "Finished" << std::endl;
}


void test2()
{
  std::cout << "-----------------------\n Test2: Atom Environments" << std::endl;
  {
    std::string smiles="CC1CC1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mol,1,0);
    TEST_ASSERT(pth.size()==1);
    TEST_ASSERT(pth[0]==0);

    pth=findAtomEnvironmentOfRadiusN(*mol,2,0);
    TEST_ASSERT(pth.size()==3);
    TEST_ASSERT(pth[0]==0);
    
    pth=findAtomEnvironmentOfRadiusN(*mol,3,0);
    TEST_ASSERT(pth.size()==4);
    TEST_ASSERT(pth[0]==0);
    
    pth=findAtomEnvironmentOfRadiusN(*mol,4,0);
    TEST_ASSERT(pth.size()==0);
    
    pth=findAtomEnvironmentOfRadiusN(*mol,1,1);
    TEST_ASSERT(pth.size()==3);

    pth=findAtomEnvironmentOfRadiusN(*mol,2,1);
    TEST_ASSERT(pth.size()==4);

    pth=findAtomEnvironmentOfRadiusN(*mol,3,1);
    TEST_ASSERT(pth.size()==0);

    delete mol;
  }

  {
    std::string smiles="CC1CC1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);
    ROMol *mH=MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mH,1,0);
    TEST_ASSERT(pth.size()==1);
    TEST_ASSERT(pth[0]==0);

    pth=findAtomEnvironmentOfRadiusN(*mH,1,0,true);
    TEST_ASSERT(pth.size()==4);

    delete mol;
    delete mH;
  }

  {
    std::string smiles="O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);
    smiles = MolToSmiles(*mol);

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mol,2,9);
    TEST_ASSERT(pth.size()==8);
    ROMol *frag=Subgraphs::pathToSubmol(*mol,pth,false);
    smiles = MolToSmiles(*frag,true,false,0,false);
    TEST_ASSERT(smiles=="C(C(C(O)C)C(C)C)C");
    delete frag;
    delete mol;
  }


  std::cout << "Finished" << std::endl;
}

void testGithubIssue103()
{
  std::cout << "-----------------------\n Testing github Issue103: stereochemistry and pathToSubmol" << std::endl;
  {
    std::string smiles="O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mol,2,12);
    TEST_ASSERT(pth.size()==5);
    ROMol *frag=Subgraphs::pathToSubmol(*mol,pth,false);
    smiles = MolToSmiles(*frag,true);
    TEST_ASSERT(smiles=="C=CC(O)C=C");
    delete frag;
    delete mol;
  }
  {
    std::string smiles="O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mol,2,12);
    TEST_ASSERT(pth.size()==5);
    ROMol *frag=Subgraphs::pathToSubmol(*mol,pth,false);
    smiles = MolToSmarts(*frag);
    TEST_ASSERT(smiles=="[#6](-[#6H](-[#8])-[#6]=[#6])=[#6]");
    delete frag;
    delete mol;
  }
  {
    std::string smiles="O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol=SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth=findAtomEnvironmentOfRadiusN(*mol,2,12);
    TEST_ASSERT(pth.size()==5);
    ROMol *frag=Subgraphs::pathToSubmol(*mol,pth,true);
    smiles = MolToSmarts(*frag);
    TEST_ASSERT(smiles=="[#6](-[#6](-[#8])-[#6]=[#6])=[#6]");
    delete frag;
    delete mol;
  }

  std::cout << "Finished" << std::endl;
}


// -------------------------------------------------------------------
int main()
{
  test1();
  test2();
  testGithubIssue103();
  return 0;
}
