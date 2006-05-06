// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Subgraphs/Subgraphs.h>


#include <iostream>
using namespace std;
using namespace RDKit;

void dumpVIV(PATH_LIST v){
  PATH_LIST::iterator i;
  PATH_TYPE::iterator j;
  for(i=v.begin();i!=v.end();i++){
    for(j=i->begin();j!=i->end();j++){
      std::cout << *j << " ";
    }
    std::cout << std::endl;
  }

}


void test1()
{
  std::cout << "-----------------------\n Test1" << std::endl;
  RWMol *mol=SmilesToMol("CCC(O)C(c1ccccc1)CC(C)N(C)C");
  CHECK_INVARIANT(mol,"");

  int nAll=0,nUnique=0;
  for(int i=1;i<18;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    nUnique += tmp.size();
  }
  CHECK_INVARIANT(nAll==1990,"");
  CHECK_INVARIANT(nUnique==907,"");

  std::cout << "Finished" << std::endl;
}


void test2()
{
  std::cout << "-----------------------\n Test2" << std::endl;
  RWMol *mol=SmilesToMol("C12C3C4C1C1C2C3N41");
  CHECK_INVARIANT(mol,"");

  int nAll=0,nUnique=0;
  for(int i=1;i<13;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    nUnique += tmp.size();
  }

  CHECK_INVARIANT(nAll==2433,"");
  CHECK_INVARIANT(nUnique==300,"");

  std::cout << "Finished" << std::endl;
}

void test3()
{
  std::cout << "-----------------------\n Test3" << std::endl;
  RWMol *mol=SmilesToMol("O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))");
  CHECK_INVARIANT(mol,"");

  int nAll=0,nUnique=0;
  for(int i=1;i<26;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    //std::cout << i << "\t" << tmp.size();
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    //std::cout << "\t" << tmp.size() << std::endl;;
    nUnique += tmp.size();
  }
  std::cout << nAll << " " << nUnique << std::endl;
  CHECK_INVARIANT(nAll==6435,"");
  CHECK_INVARIANT(nUnique==5618,"");

  std::cout << "Finished" << std::endl;
}


// -------------------------------------------------------------------
int main()
{
  test1();
  test2();
  test3();
  return 0;
}
