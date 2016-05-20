//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../RDKitBase.h"
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"

#include "RGroupDecomp.h"

using namespace RDKit;


void test1()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test1\n";
    const char* smols[] = {
      "C4(=CC=C1C(C(=NC=N1)C2=CC(=CC=C2)C(=O)N3CCN(CC3)C(=O)C)=C4)C5C=CC(=NC=5)OC",
    };
    const char* scores[] = {
      "c1ccc2c(cncn2)c1",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    RGroupDecomposite(mols, cores, options, results);
    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

/*
  >>> smis = ['Nc1ccccc1','c1cc(F)ccc1','c1cc(F)c(C)cc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1ccccc1')]
  (None, 'substitution at unmarked position')
*/
void test21()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test21\n";
    const char* smols[] = {
      "Nc1ccccc1",
      "c1cc(F)ccc1",
      "c1cc(F)c(C)cc1",
    };
    const char* scores[] = {
      "[1*]c1ccccc1",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    options.LabelledCores = true;
    options.RequireLabels = true;
    
    RGroupDecomposite(mols, cores, options, results);

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

/*
  >>> smis = ['Nc1cc(F)ccc1','c1cc(F)cc(F)c1','c1c(C)cc(C)cc1','c1c(C)c(C)ccc1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1cc([2*])ccc1')]
  (None, 'substitution at unmarked position')
*/
void test22()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test22\n";
    const char* smols[] = {
      "Nc1cc(F)ccc1",
      "c1cc(F)cc(F)c1",
      "c1c(C)cc(C)cc1",
      "c1c(C)c(C)ccc1",
    };
    const char* scores[] = {
      "[1*]c1cc([2*])ccc1",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    options.LabelledCores = true;
    options.RequireLabels = true;
    
    RGroupDecomposite(mols, cores, options, results);

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

/*
  >>> smis = ['CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4','COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('[1*]c1cccc([2*])c1')]
  >>> mols,d= RunDecomposition(ms,options)
  ((2, '*Nc1nc(nc2c1N=CN2CC)N3CCC(O)CC3'), (1, '*Cl'))
  (None, 'core matches multiple times')
*/
void test23()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test23\n";
    const char* smols[] = {
      "CCN3C=Nc2c(Nc1cccc(Cl)c1)nc(nc23)N4CCC(O)CC4",
      "COc4ccc(CNC2=NNc3ncnc(Nc1cccc(Cl)c1)c23)cc4Cl",
    };
    const char* scores[] = {
      "[1*]c1cccc([2*])c1",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    options.LabelledCores = true;
    options.RequireLabels = true;
    
    RGroupDecomposite(mols, cores, options, results);

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

/*
  >>> smis = ['CC1CC1','C1CCC1']
  >>> ms = [Chem.MolFromSmiles(x) for x in smis]
  >>> options.cores = [Chem.MolFromSmiles('C1CC1[1*]')]
  ((1, '*C'),)
  (None, 'no core matches')
*/
void test24()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "test24\n";
    const char* smols[] = {
      "CC1CC1",
      "C1CCC1",
    };
    const char* scores[] = {
      "C1CC1[1*]",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    options.LabelledCores = true;
    options.RequireLabels = true;
    
    RGroupDecomposite(mols, cores, options, results);

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//--------------------------------------------------------------------------
void testX()
{
    BOOST_LOG(rdInfoLog) << "-------------------------------------\n";
    BOOST_LOG(rdInfoLog) << "testX\n";
    const char* smols[] = {
      "",
      "",
      "",
    };
    const char* scores[] = {
      "",
      "",
      "",
    };
    RGroupDecompositionOptions options;
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> cores;
    std::vector<ROMOL_SPTR> results;

    for(int i=0; i<sizeof(smols)/sizeof(smols[0]); i++) {
        RWMol* m = SmilesToMol(smols[i]);
        TEST_ASSERT(m);
        mols.push_back(ROMOL_SPTR(m));
    }
    for(int i=0; i<sizeof(scores)/sizeof(scores[0]); i++) {
        RWMol* m = SmilesToMol(scores[i]);
        TEST_ASSERT(m);
        cores.push_back(ROMOL_SPTR(m));
    }

    //options. = ;
    
    RGroupDecomposite(mols, cores, options, results);

    BOOST_LOG(rdInfoLog) << "......... results ........\n";
    TEST_ASSERT(true);

    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//==============================================================================

int main(int argc, const char* argv[])
{
    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    BOOST_LOG(rdInfoLog) << "RGroupDecomp Unit Test \n";

    test24();

    test1();
    test21();
    test22();
    test23();
    test24();

    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    return 0;
}


