//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior
//       written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  created by Nik Stiefl May 2008
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h> 
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

#include <string>

using namespace RDKit;

void testGeneral(std::string rdbase){

  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- testing general mol2 file parsing --" << std::endl;
  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;

  {
    bool ok=false;
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/nonExistFile.mol2";
    try{
      RWMol *m = Mol2FileToMol(fName);
      delete m;
    } catch(const BadFileException &e){
      ok=true;
    }
    TEST_ASSERT(ok);
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/pyrazole_pyridine.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==5);

    // this was sf.net issue 2727976:
    TEST_ASSERT(m->getNumConformers()==1);
    TEST_ASSERT(m->getConformer().is3D());
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).x,1.5019));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).y,1.0435));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).z,0.0000));

    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/benzene.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==6);
    delete m;
  }
  {
    bool ok=false;
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/mol_noatoms.mol2";
    try {
      RWMol *m = Mol2FileToMol(fName);
      delete m;
    } catch(const FileParseException &e){
      ok=true;
    }
    TEST_ASSERT(ok);
  }
  {
    bool ok=false;
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/mol_nomol.mol2";
    try {
      RWMol *m = Mol2FileToMol(fName);
      delete m;
    } catch(const FileParseException &e){
      ok=true;
    }
    TEST_ASSERT(ok);
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/lonePairMol.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==5 && m->getNumBonds()==4);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/symmetricGuanidine.mol2";
    RWMol *m = Mol2FileToMol(fName,false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==1);
    delete m;
  }

  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/highlySymmetricGuanidine.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(4)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(12)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(20)->getFormalCharge()==1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Noxide.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(9)->getFormalCharge()==-1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Noxide.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(9)->getFormalCharge()==-1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/fusedRing.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);


    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(5)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==0);
    TEST_ASSERT(m->getAtomWithIdx(13)->getFormalCharge()==0);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/pyridiniumPhenyl.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(5)->getFormalCharge()==1);
    TEST_ASSERT(m->getAtomWithIdx(6)->getFormalCharge()==0);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/sulfonAmide.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge()==0);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/chargedAmidineRWH.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(6)->getFormalCharge()==1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/chargedAmidineEC.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(3)->getFormalCharge()==1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/chargedAmidine.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(9)->getFormalCharge()==1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/dbtranslateCharged.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==1);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/dbtranslateUncharged.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge()==0);
    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/dbtranslateUnchargedRing.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(2)->getFormalCharge()==0);
    delete m;
  }


#if 0  
  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Sulfonate.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    BOOST_LOG(rdInfoLog) <<MolToSmiles(*m)<<std::endl;
    delete m;
  }
#endif

  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- DONE general mol2 file parsing --" << std::endl;
  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;

}

void testAromaticChargedFail(std::string rdbase){

  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- testing subst arom groups in mol2 --" << std::endl;
  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;

  //this one checks on substituted aromatic groups that need to be charged, e.g. c1ccccn1C which
  //should convert to c1cccc[n+]1C

  {
    //this one is supposed to have a sanitisation error!
    //have to fix that one next ...
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/badSubstPyridine.mol2";
    RWMol *m = Mol2FileToMol(fName);
    if(m){
      delete m;
    }
  }

  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- DONE subst arom groups in mol2 --" << std::endl;
  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
}

void testIssue3399798(std::string rdbase){

  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- testing issue 3399798             --" << std::endl;
  BOOST_LOG(rdInfoLog) << "---------------------------------------" << std::endl;

  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3399798.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()==Atom::CHI_UNSPECIFIED);

    delete m;
  }

  {
    std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3399798.2.mol2";
    RWMol *m = Mol2FileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag()==Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag()!=Atom::CHI_UNSPECIFIED);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- DONE                           --" << std::endl;
  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
}


// Test the following GitHub issue: https://github.com/rdkit/rdkit/issues/114
void testIssue114(std::string rdbase) {

  BOOST_LOG(rdInfoLog) << "-----------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- testing GitHub issue #114     --" << std::endl;
  BOOST_LOG(rdInfoLog) << "-----------------------------------" << std::endl;

  std::string fName = rdbase + "/Code/GraphMol/FileParsers/test_data/EZ_mol2_issue114.mol2";
  RWMol *mol = Mol2FileToMol(fName);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  delete mol;
  
  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "-- DONE                           --" << std::endl;
  BOOST_LOG(rdInfoLog) << "------------------------------------" << std::endl;
}

  //FIX still missing chirality by 3D structure  
  //  still missing input std::string

int main(int argc,char *argv[]){
  RDLog::InitLogs();

  std::string rdbase = getenv("RDBASE");

  testGeneral(rdbase);
  testAromaticChargedFail(rdbase);
  testIssue3399798(rdbase);
  testIssue114(rdbase);
  
  return 0;
}
