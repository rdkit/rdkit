// $Id$
//
//  Created by Greg Landrum, July 2008
//

//
//  Expected test results here correspond to v1.0 of the open-source avalontoolkit
//


#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h> 
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h> 
#include <DataStructs/ExplicitBitVect.h>

#include "AvalonTools.h"

#include <string>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) << "testing canonical smiles generation" << std::endl;

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccc1"));
    TEST_ASSERT(m);
    std::string smi=AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi=="c1ccccc1");
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1cccnc1"));
    TEST_ASSERT(m);
    std::string smi=AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi=="c1ccncc1");
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("n1ccccc1"));
    TEST_ASSERT(m);
    std::string smi=AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi=="c1ccncc1");
    delete m;
  }
  {
    std::string smi=AvalonTools::getCanonSmiles("n1ccccc1",true);
    TEST_ASSERT(smi=="c1ccncc1");
  }
  {
    std::string smi=AvalonTools::getCanonSmiles("c1cccnc1",true);
    TEST_ASSERT(smi=="c1ccncc1");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void test2(){
  BOOST_LOG(rdInfoLog) << "testing coordinate generation" << std::endl;

#if 0
  {
    RWMol *m = SmilesToMol("c1cccnc1");
    TEST_ASSERT(m);

    unsigned int confId=AvalonTools::set2DCoords(*m);
    TEST_ASSERT(m->getNumConformers()==1);
    TEST_ASSERT(confId==0);
    delete m;
  }
#endif
  {
    std::string molb = AvalonTools::set2DCoords("c1cccnc1",true);
    TEST_ASSERT(molb!="");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test3(){
  BOOST_LOG(rdInfoLog) << "testing fingerprint generation" << std::endl;

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccn1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m,bv,512,false,true,0x00006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccn1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==18);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m,bv,512,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccn1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==6);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1nnccc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m,bv,512,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1nnccc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==28);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ncncc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m,bv,512,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ncncc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==25);
    delete m;
  }
  {
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP("c1cccnc1",true,bv,512,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cccnc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==18);
  }
  {
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP("c1ccccc1",true,bv,512,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==6);
  }

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1cccnc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(1024);
    AvalonTools::getAvalonFP(*m,bv,1024,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cccnc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==19);
    delete m;
  }
  {
    ExplicitBitVect bv(2048);
    AvalonTools::getAvalonFP("c1cocc1",true,bv,2048,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cocc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==48);
  }
  {
    ExplicitBitVect bv(2048);
    AvalonTools::getAvalonFP("C1=COC=C1",true,bv,2048,false,true,0x006FFF);
    BOOST_LOG(rdInfoLog) << "C1=COC=C1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits()==48);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}


void testRDK151(){
  BOOST_LOG(rdInfoLog) << "testing Jira issue RDK-151:  pyAvalonTools not generating chiral smiles from molecules" << std::endl;

  {
    std::string tSmi="C[C@H](F)Cl";
    ROMol *m = static_cast<ROMol *>(SmilesToMol(tSmi));
    TEST_ASSERT(m);
    std::string smi=AvalonTools::getCanonSmiles(tSmi,true);
    CHECK_INVARIANT(smi==tSmi,smi+"!="+tSmi);
    smi=AvalonTools::getCanonSmiles(*m);
    CHECK_INVARIANT(smi==tSmi,smi+"!="+tSmi);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSmilesFailures(){
  BOOST_LOG(rdInfoLog) << "testing handling of bad smiles strings" << std::endl;

  {
    std::string tSmi="C1C";
    std::string smi=AvalonTools::getCanonSmiles(tSmi,true);
    CHECK_INVARIANT(smi=="",smi);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSubstructFps(){
  BOOST_LOG(rdInfoLog) << "testing substructure fingerprints " << std::endl;
  {
    ExplicitBitVect bv1(512),bv2(512);
    AvalonTools::getAvalonFP("c1ccccc1",true,bv1,512,true,true,AvalonTools::avalonSSSBits);
    AvalonTools::getAvalonFP("c1ccccc1C(F)(F)F",true,bv2,512);
    TEST_ASSERT((bv1&bv2)==bv1);
    AvalonTools::getAvalonFP("c1ccccc1C(F)(F)F",true,bv1,512);
    TEST_ASSERT((bv1&bv2)==bv1);
    AvalonTools::getAvalonFP("c1cccc(C)c1C(F)(F)F",true,bv2,512);
    TEST_ASSERT((bv1&bv2)==bv1);
  }
  {
    ExplicitBitVect bv1(512),bv2(512);
    AvalonTools::getAvalonFP("c1ccccc1O",true,bv1,512,true,true,AvalonTools::avalonSSSBits);
    AvalonTools::getAvalonFP("c1ccccc1OC",true,bv2,512);
    TEST_ASSERT((bv1&bv2)==bv1);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testStruChk(){
  BOOST_LOG(rdInfoLog) << "testing structure checking " << std::endl;
  {
    int errs = 0;
    RDKit::ROMOL_SPTR m = AvalonTools::checkMol(errs, "c1ccccc1",true);
    TEST_ASSERT(errs==0);
    m = AvalonTools::checkMol(errs, "c1c(R)cccc1C1(CC-C(C)C1)C",true);
    TEST_ASSERT(errs!=0);
  }
  {
    int errs = 0;
    std::string res;
    boost::tie(res,errs)=AvalonTools::checkMolString("c1ccccc1",true);
    TEST_ASSERT(errs==0);
    TEST_ASSERT(res!="");
    boost::tie(res,errs)=AvalonTools::checkMolString("c1c(R)cccc1C1(CC-C(C)C1)C",true);
    TEST_ASSERT(errs==1);
    TEST_ASSERT(res=="");
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBadMolfile() {
  BOOST_LOG(rdInfoLog) << "testing handling bad molecules " << std::endl;
  // some tests around dealing with bad mol blocks
  {
    std::string molb="SNAP007157A\n\
  MACCS-II3194121345\n\
\n\
0    0  0  0  0";
    std::string smi=AvalonTools::getCanonSmiles(molb,false);
    CHECK_INVARIANT(smi=="",smi);

    ExplicitBitVect bv(1024);
    AvalonTools::getAvalonFP(molb,false,bv,1024);
    TEST_ASSERT(bv.getNumOnBits()==0);
    
    std::string oMolb;
    AvalonTools::set2DCoords(molb,false);
    CHECK_INVARIANT(oMolb=="",oMolb);
    
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc,char *argv[]){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  testRDK151();
  testSmilesFailures();
  testSubstructFps();
  testStruChk();
  testBadMolfile();

  return 0;
}

