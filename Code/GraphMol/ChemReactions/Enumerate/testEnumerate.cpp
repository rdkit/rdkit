//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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

#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/ChemReactions/Enumerate/CartesianProduct.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSample.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSampleAllBBs.h>
#include <GraphMol/Chemreactions/Enumerate/Enumerate.h>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>

using namespace RDKit;

void pickleTest(EnumerationStrategyBase &en, size_t len,
                bool allRGroups) {
  
  boost::shared_ptr<EnumerationStrategyBase> base(en.Clone());
      
  for (size_t i=0; i<len; ++i) {
    std::stringstream ss;
    {
      boost::archive::binary_oarchive ar(ss);
      ar & base;
    }
    boost::shared_ptr<EnumerationStrategyBase> copy;
    {
      boost::archive::binary_iarchive ar(ss);
      ar & copy;
    }
    TEST_ASSERT(base->next() == copy->next());
    TEST_ASSERT(base->currentPosition() == en.next());

  }
}

void testSamplers() {
  BBS bbs;
  bbs.resize(3);
  for(int i=0;i<10;++i)
    bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")) );
  
  for(int i=0;i<5;++i)
    bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")) );

  for(int i=0;i<6;++i)
    bbs[2].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")) );

  ChemicalReaction rxn;
  CartesianProductStrategy    cart; cart.initialize(rxn, bbs);
  RandomSampleStrategy        rand; rand.initialize(rxn, bbs);
  RandomSampleAllBBsStrategy  randBBs; randBBs.initialize(rxn, bbs);
  std::vector< boost::shared_ptr<EnumerationStrategyBase> > enumerators;
  enumerators.push_back(boost::shared_ptr<EnumerationStrategyBase>(cart.Clone()));
  enumerators.push_back(boost::shared_ptr<EnumerationStrategyBase>(rand.Clone()));
  enumerators.push_back(boost::shared_ptr<EnumerationStrategyBase>(randBBs.Clone()));

  for(size_t i=0;i<enumerators.size();++i) {
    TEST_ASSERT(enumerators[i]->getNumPermutations() == 10*5*6);
    pickleTest(*enumerators[i], 10*5*6, i==2);
  }

  //for(auto&& i: enumerators) {
  //  TEST_ASSERT(i->getNumPermutations() == 10*5*6);
  //}
}


const char *smiresults[] = {
  "C=CCNC(=S)NCc1ncc(Cl)cc1Br",
  "CC=CCNC(=S)NCc1ncc(Cl)cc1Br",
  "C=CCNC(=S)NCCc1ncc(Cl)cc1Br",
  "CC=CCNC(=S)NCCc1ncc(Cl)cc1Br",
  "C=CCNC(=S)NCCCc1ncc(Cl)cc1Br",
  "CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br"
};

void testEnumerations() {
  BBS bbs;
  bbs.resize(2);

  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("CC=CCN=C=S")) );

  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")) );
  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCCc1ncc(Cl)cc1Br")) );
  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")) );
  
  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
    "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
    "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]");

  {
    EnumerateLibrary en(*rxn, bbs);
    size_t i=0;
    for(; (bool)en; ++i) {
      std::vector<std::vector<std::string> > res = en.nextSmiles();
      TEST_ASSERT(res.size() == 1);
      TEST_ASSERT(res[0].size() == 1);
      TEST_ASSERT(res[0][0] == smiresults[i]);
    }
    TEST_ASSERT(i==6);
  }

  {
    boost::shared_ptr<EnumerateLibrary> en(new EnumerateLibrary (
        *rxn, bbs, RandomSampleStrategy()));
    size_t i=0;
    for(; i<1000; ++i) {
      // pickle and unpickle
      std::stringstream ss;
      {
        boost::archive::binary_oarchive ar(ss);
        ar & en;
      }
      boost::shared_ptr<EnumerateLibrary> copy;
      {
        boost::archive::binary_iarchive ar(ss);
        ar & copy;
      }

      for(size_t j=0;j<10;++j)
        TEST_ASSERT(en->nextSmiles() == copy->nextSmiles());
    }
    TEST_ASSERT(i==1000);
  }
  delete rxn;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  bool doLong=false;
  if(argc>1) {
    if(!strncmp(argv[1],"-l",2)){
      doLong=true;
    }
  }

  testSamplers();
  testEnumerations();
}
