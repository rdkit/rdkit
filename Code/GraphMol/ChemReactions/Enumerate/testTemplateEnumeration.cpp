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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Chemreactions/Enumerate/Template/EnumerateTemplate.h>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

using namespace RDKit;

void testRGroupTemplates() {
  BBS bbs;
  bbs.resize(2);

  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("CC=CCN=C=S")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("C1C=CC1N=C=S")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("CC=CC([O-])N=C=S")) );
  

  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")) );
  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCCc1ncc(Cl)cc1Br")) );
  bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")) );

  // double attachment point
  //bbs[1].push_back( boost::shared_ptr<ROMol>(SmilesToMol("N(C)CCCc1ncc(Cl)cc1Br")) );

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
    "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
    "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]");
  rxn->initReactantMatchers();
  ReactantTemplates templates;
  TEST_ASSERT(ReactantsToTemplates(templates, *rxn, bbs));

  RGROUPS bbIndices(2);
  MOL_SPTR_VECT reactants(2);
  
  for(size_t i=0;i<templates.m_templates[0].size(); ++i) {
    bbIndices[0] = i;
    for(size_t j=0;j<templates.m_templates[1].size(); ++j) {
      bbIndices[1] = j;
      std::string smiles = templates.smiles(bbIndices);
      ROMol *m = SmilesToMol(smiles);
      TEST_ASSERT(m!=0);//, std::string("Bad smiles:") + smiles);
      std::cout << i << " " << j << " " << smiles << std::endl;
      reactants[0] = bbs[0][i];
      reactants[1] = bbs[1][j];
      std::vector<MOL_SPTR_VECT> res = rxn->runReactants( reactants );
      std::string smi1 = MolToSmiles(*m);
      for (size_t idx=0; idx<res.size(); ++idx) {
        for (size_t prod=0; prod<res[idx].size(); ++prod) {
          std::string smi2 = MolToSmiles(*res[idx][prod].get());
          std::cerr << smi1 << " == " << smi2 << std::endl;
          TEST_ASSERT( smi1 == smi2 );
        }
      }
      delete m;
    }
  }
  delete rxn;
}

void testBlankReaction() {
  BBS bbs;
  bbs.resize(1);

  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("c1ccccc1Br")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("c1c(F)c(Cl)ccc1Br")) );
  

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
      "CC>>CC");
  
  rxn->initReactantMatchers();
  ReactantTemplates templates;
  TEST_ASSERT(!ReactantsToTemplates(templates, *rxn, bbs));
}

void testLargeAtomMapFailure() {
  BBS bbs;
  bbs.resize(1);

  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("CC")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("CCC")) );
  

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
      "[C:20]C>>[C:20]C");
  
  rxn->initReactantMatchers();
  ReactantTemplates templates;
  TEST_ASSERT(!ReactantsToTemplates(templates, *rxn, bbs));  
}

void testMondoMappings() {
  BBS bbs;
  bbs.resize(1);

  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("c1ccccc1Br")) );
  bbs[0].push_back( boost::shared_ptr<ROMol>(SmilesToMol("c1c(F)c(Cl)ccc1Br")) );
  

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
      "[c:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1");

  rxn->initReactantMatchers();
  ReactantTemplates templates;
  TEST_ASSERT(ReactantsToTemplates(templates, *rxn, bbs));

  RGROUPS bbIndices(bbs.size());
  MOL_SPTR_VECT reactants(bbs.size());

  std::set<std::string> smiles;
  // Get the Chem Set
  for (size_t i=0;i<bbs[0].size();++i) {
        reactants[0] = bbs[0][i];
        std::vector<MOL_SPTR_VECT> res = rxn->runReactants( reactants );
        
        for (size_t idx=0; idx<res.size(); ++idx) {
          for (size_t prod=0; prod<res[idx].size(); ++prod) {
            smiles.insert(MolToSmiles(*res[idx][prod].get()));
          }
        }
  }

  std::set<std::string> smiles2;

  for(size_t i=0;i<templates.m_templates[0].size(); ++i) {
    bbIndices[0] = i;
    std::string smiles = templates.smiles(bbIndices);
    std::cout << "smiles: " << smiles << std::endl;
    ROMol *m = SmilesToMol(smiles);
    std::string smi1 = MolToSmiles(*m);
    smiles2.insert(smi1);
    delete m;
  }
  TEST_ASSERT(smiles==smiles2);
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

  testRGroupTemplates();
  testBlankReaction();
  testMondoMappings();
  testLargeAtomMapFailure();
}
