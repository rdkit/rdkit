// $Id$
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
//  Created by Greg Landrum, Sept. 2006
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/SLNParse/InputFiller.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
int yysln_parse (void);
extern int yysln_debug; 

bool slnParserDoQueries;

namespace RDKit {
  namespace SLNParse {
    std::vector<RDKit::RWMol *> molList_g;


    RWMol *finalizeQueryMol(ROMol *mol,bool mergeHs){
      PRECONDITION(mol,"bad query molecule");

      // do we need to remove the Hs from the molecule?
      if(mergeHs){
        mol = MolOps::removeHs(*mol,false,true,false);
        for(ROMol::AtomIterator atomIt=mol->beginAtoms();
            atomIt!=mol->endAtoms();++atomIt){
          // if we did remove Hs from this atom (now present as explicit Hs), 
          // set a query for the H count:
          if((*atomIt)->getNumExplicitHs()){
            (*atomIt)->expandQuery(makeAtomHCountQuery((*atomIt)->getNumExplicitHs()));
          }
        }
        
      }

      // we don't want to sanitize, but we do need to get 
      // some ring info:
      VECT_INT_VECT sssr;
      MolOps::symmetrizeSSSR(*mol,sssr);
      for(ROMol::AtomIterator atomIt=mol->beginAtoms();
          atomIt!=mol->endAtoms();++atomIt){
        SLNParse::parseFinalAtomAttribs(*atomIt,true);
      }
      
      return static_cast<RWMol *>(mol);
    }

    
    RWMol *toMol(std::string inp,int func(void),bool doQueries=false){
      RWMol *res; 
      setInputCharPtr( (char *)inp.c_str());
      try {
        slnParserDoQueries=doQueries;
        func();
        if(SLNParse::molList_g.size()<=0){
          res = 0;
        } else {
          res = SLNParse::molList_g[0];
          SLNParse::molList_g.resize(0);
        }
      } catch (SLNParseException &e) {
        BOOST_LOG(rdErrorLog) << e.message() << std::endl;
        res = 0;
      }
      if(res){
        // cleanup:
        res->clearAllAtomBookmarks();
        res->clearAllBondBookmarks();

        // set up the chirality flags as soon as the molecule is finished
        // since we'll be removing Hs later and that will break things:
        adjustAtomChiralities(res);

      }
      return res;
    };
  } // end of SLNParse namespace
  
  RWMol *SLNToMol(std::string sln,bool sanitize,int debugParse){
    // FIX: figure out how to reset lexer state
    yysln_debug = debugParse;
    // strip any leading/trailing whitespace:
    boost::trim_if(sln,boost::is_any_of(" \t\r\n"));

    RWMol *res = SLNParse::toMol(sln,yysln_parse,false);
    if(res){
      for(ROMol::AtomIterator atomIt=res->beginAtoms();
          atomIt!=res->endAtoms();++atomIt){
        SLNParse::parseFinalAtomAttribs(*atomIt,false);
      }
      if(sanitize){
        // we're going to remove explicit Hs from the graph,
        // this triggers a sanitization, so we do not need to
        // worry about doing one here:
        ROMol *tmp = MolOps::removeHs(*res,false,false);
        delete res;
        res = static_cast<RWMol *>(tmp);
      }
    }
    charPtrCleanup();
    return res;
  };

  RWMol *SLNQueryToMol(std::string sln,bool mergeHs,int debugParse){
    yysln_debug = debugParse;
    // strip any leading/trailing whitespace:
    boost::trim_if(sln,boost::is_any_of(" \t\r\n"));
    RWMol *res = SLNParse::toMol(sln,yysln_parse,true);
    if(res){
      RWMol *tmp = SLNParse::finalizeQueryMol(res,mergeHs);
      delete res;
      res=tmp;
    }
    charPtrCleanup();
    return res;
  };
}
