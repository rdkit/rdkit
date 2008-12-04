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
// Created by Greg Landrum, September 2006
//
#ifndef __RD_SLNATTRIBS_H__
#define __RD_SLNATTRIBS_H__

#include <string>
#include <vector>
#include <boost/smart_ptr.hpp>

namespace RDKit{
  class Atom;
  class Bond;

  namespace SLNParse {
    typedef enum {
      AttribLowPriAnd=0,
      AttribOr,
      AttribAnd,
      AttribNot
    } AttribCombineOp;

    class AttribType {
    public:
      AttribType() :first(""),second(""),op(""),negated(false),structQuery(0) {};
      std::string first;
      std::string second;
      std::string op;
      bool negated;
      void *structQuery;
    };

    typedef std::vector< std::pair<AttribCombineOp,boost::shared_ptr<AttribType> > > AttribListType;

    //! parses the attributes provided for an atom and sets
    // the appropriate RD properties/queries.
    // NOTES: 
    //    1) Some SLN query values cannot be properly set until the molecule is fully/
    //       initialized. These are handled by parseFinalAtomAttribs()
    // 
    void parseAtomAttribs(Atom *atom,AttribListType attribs,bool doingQuery);
    void parseFinalAtomAttribs(Atom *atom,bool doingQuery);

    //! parses the attributes provided for a bond and sets
    // the appropriate RD properties/queries.
    // NOTES: 
    //    1) Some SLN query values cannot be properly set until the molecule is fully/
    //       initialized. These are handled by parseFinalBondAttribs()
    void parseBondAttribs(Bond *bond,AttribListType attribs,bool doingQuery);
    void parseFinalBondAttribs(Bond *bond,bool doingQuery);

    //! parses the attributes provided for a ctab and sets
    // the appropriate RD properties/queries.
    void parseMolAttribs(ROMol *mol,AttribListType attribs);
    
    void adjustAtomChiralities(RWMol *mol);
  }
}
#endif
