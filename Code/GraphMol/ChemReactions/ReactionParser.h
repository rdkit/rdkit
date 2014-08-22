//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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

#ifndef __RD_REACTIONPARSER_H_21Aug2006__
#define __RD_REACTIONPARSER_H_21Aug2006__

#include <string>
#include <iostream>

namespace RDKit{
  class ChemicalReaction;
  
  //! used to indicate an error in parsing reaction data
  class ChemicalReactionParserException : public std::exception {
  public:
    //! construct with an error message
    explicit ChemicalReactionParserException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    explicit ChemicalReactionParserException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~ChemicalReactionParserException () throw () {};
  private:
    std::string _msg;
  };
  
  
  
  //! Parse a text block in MDL rxn format into a ChemicalReaction 
  ChemicalReaction * RxnBlockToChemicalReaction(const std::string &rxnBlock);
  //! Parse a file in MDL rxn format into a ChemicalReaction 
  ChemicalReaction * RxnFileToChemicalReaction(const std::string &fileName);
  //! Parse a text stream in MDL rxn format into a ChemicalReaction 
  ChemicalReaction * RxnDataStreamToChemicalReaction(std::istream &rxnStream,
                                                     unsigned int &line);
  
  //! Parse a string containing "Reaction SMARTS" into a ChemicalReaction
  /*!
     Our definition of Reaction SMARTS is something that looks a lot like
     reaction SMILES, except that SMARTS queries are allowed on the reactant
     side and that atom-map numbers are required (at least for now)

     \param text          the SMARTS to convert
     \param replacements  a string->string map of replacement strings.
                          \see SmilesToMol for more information about replacements
     \param useSmiles     if set, the SMILES parser will be used instead of the SMARTS
                           parserfor the individual components 
   */ 
  ChemicalReaction * RxnSmartsToChemicalReaction(const std::string &text,
                                                 std::map<std::string,std::string> *replacements=0,
                                                 bool useSmiles=false);

  //! Parse a ROMol into a ChemicalReaction, RXN role must be set before
  /*!
     Alternative to build a reaction from a molecule (fragments) which have RXN roles
     set as atom properties: "molRxnRole" (1=reactant, 2=product, 3=agent) 

     \param mol           ROMol with RXN roles set
   */
  ChemicalReaction * RxnMolToChemicalReaction(const ROMol &mol);



  //! returns the reaction SMARTS for a reaction
  std::string ChemicalReactionToRxnSmarts(ChemicalReaction &rxn);

  //! returns the reaction SMILES for a reaction
  std::string ChemicalReactionToRxnSmiles(ChemicalReaction &rxn);

  //! returns an RXN block for a reaction
  /*!
     \param rxn            chemical reaction
     \param separateAgents flag to decide if agents were put in a seperate block,
                           otherwise they were included in the reactants block (default)
   */
  std::string ChemicalReactionToRxnBlock(const ChemicalReaction &rxn, bool separateAgents=false);
  
  //! returns a ROMol with RXN roles used to describe the reaction
  ROMol* ChemicalReactionToRxnMol(const ChemicalReaction &rxn);

}; // end of RDKit namespace

#endif
