//
//  Copyright (C) 2001-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SMILESPARSE_H_
#define _RD_SMILESPARSE_H_

#include <string>
#include <exception>
#include <map>

namespace RDKit{
  class RWMol;

  //! Construct a molecule from a SMILES string
  /*!
   \param smi           the SMILES to convert
   \param debugParse    toggles verbose debugging information from the parser
   \param sanitize      toggles H removal and sanitization of the molecule
   \param replacements  a string->string map of replacement strings. See below
                        for more information about replacements.

   \return a pointer to the new molecule; the caller is responsible for free'ing this.

   The optional replacements map can be used to do string substitution of abbreviations
   in the input SMILES. The set of substitutions is repeatedly looped through until
   the string no longer changes. It is the responsiblity of the caller to make sure
   that substitutions results in legal and sensible SMILES.

   Examples of substitutions:
   \code
     CC{Q}C with {"{Q}":"OCCO"} -> CCOCCOC
     C{A}C{Q}C with {"{Q}":"OCCO", "{A}":"C1(CC1)"} -> CC1(CC1)COCCOC
     C{A}C{Q}C with {"{Q}":"{X}CC{X}", "{A}":"C1CC1", "{X}":"N"} -> CC1CC1CCNCCNC
   \endcode

   */
  RWMol *SmilesToMol(std::string smi,int debugParse=0,bool sanitize=1,
                     std::map<std::string,std::string> *replacements=0);
  //! Construct a molecule from a SMARTS string
  /*!
   \param sma           the SMARTS to convert
   \param debugParse    toggles verbose debugging information from the parser
   \param mergeHs       toggles merging H atoms in the SMARTS into neighboring atoms
   \param replacements  a string->string map of replacement strings.
                        \see SmilesToMol for more information about replacements

   \return a pointer to the new molecule; the caller is responsible for free'ing this.
   */
  RWMol *SmartsToMol(std::string sma,int debugParse=0,bool mergeHs=false,
                     std::map<std::string,std::string> *replacements=0);

  class SmilesParseException : public std::exception {
  public:
    SmilesParseException(const char *msg) : _msg(msg) {};
    SmilesParseException(const std::string msg) : _msg(msg) {};
    const char *message () const { return _msg.c_str(); };
    ~SmilesParseException () throw () {};
  private:
    std::string _msg;
  };

}

#endif
