//
//  Copyright (C) 2002-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_FILEPARSERS_H
#define _RD_FILEPARSERS_H

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>

#include <string>
#include <iostream>
#include <vector>
#include <exception>

#include <boost/shared_ptr.hpp>

namespace RDKit{
  const int MOLFILE_MAXLINE=256;
  std::string strip(const std::string &orig);

  //-----
  // mol files
  //-----
  typedef std::vector< RWMOL_SPTR > RWMOL_SPTR_VECT;
  // \brief construct a molecule from MDL mol data in a stream
  /*! 
   *   \param inStream - stream containing the data
   *   \param line     - current line number (used for error reporting)
   *   \param sanitize - toggles sanitization and stereochemistry
   *                     perception of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param line     - current line number (used for error reporting)
   *   \param strictParsing - if not set, the parser is more lax about correctness
   *                          of the contents.
   *                          
   */
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
			    bool sanitize=true,bool removeHs=true,
                            bool strictParsing=true);
  // \overload
  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
			    bool sanitize=true,bool removeHs=true,
                            bool strictParsing=true);
  // \brief construct a molecule from an MDL mol block
  /*! 
   *   \param molBlock - string containing the mol block
   *   \param sanitize - toggles sanitization and stereochemistry
   *                     perception of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param strictParsing - if set, the parser is more lax about correctness
   *                          of the contents.
   */
  RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize=true,
                       bool removeHs=true,bool strictParsing=true);
  
  // \brief construct a molecule from an MDL mol file
  /*! 
   *   \param fName    - string containing the file name
   *   \param sanitize - toggles sanitization and stereochemistry
   *                     perception of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param strictParsing - if set, the parser is more lax about correctness
   *                          of the contents.
   */
  RWMol *MolFileToMol(std::string fName, bool sanitize=true,
                      bool removeHs=true,bool strictParsing=true);

  // \brief generates an MDL mol block for a molecule
  /*! 
   *   \param mol           - the molecule in question
   *   \param includeStereo - toggles inclusion of stereochemistry information
   *   \param confId        - selects the conformer to be used
   *   \param kekulize      - triggers kekulization of the molecule before it is written
   */
  std::string MolToMolBlock(const ROMol &mol,bool includeStereo=true,
                            int confId=-1,bool kekulize=true);
  // \brief construct a molecule from an MDL mol file
  /*! 
   *   \param mol           - the molecule in question
   *   \param fName         - the name of the file to use
   *   \param includeStereo - toggles inclusion of stereochemistry information
   *   \param confId        - selects the conformer to be used
   *   \param kekulize      - triggers kekulization of the molecule before it is written
   */
  void MolToMolFile(const ROMol &mol,std::string fName,bool includeStereo=true,
                    int confId=-1,bool kekulize=true);


  //-----
  //  TPL handling:
  //-----

  //! \brief translate TPL data (BioCad format) into a multi-conf molecule
  /*!
    \param inStream:      the stream from which to read
    \param line:          used to track the line number of errors
    \param sanitize:      toggles sanitization and stereochemistry
                          perception of the molecule
    \param skipFirstConf: according to the TPL format description, the atomic
                          coords in the atom-information block describe the first
			  conformation and the first conf block describes second 
			  conformation. The CombiCode, on the other hand, writes
			  the first conformation data both to the atom-information 
                          block and to the first conf block. We want to be able to
			  read CombiCode-style tpls, so we'll allow this mis-feature
			  to be parsed when this flag is set.
  */
  RWMol *TPLDataStreamToMol(std::istream *inStream, unsigned int &line,
			    bool sanitize=true,
                            bool skipFirstConf=false);

  //! \brief construct a multi-conf molecule from a TPL (BioCad format) file
  /*!
    \param fName:         the name of the file from which to read
    \param sanitize:      toggles sanitization and stereochemistry
                          perception of the molecule
    \param skipFirstConf: according to the TPL format description, the atomic
                          coords in the atom-information block describe the first
			  conformation and the first conf block describes second 
			  conformation. The CombiCode, on the other hand, writes
			  the first conformation data both to the atom-information 
                          block and to the first conf block. We want to be able to
			  read CombiCode-style tpls, so we'll allow this mis-feature
			  to be parsed when this flag is set.
  */
  RWMol *TPLFileToMol(std::string fName,bool sanitize=true,
                      bool skipFirstConf=false);

  std::string MolToTPLText(const ROMol &mol,
			   std::string partialChargeProp="_GasteigerCharge",
			   bool writeFirstConfTwice=false);
  void MolToTPLFile(const ROMol &mol,std::string fName,
		    std::string partialChargeProp="_GasteigerCharge",
		    bool writeFirstConfTwice=false);

  //-----
  //  MOL2 handling
  //-----

  typedef enum {
    CORINA=0  //! supports output from Corina and some dbtranslate output
  } Mol2Type;

  // \brief construct a molecule from a Tripos mol2 file
  /*! 
   *
   *   \param fName    - string containing the file name
   *   \param sanitize - toggles sanitization of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param variant  - the atom type definitions to use
   */
  RWMol *Mol2FileToMol(std::string fName,bool sanitize=true,bool removeHs=true,
                       Mol2Type variant=CORINA);

  // \brief construct a molecule from Tripos mol2 data in a stream
  /*!  
   *   \param inStream - stream containing the data
   *   \param sanitize - toggles sanitization of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param variant  - the atom type definitions to use
   */
  RWMol *Mol2DataStreamToMol(std::istream *inStream,bool sanitize=true,bool removeHs=true,
                             Mol2Type variant=CORINA);
  // \overload 
  RWMol *Mol2DataStreamToMol(std::istream &inStream,bool sanitize=true,bool removeHs=true,
                             Mol2Type variant=CORINA);

  // \brief construct a molecule from a Tripos mol2 block
  /*! 
   *   \param molBlock - string containing the mol block
   *   \param sanitize - toggles sanitization of the molecule
   *   \param removeHs - toggles removal of Hs from the molecule. H removal
   *                     is only done if the molecule is sanitized
   *   \param variant  - the atom type definitions to use
   */
  RWMol *Mol2BlockToMol(const std::string &molBlock,bool sanitize=true,bool removeHs=true,
                        Mol2Type variant=CORINA);

}

#endif
