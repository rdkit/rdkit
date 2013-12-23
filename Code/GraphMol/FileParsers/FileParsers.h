//
//  Copyright (C) 2002-2013 Greg Landrum, Rational Discovery LLC
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
   *   \param forceV3000    - force generation a V3000 mol block (happens automatically with
   *                          more than 999 atoms or bonds)
   */
  std::string MolToMolBlock(const ROMol &mol,bool includeStereo=true,
                            int confId=-1,bool kekulize=true,bool forceV3000=false);
  // \brief Writes a molecule to an MDL mol file
  /*! 
   *   \param mol           - the molecule in question
   *   \param fName         - the name of the file to use
   *   \param includeStereo - toggles inclusion of stereochemistry information
   *   \param confId        - selects the conformer to be used
   *   \param kekulize      - triggers kekulization of the molecule before it is written
   *   \param forceV3000    - force generation a V3000 mol block (happens automatically with
   *                          more than 999 atoms or bonds)
   */
  void MolToMolFile(const ROMol &mol,std::string fName,bool includeStereo=true,
                    int confId=-1,bool kekulize=true,bool forceV3000=false);


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

  RWMol *PDBBlockToMol(const char *str, bool sanitize=true,
                       bool removeHs=true, unsigned int flavor=0);

  RWMol *PDBBlockToMol(const std::string &str, bool sanitize=true,
                       bool removeHs=true, unsigned int flavor=0);
  RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize=true,
                            bool removeHs=true, unsigned int flavor=0);
  RWMol *PDBDataStreamToMol(std::istream &inStream, bool sanitize=true,
                            bool removeHs=true, unsigned int flavor=0);
  RWMol *PDBFileToMol(const std::string &fname, bool sanitize=true,
                      bool removeHs=true, unsigned int flavor=0);

  // \brief generates an PDB block for a molecule
  /*! 
   *   \param mol           - the molecule in question
   *   \param confId        - selects the conformer to be used
   *   \param flavor        - controls what gets written:
   *         flavor & 1 : Write MODEL/ENDMDL lines around each record
   *         flavor & 2 : Don't write any CONECT records
   *         flavor & 4 : Write CONECT records in both directions
   *         flavor & 8 : Don't use multiple CONECTs to encode bond order
   *         flavor & 16 : Write MASTER record
   *         flavor & 32 : Write TER record
   */
  std::string MolToPDBBlock(const ROMol &mol, int confId=-1, unsigned int flavor=0);
  // \brief Writes a molecule to an MDL mol file
  /*! 
   *   \param mol           - the molecule in question
   *   \param fName         - the name of the file to use
   *   \param confId        - selects the conformer to be used
   *   \param flavor        - controls what gets written:
   *         flavor & 1 : Write MODEL/ENDMDL lines around each record
   *         flavor & 2 : Don't write any CONECT records
   *         flavor & 4 : Write CONECT records in both directions
   *         flavor & 8 : Don't use multiple CONECTs to encode bond order
   *         flavor & 16 : Write MASTER record
   *         flavor & 32 : Write TER record
   */
  void MolToPDBFile(const ROMol &mol,const std::string &fname, int confId=-1, unsigned int flavor=0);
}

#endif
