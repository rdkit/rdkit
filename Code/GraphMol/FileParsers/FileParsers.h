//
//  Copyright (C) 2002-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
			    bool sanitize=true);
  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
			    bool sanitize=true);
  RWMol *MolBlockToMol(const std::string &inStream, bool sanitize=true);
  RWMol *MolFileToMol(std::string fName, bool sanitize=true);

  std::string MolToMolBlock(const ROMol &mol,bool includeStereo=true, int confId=-1);
  void MolToMolFile(const ROMol &mol,std::string fName,bool includeStereo=true, int confId=-1);


  //-----
  // SDF reading:
  //-----
  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream *inStream,bool sanitize=true);
  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream &inStream,bool sanitize=true);
  RWMOL_SPTR_VECT SDFileToMols(std::string fName,bool sanitize=true);

  //-----
  //  TPL handling:
  //-----

  //! \brief translate TPL data (BioCad format) into a multi-conf molecule
  /*!
    \param inStream:      the stream from which to read
    \param line:          used to track the line number of errors
    \param sanitize:      toggles sanitization of the molecule
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
    \param sanitize:      toggles sanitization of the molecule
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



}

#endif
