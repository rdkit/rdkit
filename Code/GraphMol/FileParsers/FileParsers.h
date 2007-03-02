//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
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

  //typedef std::vector<RWMol *> MOLPTR_VECT;
  typedef std::vector< RWMOL_SPTR > RWMOL_SPTR_VECT;
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
			    bool sanitize=true);
  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
			    bool sanitize=true);
  RWMol *MolBlockToMol(const std::string &inStream, bool sanitize=true);
  RWMol *MolFileToMol(std::string fName, bool sanitize=true);
  // mol writing:
  std::string MolToMolBlock(const ROMol &mol,bool includeStereo=true, int confId=-1);
  void MolToMolFile(const ROMol &mol,std::string fName,bool includeStereo=true, int confId=-1);

  // sd reading:
  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream *inStream,bool sanitize=true);
  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream &inStream,bool sanitize=true);
  RWMOL_SPTR_VECT SDFileToMols(std::string fName,bool sanitize=true);

  // TPL handling:
  RWMol *TPLDataStreamToMol(std::istream *inStream, unsigned int &line,
			    bool sanitize=true,
                            bool skipFirstConf=false);
  RWMol *TPLFileToMol(std::string fName,bool sanitize=true,
                      bool skipFirstConf=false);


}

#endif
