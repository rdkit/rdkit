//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_PDBSUPPLIER_H_
#define _RD_PDBSUPPLIER_H_
#include <GraphMol/FileParsers/MolSupplier.h>
#include <string>
#include <iostream>

namespace RDKit {

  //! lazy file parser for PDB files
  class PDBMolSupplier : public MolSupplier {
  public:
    explicit PDBMolSupplier(std::istream *inStream, bool takeOwnership=true,
                            bool sanitize=true, bool removeHs=true,
                            unsigned int flavor=0);
    explicit PDBMolSupplier(const std::string &fname,
                            bool sanitize=true, bool removeHs=true,
                            unsigned int flavor=0);

    virtual ~PDBMolSupplier() {
      if (df_owner && dp_inStream)
        delete dp_inStream;
    };

    virtual void init();
    virtual void reset();
    virtual ROMol *next();
    virtual bool atEnd();

  protected:
    bool df_sanitize,df_removeHs;
    unsigned int d_flavor;
  };
}

#endif  // _RD_PDBSUPPLIER_H_

