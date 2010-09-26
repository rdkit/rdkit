// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>

#include "Depictor.h"
namespace RDKit {

  //*************************************************************************************
  //
  //  Adds 2D coordinates to a molecule.
  //
  //  ARGUMENTS:
  //    mol:    the molecule to be altered
  //    useDLL: (OPTIONAL) toggles use of the Depict DLL to do the depiction.
  //            This option only makes sense on Windows, where it is the default.
  //
  //  RETURNS:
  //   1 on success, 0 otherwise
  //
  //  NOTES:
  //   - To use the DLL, it's essential that the COMBICHEM_ROOT and COMBICHEM_RELEASE
  //     environment variables be set.  If this isn't done, this whole process
  //     will fail.
  //
  //*************************************************************************************
  int Add2DCoordsToMol(ROMol &mol,bool useDLL){
    int res;
#ifdef WIN32_DLLBUILD
    if(useDLL){
      res = Add2DCoordsToMolDLL(mol);
    } else {
      res = 0;
    }
#else
    res = 0;
#endif
    return res;
  }

  
}
