//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_BCUT_H
#define RDKIT_BCUT_H
#ifdef RDK_HAS_EIGEN3
#include <RDGeneral/export.h>
#include <vector>
#include <string>
namespace RDKit {
  class ROMol;
  namespace Descriptors {
    const std::string BCUT2DVersion = "1.0.0";
    RDKIT_DESCRIPTORS_EXPORT
    std::pair<double,double> BCUT2D(const ROMol &m, const std::vector<double> &atom_props);
    RDKIT_DESCRIPTORS_EXPORT
    std::pair<double,double> BCUT2D(const ROMol &m, const std::string &atom_double_prop);
    RDKIT_DESCRIPTORS_EXPORT    
    std::vector<double> BCUT2D(const ROMol &m);
  }
}
 
#endif
#endif
