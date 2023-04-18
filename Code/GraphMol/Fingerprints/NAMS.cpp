//
//  Copyright (C) 2023 Rocco Moretti and other RDKit contributors
//  Copyright (C) 2013, Andre Falcao and Ana Teixeira, University of Lisbon - LaSIGE
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Portions of this implementation include or are based on NAMS by Andre Falcao and Ana Teixeira 
//  (https://github.com/aofalcao/nams-docker)
//  NAMS is free software: it can be redistributed and/or modifed
//  under the terms of the MIT License as published on the official site of Open Source Initiative
//
// Please cite the authors in any work or product based on this material:
//
// AL Teixeira, AO Falcao. 2013. A non-contiguous atom matching structural similarity function. J. Chem. Inf. Model. DOI: 10.1021/ci400324u.
// (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

#include <GraphMol/Fingerprints/NAMS.h>
#include <GraphMol/ROMol.h>

#include <iostream>
#include <sstream>

namespace RDKit {
namespace NAMS {

// This should dump information to cout in approx
std::string
NAMSMolInfo::dump(int cid) const {
  std::stringstream ss;

  ss << cid << ' ' << int(molwt*10) << ' ' << smiles << '\n';
  ss << natoms << ' ' << nbonds << ' ' << abas.size() << '\n';
  for ( auto const & level: abas ) {
    for ( int x: level ) {
      ss << x << ' ';
    }
    ss << '\n';
  }
  for ( auto const & level: mat_aba_typs ) {
    for ( int x: level ) {
      ss << x << ' ';
    }
    ss << '\n';
  }
  for ( auto const & level: mat_levels ) {
    for ( int x: level ) {
      ss << x << ' ';
    }
    ss << '\n';
  }

  return ss.str();
}

NAMSMolInfo * getNAMSMolInfo(const ROMol &mol) {
  return new NAMSMolInfo; 
}

double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2) {
  // For debugging right now
  std::cerr << "--%%----------------------------------\n";
  std::cerr << molinfo1.dump(1);
  std::cerr << molinfo2.dump(2);
  std::cerr << "--%%----------------------------------\n";
  return 0;
}

} // namespace NAMS
} // namespace RDKit
