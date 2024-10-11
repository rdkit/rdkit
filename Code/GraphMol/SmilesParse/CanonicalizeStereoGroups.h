//
//  Copyright (C) 2002-2024  RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_CANONICALIZESTEROGROUPS_H
#define RD_CANONICALIZESTEROGROUPS_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <limits>

namespace RDKit {

enum class StereoGroupAbsOptions {
  OnlyIncludeWhenOtherGroupsExist =
      0,  // when stereo groups are canonicalized, ABS groups are added only
          // when other Stereo Groups are added this prevents a mol from having
          // JUST an ABS stereo group
  NeverInclude,  // Abs groups are never included, Just AND and OR groups after
                 // stereo group canonicalization
  AlwaysInclude  // when stereo groups are canonicalized, ABS groups are added
                 // only when other Stereo Groups are added this prevents a mol
                 // from having JUST an ABS stereo group
};

class RigorousEnhancedStereoException : public std::runtime_error {
 public:
  explicit RigorousEnhancedStereoException(std::string message)
      : std::runtime_error(message) {};
};

RDKIT_SMILESPARSE_EXPORT void canonicalizeStereoGroups(
    std::unique_ptr<ROMol> &mol,
    StereoGroupAbsOptions outputAbsoluteGroups =
        StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist);

}  // namespace RDKit

#endif