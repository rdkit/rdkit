//
//  Copyright (C) 2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx).
//
// This base class redefines alignString for the non-FreeType
// DrawText classes.

#ifndef RDKIT_DRAWTEXTNOTFT_H
#define RDKIT_DRAWTEXTNOTFT_H

#include <RDGeneral/export.h>
#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {
namespace MolDraw2D_detail {

class RDKIT_MOLDRAW2D_EXPORT DrawTextNotFT : public DrawText {
 public:
  DrawTextNotFT(double max_fnt_sz, double min_fnt_sz);
  virtual ~DrawTextNotFT();

  void alignString(
      TextAlignType align, const std::vector<TextDrawType> &draw_modes,
      std::vector<std::shared_ptr<StringRect>> &rects) const override;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTNOTFT_H
