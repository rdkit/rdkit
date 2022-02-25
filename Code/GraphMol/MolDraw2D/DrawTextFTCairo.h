//
//  Copyright (C) 2020-2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx).
//

#ifndef RDKIT_DRAWTEXTFTCAIRO_H
#define RDKIT_DRAWTEXTFTCAIRO_H

#include <cairo.h>

#include <GraphMol/MolDraw2D/DrawTextFT.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************

class DrawTextFTCairo : public DrawTextFT {
 public:
  ~DrawTextFTCairo() override = default;
  DrawTextFTCairo(double max_fnt_sz, double min_fnt_sz,
                  const std::string &font_file, cairo_t *dp_cr);
  DrawTextFTCairo(const DrawTextFTCairo &) = delete;
  DrawTextFTCairo(DrawTextFTCairo &&) = delete;
  DrawTextFTCairo &operator=(const DrawTextFTCairo &) = delete;
  DrawTextFTCairo &operator=(DrawTextFTCairo &&) = delete;

  int MoveToFunctionImpl(const FT_Vector *to) override;
  int LineToFunctionImpl(const FT_Vector *to) override;
  int ConicToFunctionImpl(const FT_Vector *control,
                          const FT_Vector *to) override;
  int CubicToFunctionImpl(const FT_Vector *controlOne,
                          const FT_Vector *controlTwo,
                          const FT_Vector *to) override;
  void setCairoContext(cairo_t *cr);

  // adds x_trans_ and y_trans_ to coords returns x advance distance
  double extractOutline() override;

  cairo_t *dp_cr_;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTFTCAIRO_H
