//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: Greg Landrum
//

#ifndef RDKIT_DRAWTEXTFTJS_H
#define RDKIT_DRAWTEXTFTJS_H

#include <iosfwd>
#include <emscripten.h>
#include <emscripten/val.h>
#include <GraphMol/MolDraw2D/DrawTextFT.h>

namespace RDKit {

// ****************************************************************************
class DrawTextFTJS : public DrawTextFT {
 public:
  DrawTextFTJS(double max_fnt_sz, double min_fnt_sz,
               const std::string &font_file, emscripten::val &context);

  int MoveToFunctionImpl(const FT_Vector *to) override;
  int LineToFunctionImpl(const FT_Vector *to) override;
  int ConicToFunctionImpl(const FT_Vector *control,
                          const FT_Vector *to) override;
  int CubicToFunctionImpl(const FT_Vector *controlOne,
                          const FT_Vector *controlTwo,
                          const FT_Vector *to) override;

 protected:
  // adds x_trans_ and y_trans_ to coords returns x advance distance
  virtual double extractOutline() override;

 private:
  emscripten::val &context_;
};

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTFTSVG_H
