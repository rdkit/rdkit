//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 08/05/2020.
//

#ifndef RDKIT_DRAWTEXTFTQT_H
#define RDKIT_DRAWTEXTFTQT_H

#include <GraphMol/MolDraw2D/DrawTextFT.h>
#include <GraphMol/MolDraw2D/DrawTextQt.h>
class QPainter;
class QPainterPath;

namespace RDKit {

// ****************************************************************************

class DrawTextFTQt : public DrawTextFT {
 public:
  DrawTextFTQt(double max_fnt_sz, double min_fnt_sz,
               const std::string &font_file, QPainter &qp);

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
  QPainter &d_qp;
  std::unique_ptr<QPainterPath> dp_qpp;
};

}  // namespace RDKit

#endif  // RDKIT_DRAWTEXTFTQT_H
