//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_DRAWTEXTQT_H
#define RDKIT_DRAWTEXTQT_H

#include <RDGeneral/export.h>
#include <GraphMol/MolDraw2D/DrawTextNotFT.h>

class QPainter;

namespace RDKit {

namespace MolDraw2D_detail {

// ****************************************************************************
class RDKIT_MOLDRAW2DQT_EXPORT DrawTextQt : public DrawTextNotFT {
 public:
  DrawTextQt(double max_fnt_sz, double min_fnt_sz, QPainter *qp);
  DrawTextQt(const DrawTextQt &rhs) = delete;
  DrawTextQt(DrawTextQt &&rhs) = delete;
  DrawTextQt &operator=(const DrawTextQt &rhs) = delete;
  DrawTextQt &operator=(DrawTextQt &&rhs) = delete;

#if 0
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const override;
#endif
  void drawChar(char c, const Point2D &cds) override;

 private:
  QPainter *d_qp;

  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes,
                      std::vector<char> &draw_chars) const override;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit
#endif  // RDKIT_DRAWTEXTQT_H
