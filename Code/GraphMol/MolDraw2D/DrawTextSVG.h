//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//
// A concrete class derived from DrawText that uses SVG
// to draw text onto a picture.

#ifndef RDKIT_DRAWTEXTSVG_H
#define RDKIT_DRAWTEXTSVG_H

#include <iosfwd>

#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************

class DrawTextSVG : public DrawText {

 public:
   DrawTextSVG(std::ostream &oss, std::string &d_act_class);

#if 0
   void getStringSize(const std::string &label, double &label_width,
                      double &label_height) const override;
#endif
   //! drawString centres the string on cds.
   void drawString(const std::string &str, const Point2D &cds,
                   TextAlignType align) override;
   void drawChar(char c, const Point2D &cds) override;

 protected:
  void alignString(const Point2D &in_cds, TextAlignType align,
                   const std::vector<std::shared_ptr<StringRect> > &rects,
                   const std::vector<TextDrawType> &draw_modes, Point2D &out_cds) const override;

 private:
  std::ostream &oss_;
  std::string &d_active_class_;

  // return a vector of StringRects, one for each char in text, with
  // super- and subscripts taken into account.  Sizes in pixel coords,
  // i.e. scaled by fontScale().
  void getStringRects(const std::string &text,
                      std::vector<std::shared_ptr<StringRect>> &rects,
                      std::vector<TextDrawType> &draw_modes) const override;

};

} // namespace RDKit

#endif  // RDKIT_DRAWTEXTSVG_H
