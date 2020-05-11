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

   void getStringSize(const std::string &label, double &label_width,
                      double &label_height) const override;
   //! drawString centres the string on cds.
   void drawString(const std::string &str, const Point2D &cds,
                   MolDraw2D::AlignType align) override;
   void drawChar(char c, const Point2D &cds) override;

  protected:
   void alignString(const std::string &str, const Point2D &in_cds,
                    MolDraw2D::AlignType align, Point2D &out_cds) override;

 private:
   std::ostream &oss_;
   std::string &d_active_class_;

};

} // namespace RDKit

#endif  // RDKIT_DRAWTEXTSVG_H
