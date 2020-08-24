//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: Greg Landrum
// A concrete class derived from DrawText that uses the JS Canvas
// to draw text onto a picture.
#ifdef __EMSCRIPTEN__
#include <sstream>
#include <boost/algorithm/string.hpp>

#include <GraphMol/MolDraw2D/DrawTextJS.h>
#include <GraphMol/MolDraw2D/MolDraw2DJS.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

namespace RDKit {

std::string DrawColourToSVG(const RDKit::DrawColour &col);

// ****************************************************************************
DrawTextJS::DrawTextJS(double max_fnt_sz, double min_fnt_sz,
                       emscripten::val &canvas)
    : DrawText(max_fnt_sz, min_fnt_sz), canvas_(canvas) {}

namespace {
void escape_xhtml(std::string &data) {
  boost::algorithm::replace_all(data, "&", "&amp;");
  boost::algorithm::replace_all(data, "\"", "&quot;");
  boost::algorithm::replace_all(data, "\'", "&apos;");
  boost::algorithm::replace_all(data, "<", "&lt;");
  boost::algorithm::replace_all(data, ">", "&gt;");
}
}  // namespace

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextJS::drawChar(char c, const Point2D &cds) {}

// ****************************************************************************
void DrawTextJS::getStringRects(const std::string &text,
                                std::vector<std::shared_ptr<StringRect>> &rects,
                                std::vector<TextDrawType> &draw_modes,
                                std::vector<char> &draw_chars) const {}

}  // namespace RDKit
#endif  // __EMSCRIPTEN__
