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

#include <sstream>
#include <boost/algorithm/string.hpp>

#include <GraphMol/MolDraw2D/DrawTextSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

using namespace std;

namespace RDKit {

string DrawColourToSVG(const RDKit::DrawColour &col);

// ****************************************************************************
DrawTextSVG::DrawTextSVG(double max_fnt_sz, ostream &oss, string &d_act_class) :
    DrawText(max_fnt_sz), oss_(oss), d_active_class_(d_act_class) {
  cout << "DrawTextSVG" << endl;
}

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
void DrawTextSVG::drawChar(char c, const Point2D &cds) {

  unsigned int fontSz = fontSize();
  std::string col = DrawColourToSVG(colour());

  oss_ << "<text";
  oss_ << " x='" << cds.x;
  oss_ << "' y='" << cds.y << "'";
  if (!d_active_class_.empty()) {
    oss_ << " class='" << d_active_class_ << "'";
  }
  oss_ << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  oss_ << " >";
  oss_ << c;
  oss_ << "</text>" << endl;

}

// ****************************************************************************
void DrawTextSVG::getStringRects(const string &text,
                                 vector<shared_ptr<StringRect> > &rects,
                                 vector<TextDrawType> &draw_modes,
                                 vector<char> &draw_chars) const {

  cout << "DrawTextSVG::getStringRects" << endl;

  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  double running_x = 0.0;
  double act_font_size = fontSize();
  double char_height;
  for (size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_chars.push_back(text[i]);

    double char_width =
        act_font_size *
        static_cast<double>(MolDraw2D_detail::char_widths[(int)text[i]]) /
        MolDraw2D_detail::char_widths[(int)'W'];
    // Absent a proper set of font metrics (we don't know what font we'll be
    // using, for starters) this is something of an empirical bodge.
    if(text[i] == '+') {
      char_height = 0.6 * act_font_size;
    } else if(text[i] == '-') {
      char_height = 0.4 * act_font_size;
    } else {
      char_height = 0.8 * act_font_size;
    }
    double cscale = selectScaleFactor(text[i], draw_mode);
    char_height *= cscale;
    char_width *= cscale;
    cout << draw_chars.back() << " : " << char_width << " : " << char_height
         << " and " << running_x << "  cscale = " << cscale << endl;
    Point2D offset(char_width / 2, char_height / 2);
    if(text[i] == '+' || text[i] == '-') {
      offset.y /= 2.0;
    }
    Point2D g_centre(char_width / 2, char_height / 2);
    rects.push_back(shared_ptr<StringRect>(new StringRect(offset, g_centre, char_width, char_height)));
    rects.back()->trans_.x += running_x;
    draw_modes.push_back(draw_mode);
    cout << "SVG rect : " << text[i] << " : " << rects.back()->trans_
         << " :: " << rects.back()->width_ << " by " << rects.back()->height_
         << "  offset = " << rects.back()->offset_ << endl;

    running_x += char_width;
  }
  for(auto r: rects) {
    r->g_centre_.y = act_font_size - r->g_centre_.y;
    r->offset_.y = act_font_size / 2.0;
  }

  adjustStringRectsForSuperSubScript(draw_modes, draw_chars, rects);

}

} // namespace RDKit
