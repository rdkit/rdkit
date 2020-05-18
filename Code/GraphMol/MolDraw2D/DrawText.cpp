//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx) on 29/04/2020.
//


#include <GraphMol/MolDraw2D/DrawText.h>

using namespace std;

namespace RDKit {

// ****************************************************************************
DrawText::DrawText() : font_scale_(1.0) {

  cout << "DrawText : " << FONT_SIZE << endl;

}

// ****************************************************************************
DrawColour const &DrawText::colour() const {
  return colour_;
}

// ****************************************************************************
void DrawText::setColour(const DrawColour &col) {
  colour_ = col;
}

// ****************************************************************************
double DrawText::fontSize() const {
  return fontScale() * FONT_SIZE;
}

// ****************************************************************************
double DrawText::fontScale() const {
  return font_scale_;
}

// ****************************************************************************
void DrawText::setFontScale(double new_scale) { font_scale_ = new_scale; }

// ****************************************************************************
void DrawText::drawString(const string &str, const Point2D &cds,
                          TextAlignType align) {

  cout << "XXXXXXXXXXXXXXXXX" << endl;
  cout << "DrawText : " << str << " at " << cds.x << ", " << cds.y
       << " align = " << align << " fontScale = " << fontScale() << endl;

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  getStringRects(str, rects, draw_modes);

  Point2D a_cds;
  alignString(cds, align, rects, draw_modes, a_cds);

  double full_scale = fontScale();
  double subs_scale = SUBS_SCALE * fontScale();
  double super_scale = SUPER_SCALE * fontScale();
  for(size_t i = 0, j = 0; i < str.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_modes[j], i)) {
      continue;
    }

    switch(draw_modes[j]) {
      case TextDrawType::TextDrawNormal:
        setFontScale(full_scale); break;
      case TextDrawType::TextDrawSubscript:
        setFontScale(subs_scale); break;
      case TextDrawType::TextDrawSuperscript:
        if(str[i] == '-' || str[i] == '+') {
          setFontScale(subs_scale);
        } else {
          setFontScale(super_scale);
        }
        break;
    }

    Point2D draw_cds;
    Point2D tl, tr, br, bl;
    rects[j]->calcCorners(tl, tr, br, bl);
    cout << str[i] << " : " << tl << " to " << br << endl;
    draw_cds.x = a_cds.x + rects[j]->centre_.x - 0.5 * rects[j]->width_;
    draw_cds.y = a_cds.y + rects[j]->centre_.y - 0.5 * rects[j]->height_;
    drawChar(str[i], draw_cds);
    ++j;
  }
  setFontScale(full_scale);

}

// ****************************************************************************
StringRect DrawText::getStringRect(const std::string &label,
                                   OrientType orient) const {

  StringRect string_rect;
  if (label.empty()) {
    return string_rect;
  }
  double x_min, y_min, x_max, y_max;
  x_min = y_min = numeric_limits<double>::max();
  x_max = y_max = -numeric_limits<double>::max();
  double running_y = 0;
  vector<string> label_pieces = atomLabelToPieces(label, orient);
  for (auto lab : label_pieces) {
    double t_x_min, t_y_min, t_x_max, t_y_max;
    vector<shared_ptr<StringRect>> rects;
    vector<TextDrawType> draw_modes; // not needed for bounding box
    getStringRects(lab, rects, draw_modes);

    TextAlignType align = TextAlignType::MIDDLE;
    if(orient == OrientType::E) {
      align = TextAlignType::START;
    } else if(orient == OrientType::W) {
      align = TextAlignType::END;
    }
    Point2D a_cds;
    alignString(Point2D(0,0), align, rects, draw_modes, a_cds);

    t_x_min = t_y_min = numeric_limits<double>::max();
    t_x_max = t_y_max = -numeric_limits<double>::max();
    for(auto rect: rects) {
      rect->centre_ -= a_cds;
      cout << "rect : " << rect->centre_.x << ", " << rect->centre_.y
           << " dims " << rect->width_ << " by " << rect->height_ << endl;
      x_min = min(x_min, rect->centre_.x - rect->width_ / 2.0);
      y_min = min(y_min, rect->centre_.y - rect->height_ / 2.0);
      x_max = max(x_max, rect->centre_.x + rect->width_ / 2.0);
      y_max = max(y_max, rect->centre_.y + rect->height_ / 2.0);
    }

    t_y_max += running_y;
    x_min = min(x_min, t_x_min);
    y_min = min(y_min, t_y_min);
    x_max = max(x_max, t_x_max);
    y_max = max(y_max, t_y_max);
  }

  string_rect.width_ = x_max - x_min;
  string_rect.height_ = y_max - y_min;
  string_rect.centre_.x = string_rect.width_ / 2.0;
  string_rect.centre_.y = string_rect.height_ / 2.0;

  cout << "label : " << label << " at " << string_rect.centre_
       << " dims = " << string_rect.width_ << " by "
       << string_rect.height_ << endl;

  return string_rect;

}


// ****************************************************************************
void DrawText::alignString(const Point2D &in_cds, TextAlignType align,
                           const vector<shared_ptr<StringRect> > &rects,
                           const vector<TextDrawType> &draw_modes,
                           Point2D &out_cds) const {

  cout << "DrawText::alignString" << endl;

  double full_width = 0.0;
  for(const auto r: rects) {
    full_width += r->width_;
  }

  double align_h = 0.0;
  double align_w = 0.0;
  if(align == TextAlignType::START || align == TextAlignType::END) {
    size_t align_char = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        align_char = i;
        if (align == TextAlignType::START) {
          break;
        }
      }
    }
//    cout << "Align char : " << align_char << endl;
    align_h = rects[align_char]->height_;
    align_w = rects[align_char]->width_;
  } else {
    double y_min = numeric_limits<double>::max();
    double y_max = -y_min;
    for(const auto r: rects) {
      y_max = max(y_max, r->centre_.y + r->height_ / 2);
      y_min = min(y_min, r->centre_.y - r->height_ / 2);
    }
    align_h = y_max - y_min;
    align_w = full_width;
  }
//  cout << "Align height : " << align_h << endl;
  out_cds.y = in_cds.y + align_h / 2;

  switch(align) {
    case TextAlignType::START:
      out_cds.x = in_cds.x - align_w / 2;
      break;
    case TextAlignType::MIDDLE:
      out_cds.x = in_cds.x - align_w / 2;
      break;
    case TextAlignType::END:
      out_cds.x = in_cds.x - full_width + align_w / 2;
      break;
  }

}

// ****************************************************************************
void DrawText::adjustStringRectsForSuperSubScript(const vector<TextDrawType> &draw_modes,
                                                  const vector<char> &to_draw,
                                                  vector<shared_ptr<StringRect> > &rects) const {

  double last_height = -1.0;
  for(size_t i = 0; i < draw_modes.size(); ++i) {
    switch(draw_modes[i]) {
      case TextDrawType::TextDrawSuperscript:
        // superscripts may come before or after a letter.  If before,
        // spin through to first non-superscript for last_height
        if(last_height < 0.0) {
          for (size_t j = i + 1; j < draw_modes.size(); ++j) {
            if (draw_modes[j] == TextDrawType::TextDrawNormal) {
              last_height = rects[j]->height_;
              break;
            }
          }
        }
        // adjust y up by last_height / 2, unless it's a + which tends
        // to stick above.
        if(to_draw[i] == '+') {
          rects[i]->centre_.y -= 0.33 * last_height;
        } else {
          rects[i]->centre_.y -= 0.5 * last_height;
        }
        // if the last char was a subscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSubscript) {
          rects[i]->centre_.x = rects[i - 1]->centre_.x;
        }
        break;
      case TextDrawType::TextDrawSubscript:
        // adjust y down by last_height / 2
        rects[i]->centre_.y += 0.5 * last_height;
        // if the last char was a superscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSuperscript) {
          rects[i]->centre_.x = rects[i - 1]->centre_.x;
        }
        break;
      case TextDrawType::TextDrawNormal:
        last_height = rects[i]->height_;
        break;
    }
  }

}

// ****************************************************************************
void DrawText::getStringSize(const std::string &label, double &label_width,
                             double &label_height) const {

  cout << "DrawText::getStringSize" << endl;
  StringRect rect = getStringRect(label, OrientType::W);
  label_width = rect.width_;
  label_height = rect.height_;

}

// ****************************************************************************
bool setStringDrawMode(const string &instring,
                       TextDrawType &draw_mode,
                       size_t &i) {

  string bit1 = instring.substr(i, 5);
  string bit2 = instring.substr(i, 6);

  // could be markup for super- or sub-script
  if (string("<sub>") == bit1) {
    draw_mode = TextDrawType::TextDrawSubscript;
    i += 4;
    return true;
  } else if (string("<sup>") == bit1) {
    draw_mode = TextDrawType::TextDrawSuperscript;
    i += 4;
    return true;
  } else if (string("</sub>") == bit2) {
    draw_mode = TextDrawType::TextDrawNormal;
    i += 5;
    return true;
  } else if (string("</sup>") == bit2) {
    draw_mode = TextDrawType::TextDrawNormal;
    i += 5;
    return true;
  }

  return false;
}

// ****************************************************************************
vector<string> atomLabelToPieces(const string &label, OrientType orient) {

  // cout << "ZZZZZZZZZZ\nsplitting " << label << " : " << orient << endl;
  vector<string> label_pieces;
  if (label.empty()) {
    return label_pieces;
  }

  // if we have the mark-up <lit>XX</lit> the symbol is to be used
  // without modification
  if (label.substr(0, 5) == "<lit>") {
    string lit_sym = label.substr(5);
    size_t idx = lit_sym.find("</lit>");
    if (idx != string::npos) {
      lit_sym = lit_sym.substr(0, idx);
    }
    label_pieces.emplace_back(lit_sym);
    return label_pieces;
  }

  string next_piece;
  size_t i = 0;
  while (true) {
    if (i == label.length()) {
      if (!next_piece.empty()) {
        label_pieces.emplace_back(next_piece);
        break;
      }
    }
    if (label.substr(i, 2) == "<s" || label[i] == ':' || isupper(label[i])) {
      // save the old piece, start a new one
      if (!next_piece.empty()) {
        label_pieces.emplace_back(next_piece);
        next_piece.clear();
      }
    }
    next_piece += label[i++];
  }
  if (label_pieces.size() < 2) {
    return label_pieces;
  }

  // now some re-arrangement to make things look nicer
  // if there's an atom map (:nn) it needs to go after the
  // first atomic symbol which, because of <sup> might not be the first
  // piece.
  for (size_t j = 0; j < label_pieces.size(); ++j) {
    if (label_pieces[j][0] == ':') {
      if (label_pieces[0].substr(0, 5) == "<sup>") {
        label_pieces[1] += label_pieces[j];
      } else {
        label_pieces[0] += label_pieces[j];
      }
      label_pieces[j].clear();
      break;
    }
  }

  // if there's isotope info, and orient is W we want it after the first
  // symbol.  It will be the first piece as getAtomSymbol puts it
  // together.
  if (orient == OrientType::W && label_pieces[0].substr(0, 5) == "<sup>" &&
      isdigit(label_pieces[0][6])) {
    label_pieces[1] = label_pieces[0] + label_pieces[1];
    label_pieces[0].clear();
    label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                       label_pieces.end());
  }

  // if there's a charge, it always needs to be at the end.
  string charge_piece;
  for (size_t j = 0; j < label_pieces.size(); ++j) {
    if (label_pieces[j].substr(0, 6) == "<sup>+" ||
        label_pieces[j].substr(0, 6) == "<sup>-") {
      charge_piece += label_pieces[j];
      label_pieces[j].clear();
    }
  }

  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());
  // if orient is W charge goes to front, otherwise to end.
  if (!charge_piece.empty()) {
    if (orient == OrientType::W) {
      label_pieces.insert(label_pieces.begin(), charge_piece);
    } else {
      label_pieces.emplace_back(charge_piece);
    }
  }

  // if there's a <sub> piece, attach it to the one before.
  for (size_t j = 1; j < label_pieces.size(); ++j) {
    if (label_pieces[j].substr(0, 5) == "<sub>") {
      label_pieces[j - 1] += label_pieces[j];
      label_pieces[j].clear();
      break;
    }
  }
  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());

  // if there's a <sup>[+-.] piece, attach it to the one after.
  if (label_pieces.size() > 1) {
    for (size_t j = 0; j < label_pieces.size() - 1; ++j) {
      if (label_pieces[j].substr(0, 5) == "<sup>") {
        if (orient == OrientType::W &&
            (label_pieces[j][5] == '+' || label_pieces[j][5] == '-' ||
             label_pieces[j][5] == '.')) {
          label_pieces[j + 1] = label_pieces[j + 1] + label_pieces[j];
        } else {
          label_pieces[j + 1] = label_pieces[j] + label_pieces[j + 1];
        }
        label_pieces[j].clear();
        break;
      }
    }
  }
  // and if orient is N or S and the last piece is a charge, attach it to the
  // one before
  if (orient == OrientType::N || orient == OrientType::S) {
    if (label_pieces.back().substr(0, 6) == "<sup>+" ||
        label_pieces.back().substr(0, 6) == "<sup>-") {
      label_pieces[label_pieces.size() - 2] += label_pieces.back();
      label_pieces.back().clear();
    }
  }
  label_pieces.erase(remove(label_pieces.begin(), label_pieces.end(), ""),
                     label_pieces.end());

  // cout << "Final pieces : ";
  // for(auto l: label_pieces) {
  //   cout << l << endl;
  // }
  // cout << endl;

  return label_pieces;
}


std::ostream& operator<<(std::ostream &oss, TextAlignType tat) {
  switch(tat) {
    case TextAlignType::START: oss << "START"; break;
    case TextAlignType::MIDDLE: oss << "MIDDLE"; break;
    case TextAlignType::END: oss << "END"; break;
  }
  return oss;
}
std::ostream& operator<<(std::ostream &oss, TextDrawType tdt) {
  switch(tdt) {
    case TextDrawType::TextDrawNormal: oss << "TextDrawNormal"; break;
    case TextDrawType::TextDrawSuperscript: oss << "TextDrawSuperscript"; break;
    case TextDrawType::TextDrawSubscript: oss << "TextDrawSubscript"; break;
  }
  return oss;
}
std::ostream& operator<<(std::ostream &oss, const OrientType &o) {
  switch(o) {
    case OrientType::C: oss << "C"; break;
    case OrientType::N: oss << "N"; break;
    case OrientType::S: oss << "S"; break;
    case OrientType::E: oss << "E"; break;
    case OrientType::W: oss << "W"; break;
  }
  return oss;
}

} // namespace RDKit