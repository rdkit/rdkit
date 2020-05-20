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
DrawText::DrawText(double max_fnt_sz)
    : colour_(DrawColour(0.0, 0.0, 0.0)), font_scale_(1.0),
      max_font_size_(max_fnt_sz) {

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
double DrawText::maxFontSize() const {
  return max_font_size_;
}

// ****************************************************************************
void DrawText::setMaxFontSize(double new_max) {
  max_font_size_ = new_max;
}

// ****************************************************************************
double DrawText::fontScale() const {
  return font_scale_;
}

// ****************************************************************************
void DrawText::setFontScale(double new_scale) {

  font_scale_ = new_scale;
  double nfs = fontSize();
  if(nfs > max_font_size_) {
    font_scale_ = max_font_size_ / FONT_SIZE;
  }
  cout << "New font size : " << fontSize() << endl;

}

// ****************************************************************************
void DrawText::drawString(const string &str, const Point2D &cds,
                          TextAlignType talign) {

  cout << "XXXXXXXXXXXXXXXXX" << endl;
  cout << "DrawText::drawString 1 : " << str << " at " << cds.x << ", " << cds.y
       << " talign = " << talign << " fontScale = " << fontScale() << endl;

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> draw_chars;
  getStringRects(str, rects, draw_modes, draw_chars);
  alignString(talign, draw_modes, rects);
  drawRects(cds, rects, draw_modes, draw_chars);

}

// ****************************************************************************
void DrawText::drawString(const string &label, const Point2D &cds,
                          OrientType orient) {

  cout << "XXXXXXXXXXXXXXXXX" << endl;
  cout << "DrawText::drawString 1 : " << label << " at " << cds.x << ", " << cds.y
       << " orient = " << orient << " fontScale = " << fontScale() << endl;

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> draw_chars;
  getStringRects(label, orient, rects, draw_modes, draw_chars);
  drawRects(cds, rects, draw_modes, draw_chars);

}

// ****************************************************************************
StringRect DrawText::getStringRect(const string &label,
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
    vector<char> draw_chars;
    getStringRects(lab, rects, draw_modes, draw_chars);

    TextAlignType talign = TextAlignType::MIDDLE;
    if(orient == OrientType::E) {
      talign = TextAlignType::START;
    } else if(orient == OrientType::W) {
      talign = TextAlignType::END;
    }
    alignString(talign, draw_modes, rects);

    t_x_min = t_y_min = numeric_limits<double>::max();
    t_x_max = t_y_max = -numeric_limits<double>::max();
    for(auto rect: rects) {
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
void DrawText::alignString(TextAlignType talign,
                           const vector<TextDrawType> &draw_modes,
                           vector<shared_ptr<StringRect> > &rects) const {

  cout << "DrawText::alignString : " << talign << endl;
  cout << "Rects : " << endl;
  for(auto r: rects) {
    cout << r->centre_ << " dims " << r->width_ << " by " << r->height_ << endl;
  }

  double full_width = 0.0;
  for(const auto r: rects) {
    full_width += r->width_;
  }

  double align_h = 0.0;
  double align_w = 0.0;
  if(talign == TextAlignType::START || talign == TextAlignType::END) {
    size_t align_char = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        align_char = i;
        if (talign == TextAlignType::START) {
          break;
        }
      }
    }
    cout << "Align char : " << align_char << endl;
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
  cout << "Align width : " << align_w << endl;
  cout << "Align height : " << align_h << endl;
  for(auto r: rects) {
    r->centre_.y += align_h / 2;
    switch(talign) {
      case TextAlignType::START: case TextAlignType::MIDDLE:
        r->centre_.x -= align_w / 2;
        break;
      case TextAlignType::END:
        r->centre_.x -= full_width + align_w / 2;
        break;
    }
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
double DrawText::selectScaleFactor(char c, TextDrawType draw_type) const {

  switch(draw_type) {
    case TextDrawType::TextDrawNormal:
      return 1.0;
    case TextDrawType::TextDrawSubscript:
      return SUBS_SCALE;
    case TextDrawType::TextDrawSuperscript:
      if(c == '-' || c == '+') {
        return SUBS_SCALE;
      } else {
        return SUPER_SCALE;
      }
  }
  return 1.0;
}

// ****************************************************************************
void DrawText::getStringSize(const std::string &label, double &label_width,
                             double &label_height) const {

  cout << "DrawText::getStringSize" << endl;
  StringRect rect = getStringRect(label, OrientType::E);
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
void DrawText::getStringRects(const string &text, OrientType orient,
                              vector<shared_ptr<StringRect> > &rects,
                              vector<TextDrawType> &draw_modes,
                              vector<char> &draw_chars) {

  vector<string> text_bits = atomLabelToPieces(text, orient);
  if(orient == OrientType::W) {
    // stick the pieces together again backwards and draw as one so there
    // aren't ugly splits in the string.
    string new_lab;
    for (auto i = text_bits.rbegin(); i != text_bits.rend(); ++i) {
      new_lab += *i;
    }
    getStringRects(new_lab, rects, draw_modes, draw_chars);
    alignString(TextAlignType::END, draw_modes, rects);
  } else if(orient == OrientType::E) {
    // likewise, but forwards
    string new_lab;
    for (auto lab: text_bits) {
      new_lab += lab;
    }
    getStringRects(new_lab, rects, draw_modes, draw_chars);
    alignString(TextAlignType::START, draw_modes, rects);
  } else {
    double running_y = 0;
    for(size_t i = 0; i < text_bits.size(); ++i) {
      vector<shared_ptr<StringRect>> t_rects;
      vector<TextDrawType> t_draw_modes;
      vector<char> t_draw_chars;
      getStringRects(text_bits[i], t_rects, t_draw_modes, t_draw_chars);
      alignString(TextAlignType::START, t_draw_modes, t_rects);
      double max_height = -numeric_limits<double>::max();
      for(auto r: t_rects) {
        max_height = max(r->height_, max_height);
        r->centre_.y += running_y;
      }
      rects.insert(rects.end(), t_rects.begin(), t_rects.end());
      draw_modes.insert(draw_modes.end(), t_draw_modes.begin(),
                        t_draw_modes.end());
      draw_chars.insert(draw_chars.end(), t_draw_chars.begin(),
                        t_draw_chars.end());
      running_y += max_height;
    }
  }

}

// ****************************************************************************
void DrawText::drawRects(const Point2D &a_cds,
                         const std::vector<std::shared_ptr<StringRect> > &rects,
                         const std::vector<TextDrawType> &draw_modes,
                         const std::vector<char> &draw_chars) {

  double full_scale = fontScale();
  for(size_t i = 0; i < rects.size(); ++i) {
    cout << "rect " << i << " :: " << rects[i]->centre_ << " : "
        << rects[i]->width_ << " by " << rects[i]->height_ << endl;
    Point2D draw_cds;
    draw_cds.x = a_cds.x + rects[i]->centre_.x - rects[i]->width_ / 2.0;
    draw_cds.y = a_cds.y + rects[i]->centre_.y - rects[i]->height_ / 2.0;
    setFontScale(full_scale * selectScaleFactor(draw_chars[i], draw_modes[i]));
    drawChar(draw_chars[i], draw_cds);
    setFontScale(full_scale);
  }

}

// ****************************************************************************
vector<string> atomLabelToPieces(const string &label, OrientType orient) {

  cout << "ZZZZZZZZZZ\nsplitting " << label << " : " << orient << endl;
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
            (label_pieces[j][5] == '+' || label_pieces[j][5] == '-')) {
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

  cout << "Final pieces : ";
  for(auto l: label_pieces) {
    cout << l << endl;
  }
  cout << endl;

  return label_pieces;
}


ostream& operator<<(ostream &oss, const TextAlignType &tat) {
  switch(tat) {
    case TextAlignType::START: oss << "START"; break;
    case TextAlignType::MIDDLE: oss << "MIDDLE"; break;
    case TextAlignType::END: oss << "END"; break;
  }
  return oss;
}
ostream& operator<<(ostream &oss, const TextDrawType &tdt) {
  switch(tdt) {
    case TextDrawType::TextDrawNormal: oss << "TextDrawNormal"; break;
    case TextDrawType::TextDrawSuperscript: oss << "TextDrawSuperscript"; break;
    case TextDrawType::TextDrawSubscript: oss << "TextDrawSubscript"; break;
  }
  return oss;
}
ostream& operator<<(ostream &oss, const OrientType &o) {
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