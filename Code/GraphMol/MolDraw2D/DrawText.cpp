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
  cout << "DrawText::drawString 2 : " << label << " at " << cds.x << ", " << cds.y
       << " orient = " << orient << " fontScale = " << fontScale() << endl;

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> draw_chars;
  getStringRects(label, orient, rects, draw_modes, draw_chars);
  drawRects(cds, rects, draw_modes, draw_chars);

}

// ****************************************************************************
void DrawText::adjustLineForString(const string &label, OrientType orient,
                                   const Point2D &end1, Point2D &end2) const {

  cout << "Intersection between " << label << " : " << orient
       << " and " << end1 << " and " << end2 << endl;

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> draw_chars;
  Point2D lab_pos = end2;

  getStringRects(label, orient, rects, draw_modes, draw_chars);

  for(auto r: rects) {
    r->trans_.x += lab_pos.x;
    r->trans_.y += lab_pos.y;

    Point2D tl, tr, bl, br;
    r->calcCorners(tl, tr, br, bl);
    unique_ptr<Point2D> ip(new Point2D);

    // if it's a wide label, such as C:7, the bond can intersect
    // more than 1 side of the rectangle, so check them all.
    if (doLinesIntersect(end2, end1, tl, tr, ip.get())) {
      end2 = *ip;
    }
    if (doLinesIntersect(end2, end1, tr, br, ip.get())) {
      end2 = *ip;
    }
    if (doLinesIntersect(end2, end1, br, bl, ip.get())) {
      end2 = *ip;
    }
    if (doLinesIntersect(end2, end1, bl, tl, ip.get())) {
      end2 = *ip;
    }
  }
  cout << "End of Intersection between" << endl;

}

// ****************************************************************************
void DrawText::drawStringRects(const string &label, OrientType orient,
                               const Point2D &cds, MolDraw2D &mol_draw) const {

  cout << "CCCCCCCCCCCCCCCCCCCCC" << endl;
  cout << "drawStringRects" << endl;
  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> draw_chars;

  size_t i = 0;
  getStringRects(label, orient, rects, draw_modes, draw_chars);
  for(auto r: rects) {
    r->trans_.x += cds.x;
    r->trans_.y += cds.y;
//    if(draw_chars[i] == '+') {
//      r->trans_.y -= 28.86;
//    } else if(draw_chars[i] == '-') {
//      r->trans_.y -= 22.58;
//    }
    Point2D tl, tr, br, bl;
    r->calcCorners(tl, tr, br, bl);
    cout << draw_chars[i] << endl;
    cout << "coords : " << cds << endl;
    cout << "trans : " << r->trans_ << endl;
    cout << "g_centre : " << r->g_centre_ << endl;
    cout << "offset : " << r->offset_ << endl;
    cout << "top left : " << tl << " : " << tl - cds << endl;
    cout << "bottom right : " << br << " : " << br - cds << endl;

    tl = mol_draw.getAtomCoords(make_pair(tl.x, tl.y));
    tr = mol_draw.getAtomCoords(make_pair(tr.x, tr.y));
    br = mol_draw.getAtomCoords(make_pair(br.x, br.y));
    bl = mol_draw.getAtomCoords(make_pair(bl.x, bl.y));

    mol_draw.setColour(DrawColour(1.0, 0.0, 0.0));
    mol_draw.drawLine(tl, tr);
    mol_draw.setColour(DrawColour(0.0, 1.0, 0.0));
    mol_draw.drawLine(tr, br);
    mol_draw.setColour(DrawColour(0.0, 0.0, 1.0));
    mol_draw.drawLine(br, bl);
    mol_draw.setColour(DrawColour(0.0, 1.0, 1.0));
    mol_draw.drawLine(bl, tl);
    ++i;
  }

}

// ****************************************************************************
void DrawText::getStringExtremes(const string &label, OrientType orient,
                                 double &x_min, double &y_min,
                                 double &x_max, double &y_max) const {

  if (label.empty()) {
    x_min = x_max = 0.0;
    y_min = y_max = 0.0;
    return;
  }

  x_min = y_min = numeric_limits<double>::max();
  x_max = y_max = -numeric_limits<double>::max();

  vector<shared_ptr<StringRect>> rects;
  vector<TextDrawType> draw_modes;
  vector<char> to_draw;
  getStringRects(label, orient, rects, draw_modes, to_draw);
  cout << "label : " << label << " mins : " << x_min << ", " << y_min
       << "  maxes : " << x_max << ", " << y_max << " fontscale = " << fontScale() << endl;
  for(auto r: rects) {
    Point2D tl, tr, br, bl;
    r->calcCorners(tl, tr, br, bl);
    x_min = min(bl.x, x_min);
    y_min = min(bl.y, y_min);
    x_max = max(tr.x, x_max);
    y_max = max(tr.y, y_max);
    cout << "label : " << label << " mins : " << x_min << ", " << y_min
         << "  maxes : " << x_max << ", " << y_max << " fontscale = " << fontScale() << endl;
  }

  cout << "label : " << label << " mins : " << x_min << ", " << y_min
       << "  maxes : " << x_max << ", " << y_max << " fontscale = " << fontScale() << endl;

}

// ****************************************************************************
void DrawText::alignString(TextAlignType talign,
                           const vector<TextDrawType> &draw_modes,
                           vector<shared_ptr<StringRect> > &rects) const {

  // string comes in with rects aligned with first char with its
  // left hand and bottom edges at 0 on y and x respectively.
  // Adjust relative to that so that the relative alignment point is at
  // (0,0).
  cout << "DrawText::alignString : " << talign << endl;
  cout << "Rects : " << endl;
  for(auto r: rects) {
    cout << r->trans_ << " dims " << r->width_ << " by " << r->height_ << endl;
  }

  if(talign == TextAlignType::MIDDLE) {
    size_t num_norm = count(draw_modes.begin(), draw_modes.end(),
                            TextDrawType::TextDrawNormal);
    if (num_norm == 1) {
      talign = TextAlignType::START;
    }
  }
  cout << "DrawText::alignString amended align : " << talign << endl;

  Point2D align_trans, align_offset;
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
    cout << "Align char : " << align_char << "  offset = " << rects[align_char]->offset_ << endl;
    align_trans = rects[align_char]->trans_;
    align_offset = rects[align_char]->offset_;
  } else {
    // centre on the middle of the Normal text.  The super- or subscripts
    // should be at the ends.
    double full_width = 0.0;
    double y_min = numeric_limits<double>::max();
    double y_max = -y_min;
    align_offset.x = align_offset.y = 0.0;
    int num_norm = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        full_width += rects[i]->width_;
        y_max = max(y_max, rects[i]->trans_.y + rects[i]->height_ / 2);
        y_min = min(y_min, rects[i]->trans_.y - rects[i]->height_ / 2);
        align_offset += rects[i]->offset_;
        ++num_norm;
      }
    }
    align_trans.x = full_width / 2.0;
    align_trans.y = (y_max + y_min) / 2.0;
    align_offset /= num_norm;
  }

  for(auto r: rects) {
    r->trans_ -= align_trans;
    r->offset_ = align_offset;
  }

  for(auto r: rects) {
    cout << "Aligned : " << r->trans_ << " dims " << r->width_
         << " by " << r->height_ << endl;
  }

}

// ****************************************************************************
void DrawText::adjustStringRectsForSuperSubScript(const vector<TextDrawType> &draw_modes,
                                                  const vector<char> &to_draw,
                                                  vector<shared_ptr<StringRect> > &rects) const {

  double last_char = -1;
  for(size_t i = 0; i < draw_modes.size(); ++i) {
    switch(draw_modes[i]) {
      case TextDrawType::TextDrawSuperscript:
        cout << "Adjusting superscript " << to_draw[i] << " : "
             << rects[i]->trans_ << " :: " << rects[i]->width_
             << " by " << rects[i]->height_ << endl;
        // superscripts may come before or after a letter.  If before,
        // spin through to first non-superscript for last_height
        if(last_char < 0) {
          for (size_t j = i + 1; j < draw_modes.size(); ++j) {
            if (draw_modes[j] == TextDrawType::TextDrawNormal) {
              last_char = j;
              break;
            }
          }
        }
        cout << "last_char " << last_char << " : " << to_draw[last_char]
             << " : " << rects[last_char]->trans_
             << " :: "
             << rects[last_char]->width_ << " by " << rects[last_char]->height_ << endl;
        // adjust up by last height / 2.
        rects[i]->rect_corr_ = rects[last_char]->height_;
        rects[i]->trans_.y -= rects[i]->rect_corr_ / 2.0;
        cout << "Corrector : " << rects[i]->rect_corr_ << endl;
        cout << "After adjusting superscript " << to_draw[i] << " : "
             << rects[i]->trans_ << " :: " << rects[i]->width_
             << " by " << rects[i]->height_ << endl;
        // if the last char was a subscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSubscript) {
          double move_by = rects[i]->trans_.x - rects[i-1]->trans_.x;
          if(move_by > 0.0) {
            for (size_t j = 0; j < i; ++j) {
              rects[j]->trans_.x += move_by;
            }
          } else {
            for (size_t j = i; j < draw_modes.size(); ++j) {
              rects[j]->trans_.x += move_by;
            }
          }
          cout << "correct advance for " << to_draw[i] << " : "
               << i << " of " << draw_modes.size() - 1 << endl;
          rects[i]->trans_.x = rects[i-1]->trans_.x;
        }
        break;
      case TextDrawType::TextDrawSubscript:
        cout << "Adjusting superscript " << to_draw[i] << " : "
             << rects[i]->trans_ << " :: " << rects[i]->width_
             << " by " << rects[i]->height_ << endl;
        // adjust y down by last height / 2
        rects[i]->rect_corr_ = -rects[last_char]->height_;
        rects[i]->trans_.y -= rects[i]->rect_corr_ / 2.0;
        // if the last char was a superscript, remove the advance
        if (i && draw_modes[i - 1] == TextDrawType::TextDrawSuperscript) {
          cout << "Correct advance" << endl;
          rects[i]->trans_.x = rects[i-1]->trans_.x;
        }
        break;
      case TextDrawType::TextDrawNormal:
        last_char = i;
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
  double x_min, y_min, x_max, y_max;
  getStringExtremes(label, OrientType::E, x_min, y_min, x_max, y_max);
  label_width = x_max - x_min;
  label_height = y_max - y_min;

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
                              vector<char> &draw_chars) const {

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
      alignString(TextAlignType::MIDDLE, t_draw_modes, t_rects);
      double max_height = -numeric_limits<double>::max();
      for(auto r: t_rects) {
        max_height = max(r->height_, max_height);
        r->trans_.y += running_y;
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
    cout << "raw draw rect " << i << " :: " << rects[i]->trans_ << " : "
        << rects[i]->width_ << " by " << rects[i]->height_ << endl;
    Point2D draw_cds;
    draw_cds.x = a_cds.x + rects[i]->trans_.x - rects[i]->offset_.x;
    draw_cds.y = a_cds.y - rects[i]->trans_.y + rects[i]->offset_.y; // opposite sign convention
    draw_cds.y -= rects[i]->rect_corr_;
//    if(draw_chars[i] == '+') {
//      draw_cds.y -= 30;
//    } else if(draw_chars[i] == '-') {
//      draw_cds.y -= 30;
//    }
    cout << "Drawing char at : " << draw_cds.x << ", " << draw_cds.y << endl;
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

  // if the orientation is E, any charge flag needs to be at the end.
  if(orient == OrientType::E) {
    for (size_t i = 0; i < label_pieces.size(); ++i) {
      if(label_pieces[i] == "<sup>+</sup>" || label_pieces[i] == "<sup>-</sup>") {
        label_pieces.push_back(label_pieces[i]);
        label_pieces[i].clear();
        break;
      }
    }
  }

  // Now group some together.  This relies on the order that
  // getAtomLabel built them in the first place.  Each atom symbol
  // needs to be flanked by any <sub> and <super> pieces.
  vector<string> final_pieces;
  string curr_piece;
  bool had_symbol = false;
  for(const auto p: label_pieces) {
    if(p.empty()) {
      continue;
    }
    if(!isupper(p[0])) {
      curr_piece += p;
    } else {
      if(had_symbol) {
        final_pieces.push_back(curr_piece);
        curr_piece = p;
        had_symbol = true;
      } else {
        curr_piece += p;
        had_symbol = true;
      }
    }
  }
  if(!curr_piece.empty()) {
    final_pieces.push_back(curr_piece);
  }

   cout << "Final pieces : " << endl;
   for(auto l: final_pieces) {
     cout << l << endl;
   }
   cout << endl;

  return final_pieces;

#if 0
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
#endif


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