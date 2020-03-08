//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw a
// molecule into an SVG file

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DSVG_H
#define MOLDRAW2DSVG_H

#include <iostream>
#include <sstream>
#include "MolDraw2D.h"

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2DSVG : public MolDraw2D {
 public:
  // initialize to use a particular ostream
  MolDraw2DSVG(int width, int height, std::ostream &os, int panelWidth = -1,
               int panelHeight = -1)
      : MolDraw2D(width, height, panelWidth, panelHeight), d_os(os) {
    initDrawing();
  };
  // initialize to use the internal stringstream
  MolDraw2DSVG(int width, int height, int panelWidth = -1, int panelHeight = -1)
      : MolDraw2D(width, height, panelWidth, panelHeight), d_os(d_ss) {
    initDrawing();
  };

  void setColour(const DrawColour &col) override;

  // not sure if this goes here or if we should do a dtor since initDrawing() is
  // called in the ctor,
  // but we'll start here
  void finishDrawing();

  void drawLine(const Point2D &cds1, const Point2D &cds2) override;
  void drawString(const std::string &str, const Point2D &cds) override;
  void drawString(const std::string &str, const Point2D &cds,
                  AlignType align) override;
  void alignString(const std::string &str,
                   const std::string &align_char, int align,
                   const Point2D &in_cds, Point2D &out_cds) const override;

  void drawPolygon(const std::vector<Point2D> &cds) override;
  void drawEllipse(const Point2D &cds1, const Point2D &cds2) override;
  void clearDrawing() override;

  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16,
                    double vertOffset = 0.05) override;

  // using the current scale, work out the size of the label in molecule
  // coordinates
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const override;

  // this only makes sense if the object was initialized without a stream
  std::string getDrawingText() const { return d_ss.str(); };

  // adds additional tags to the atoms and bonds in the SVG. This should be
  // invoked *after* the molecule has been drawn
  using MolDraw2D::tagAtoms;  // Avoid overload warning.
  void tagAtoms(const ROMol &mol) override { tagAtoms(mol, 0.2); }
  void tagAtoms(const ROMol &mol, double radius,
                const std::map<std::string, std::string> &events = {});

  // adds metadata describing the molecule to the SVG. This allows
  // molecules to be re-built from SVG with MolFromSVG
  void addMoleculeMetadata(const ROMol &mol, int confId = -1) const;
  void addMoleculeMetadata(const std::vector<ROMol *> &mols,
                           const std::vector<int> confIds = {}) const;

 private:
  std::ostream &d_os;
  std::stringstream d_ss;
  std::string d_activeClass;

  void drawChar(char c, const Point2D &cds) override;
  void initDrawing();

 protected:
  void drawBond(
      const ROMol &mol, const Bond *bond, int at1_idx, int at2_idx,
      const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::vector<int> *highlight_bonds = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::vector<std::pair<DrawColour, DrawColour> > *bond_colours = nullptr) override;
};
}  // namespace RDKit
#endif  // MOLDRAW2DSVG_H
