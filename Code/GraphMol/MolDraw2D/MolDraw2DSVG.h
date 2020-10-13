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
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/DrawTextSVG.h>

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2DSVG : public MolDraw2D {
 public:
  // initialize to use a particular ostream
  MolDraw2DSVG(int width, int height, std::ostream &os, int panelWidth = -1,
               int panelHeight = -1, bool noFreetype = false)
      : MolDraw2D(width, height, panelWidth, panelHeight), d_os(os) {
    initDrawing();
    initTextDrawer(noFreetype);
  };

  // initialize to use the internal stringstream
  MolDraw2DSVG(int width, int height, int panelWidth = -1, int panelHeight = -1,
               bool noFreetype = false)
      : MolDraw2D(width, height, panelWidth, panelHeight), d_os(d_ss) {
    initDrawing();
    initTextDrawer(noFreetype);
  };

  void setColour(const DrawColour &col) override;

  // not sure if this goes here or if we should do a dtor since initDrawing() is
  // called in the ctor,
  // but we'll start here
  void finishDrawing();

  void drawLine(const Point2D &cds1, const Point2D &cds2) override;
  void drawPolygon(const std::vector<Point2D> &cds) override;
  void drawEllipse(const Point2D &cds1, const Point2D &cds2) override;
  void clearDrawing() override;

  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16,
                    double vertOffset = 0.05) override;

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

  void initDrawing() override;
  void initTextDrawer(bool noFreetype) override;

 protected:
  void drawBond(const ROMol &mol, const Bond *bond, int at1_idx, int at2_idx,
                const std::vector<int> *highlight_atoms = nullptr,
                const std::map<int, DrawColour> *highlight_atom_map = nullptr,
                const std::vector<int> *highlight_bonds = nullptr,
                const std::map<int, DrawColour> *highlight_bond_map = nullptr,
                const std::vector<std::pair<DrawColour, DrawColour>>
                    *bond_colours = nullptr) override;
  void drawAtomLabel(int atom_num, const DrawColour &draw_colour) override;
  void drawAnnotation(const std::string &note,
                      const std::shared_ptr<StringRect> &note_rect) override;
};

}  // namespace RDKit
#endif  // MOLDRAW2DSVG_H
