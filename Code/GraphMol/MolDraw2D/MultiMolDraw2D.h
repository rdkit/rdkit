//
//   Copyright (C) 2016 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKITMULTIMOLDRAW2D_H
#define RDKITMULTIMOLDRAW2D_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

namespace RDKit {

template <typename T>
class MultiMolDraw2D {
 public:
  //!
  MultiMolDraw2D(unsigned int nRows, unsigned int nCols, int width, int height,
                 bool globalScaling = true);
  virtual ~MultiMolDraw2D() {}
  virtual void drawMolecules(
      const std::vector<ROMOL_SPTR> &mols,
      const std::vector<std::string> *legends = NULL,
      const std::vector<std::vector<int> > *highlight_atoms = NULL,
      const std::vector<std::vector<int> > *highlight_bonds = NULL,
      const std::vector<std::map<int, DrawColour> > *highlight_atom_maps = NULL,
      const std::vector<std::map<int, DrawColour> > *highlight_bond_maps = NULL,
      const std::vector<std::map<int, double> > *highlight_radii = NULL,
      const std::vector<int> *confIds = NULL);

  virtual int width() const { return width_; }
  virtual int height() const { return height_; }
  virtual int nRows() const { return nRows_; }
  virtual int nCols() const { return nCols_; }

  MolDrawOptions &drawOptions() { return options_; }
  const MolDrawOptions &drawOptions() const { return options_; }

 private:
  unsigned int nRows_, nCols_;
  int width_, height_;
  bool globalScaling_;
  MolDrawOptions options_;

  std::vector<std::shared_ptr<T> > drawers_;
};
}

#endif  // RDKITMOLDRAW2D_H
