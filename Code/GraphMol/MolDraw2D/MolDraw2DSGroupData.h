//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKITMOLDRAW2DSGROUPDATA_H
#define RDKITMOLDRAW2DSGROUPDATA_H

#include <string>
#include <vector>

#include <Geometry/point.h>

namespace RDKit {

class ROMol;

namespace MolDraw2D_detail {

//! Holds the text and position of a DAT SGroup label.
struct SGroupDataLabel {
  std::string text;     ///< label text
  RDGeom::Point2D pos;  ///< position in molecule coordinates
  bool positioned;      ///< true if pos came from FIELDDISP; false if pos is
                        ///< the associated atom's conformer position (fallback)
  int atomIdx;          ///< index of the associated atom (-1 if none)
};

//! Returns the text and positions of DAT SGroup labels for a molecule,
//! using the same placement logic as the drawing code.
/*!
  \param mol    the molecule
  \param rotate optional rotation angle in degrees (default 0.0)
  \return a vector of SGroupDataLabel objects, one per rendered DAT SGroup
*/
RDKIT_MOLDRAW2D_EXPORT std::vector<SGroupDataLabel> getSGroupDataLabels(
    const ROMol &mol, double rotate = 0.0);

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif
