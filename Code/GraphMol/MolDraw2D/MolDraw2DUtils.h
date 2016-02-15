//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef MOLDRAW2DUTILS_H
#define MOLDRAW2DUTILS_H
#include <GraphMol/RWMol.h>

// ****************************************************************************

namespace RDKit {
namespace MolDraw2DUtils {
void prepareMolForDrawing(RWMol &mol, bool kekulize = true,
                          bool addChiralHs = true, bool wedgeBonds = true,
                          bool forceCoords = false);
}
}
#endif  // MOLDRAW2DUTILS_H
