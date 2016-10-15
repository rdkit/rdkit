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
class MolDraw2D;
namespace MolDraw2DUtils {
//! Does some cleanup operations on the molecule to prepare it to draw nicely
/*
The operations include: kekulization, addition of chiral Hs (so that we can draw
wedges to them), wedging of bonds at chiral centers, and generation of a 2D
conformation if the molecule does not already have a conformation

\param mol: the molecule to be modified
\param kekulize: toggles kekulization (this can fail, see below)
\param addChiralHs: adds Hs to the graph on chiral atoms
\param wedgeBonds: calls WedgeMolBonds()
\param forceCoords: generates a 2D conformation even if one is present already

NOTE: the kekulization step can fail, throwing a MolSanitizeExecption. If this
happens the molecule will be in an inconsistent, partially kekulized, state.
This isn't normally a problem for molecules that have been sanitized, but can be
problematic if the molecules have been modified post santitization.
*/
void prepareMolForDrawing(RWMol &mol, bool kekulize = true,
                          bool addChiralHs = true, bool wedgeBonds = true,
                          bool forceCoords = false);

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const char *json);
void updateDrawerParamsFromJSON(MolDraw2D &drawer, const std::string &json);
}
}
#endif  // MOLDRAW2DUTILS_H
