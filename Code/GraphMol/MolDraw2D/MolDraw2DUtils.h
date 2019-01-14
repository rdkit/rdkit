//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DUTILS_H
#define MOLDRAW2DUTILS_H
#include <GraphMol/RWMol.h>

#include <boost/tuple/tuple.hpp>

// ****************************************************************************

namespace RDKit {
typedef boost::tuple<float, float, float> DrawColour;
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
RDKIT_MOLDRAW2D_EXPORT void prepareMolForDrawing(RWMol &mol,
                                                 bool kekulize = true,
                                                 bool addChiralHs = true,
                                                 bool wedgeBonds = true,
                                                 bool forceCoords = false);

//! prepare a molecule for drawing and draw it
/*
  \param mol: the molecule to draw
  \param legend: (optional) the legend (to be drawn under the molecule)
  \param highlight_atoms: (optional) vector of atom ids to highlight
  \param highlight_atoms: (optional) vector of bond ids to highlight
  \param highlight_atom_map: (optional) map from atomId -> DrawColour
            providing the highlight colors. If not provided the default
            higlight colour from \c drawOptions() will be used.
  \param highlight_bond_map: (optional) map from bondId -> DrawColour
            providing the highlight colors. If not provided the default
            highlight colour from \c drawOptions() will be used.
  \param highlight_radii: (optional) map from atomId -> radius (in molecule
            coordinates) for the radii of atomic highlights. If not provided
            the default value from \c drawOptions() will be used.
  \param confId: (optional) conformer ID to be used for atomic coordinates

*/
RDKIT_MOLDRAW2D_EXPORT void prepareAndDrawMolecule(
    MolDraw2D &drawer, const ROMol &mol, const std::string &legend = "",
    const std::vector<int> *highlight_atoms = nullptr,
    const std::vector<int> *highlight_bonds = nullptr,
    const std::map<int, DrawColour> *highlight_atom_map = nullptr,
    const std::map<int, DrawColour> *highlight_bond_map = nullptr,
    const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const char *json);
RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const std::string &json);
}  // namespace MolDraw2DUtils
}  // namespace RDKit
#endif  // MOLDRAW2DUTILS_H
