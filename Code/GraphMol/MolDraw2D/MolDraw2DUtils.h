//
//  Copyright (C) 2016-2021 Greg Landrum and other RDKit contributors
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

#include <tuple>

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
\param wavyBonds: calls addWavyBondsForStereoAny() and clears other markers that
   double bond stereo is unknown

NOTE: the kekulization step can fail, throwing a MolSanitizeExecption. If this
happens the molecule will be in an inconsistent, partially kekulized, state.
This isn't normally a problem for molecules that have been sanitized, but can be
problematic if the molecules have been modified post santitization.
*/
RDKIT_MOLDRAW2D_EXPORT void prepareMolForDrawing(
    RWMol &mol, bool kekulize = true, bool addChiralHs = true,
    bool wedgeBonds = true, bool forceCoords = false, bool wavyBonds = false);

//! prepare a molecule for drawing and draw it
/*
  \param mol: the molecule to draw
  \param legend: (optional) the legend (to be drawn under the molecule)
  \param highlight_atoms: (optional) vector of atom ids to highlight
  \param highlight_atoms: (optional) vector of bond ids to highlight
  \param highlight_atom_map: (optional) map from atomId -> DrawColour
            providing the highlight colors. If not provided the default
            highlight colour from \c drawOptions() will be used.
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
    const std::map<int, double> *highlight_radii = nullptr, int confId = -1,
    bool kekulize = true, bool addChiralHs = true, bool wedgeBonds = true,
    bool forceCoords = false, bool wavyBonds = false);

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const char *json);
RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const std::string &json);
RDKIT_MOLDRAW2D_EXPORT void updateMolDrawOptionsFromJSON(MolDrawOptions &opts,
                                                         const char *json);
RDKIT_MOLDRAW2D_EXPORT void updateMolDrawOptionsFromJSON(
    MolDrawOptions &opts, const std::string &json);

struct ContourParams {
  bool setScale = true;           // assumes the grid is drawn first
  bool dashNegative = true;       // use dashed lines for negative contours
  bool fillGrid = false;          // shade the grid
  double gridResolution = 0.15;   // spacing between elements of the grid
  double contourWidth = 1.0;      // linewidth for drawing contours
  double extraGridPadding = 0.0;  // extra padding (in molecule coordinates)
  DrawColour contourColour = {0.5, 0.5, 0.5,
                              0.5};  // color for drawing contours
  std::vector<DrawColour> colourMap = {
      {0.557, 0.004, 0.322, 0.5},
      {1, 1, 1, 0.5},
      {0.153, 0.392, 0.098, 0.5}};  // similarity map color scheme
  bool drawAsLines =
      true;  // draws the contours as continuous lines instead of line segments.
  double coordScaleForQuantization =
      1000.;  // caling factor used to convert coordinates to ints when forming
              // the continuous lines
  double isovalScaleForQuantization =
      1e6;  // scaling factor used to convert isovalues to ints when forming the
            // continuous lines
  bool useFillThreshold =
      false;  // use a magnitude threshold to determine if a grid box is filled
  double fillThreshold = 0.01;  // magnitude threshold for filling grid boxes
  bool fillThresholdIsFraction = true;  // if true, the fill threshold is a
                                        // fraction of the range of the data
};

//! Generates and draws contours for data on a grid
/*
  \param drawer: the MolDraw2D object to use
  \param grid: the data to be contoured
  \param xcoords: x positions of the grid points
  \param ycoords: y positions of the grid points
  \param nContours: the number of contours to draw
  \param levels: the contours to use
  \param ps: additional parameters controlling the contouring.
  \param mol: molecule to be used to adjust the scale of the drawing.
  If the \c levels argument is empty, the contour levels will be determined
  automatically from the max and min values on the grid and \c levels will
  be updated to include the contour levels.

  If \c ps.fillGrid is set, the data on the grid will also be drawn using
  the color scheme in \c ps.colourMap

  if the \c mol argument is given, it will be used to adjust the scale of
  drawing.  This is because a common use is to draw the molecule onto
  the contour, and it makes sense if it fits.

*/
RDKIT_MOLDRAW2D_EXPORT void contourAndDrawGrid(
    MolDraw2D &drawer, const double *grid, const std::vector<double> &xcoords,
    const std::vector<double> &ycoords, size_t nContours,
    std::vector<double> &levels, const ContourParams &ps = ContourParams(),
    const ROMol *mol = nullptr);
//! \overload
RDKIT_MOLDRAW2D_EXPORT inline void contourAndDrawGrid(
    MolDraw2D &drawer, const double *grid, const std::vector<double> &xcoords,
    const std::vector<double> &ycoords, size_t nContours = 10,
    const ContourParams &ps = ContourParams(), const ROMol *mol = nullptr) {
  std::vector<double> levels;
  contourAndDrawGrid(drawer, grid, xcoords, ycoords, nContours, levels, ps,
                     mol);
};

//! Generates and draws contours for a set of gaussians
/*
  \param drawer: the MolDraw2D object to use
  \param locs: locations of the gaussians
  \param heights: the heights (or weights) of the gaussians
  \param widths: the standard deviations of the gaussians
  \param nContours: the number of contours to draw
  \param levels: the contours to use
  \param ps: additional parameters controlling the contouring.
  \param mol: molecule to be used to adjust the scale of the drawing.

  The values are calculated on a grid with spacing \c ps.gridResolution.
  If \c ps.setScale  is set, the grid size will be calculated based on the
  locations of the gaussians and \c ps.extraGridPadding. Otherwise the current
  size of the viewport will be used.

  If the \c levels argument is empty, the contour levels will be determined
  automatically from the max and min values on the grid and \c levels will
  be updated to include the contour levels.

  If \c ps.fillGrid is set, the data on the grid will also be drawn using
  the color scheme in \c ps.colourMap

   if the \c mol argument is given, it will be used to adjust the scale of
  drawing.  This is because a common use is to draw the molecule onto
  the contour, and it makes sense if it fits.

*/
RDKIT_MOLDRAW2D_EXPORT void contourAndDrawGaussians(
    MolDraw2D &drawer, const std::vector<Point2D> &locs,
    const std::vector<double> &heights, const std::vector<double> &widths,
    size_t nContours, std::vector<double> &levels,
    const ContourParams &ps = ContourParams(), const ROMol *mol = nullptr);
//! \overload
RDKIT_MOLDRAW2D_EXPORT inline void contourAndDrawGaussians(
    MolDraw2D &drawer, const std::vector<Point2D> &locs,
    const std::vector<double> &heights, const std::vector<double> &widths,
    size_t nContours = 10, const ContourParams &ps = ContourParams(),
    const ROMol *mol = nullptr) {
  std::vector<double> levels;
  contourAndDrawGaussians(drawer, locs, heights, widths, nContours, levels, ps,
                          mol);
};

//! Draw a molecule to a MolDraw2D object according to ACS 1996 guidelines
/*
  The ACS1996 guidelines, as described at
  https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing
  A number of values in drawer.drawOptions() are changed.
  This is designed to be used with a flexiCanvas, i.e. a MolDraw2D object
  created with width and height -1, because it works to a fixed scale.
  It will issue a warning if the dimensions are otherwise and the picture may
  look sub-optimal.
 */
RDKIT_MOLDRAW2D_EXPORT void drawMolACS1996(
    MolDraw2D &drawer, const ROMol &mol, const std::string &legend,
    const std::vector<int> *highlight_atoms,
    const std::vector<int> *highlight_bonds,
    const std::map<int, DrawColour> *highlight_atom_map = nullptr,
    const std::map<int, DrawColour> *highlight_bond_map = nullptr,
    const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

//! Set the draw options to produce something as close as possible to
//! the ACS 1996 guidelines as described at
//! https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing
/*
 \param MolDrawOptions opt - the options what will be changed
 \param float meanBondLength - mean bond length of the molecule

 Works best if the MolDraw2D object is created with width and height -1 (a
 flexiCanvas).
 The mean bond length may be calculated with MolDraw2DUtils::meanBondLength.
 It is used to calculate the offset for the lines in multiple bonds.

 Options changed are:
   bondLineWidth = 0.6
   scaleBondWidth = false
   scalingFactor = 14.4 / meanBondLen
   multipleBondOffset = 0.18
   highlightBondWidthMultiplier = 32
   setMonochromeMode - black and white
   fixedFontSize = 10
   additionalAtomLabelPadding = 0.066
   fontFile - if it isn't set already, then if RDBASE is set and the file
              exists, uses $RDBASE/Fonts/Data/FreeSans.ttf.  Otherwise uses
              BuiltinRobotoRegular.
 */
RDKIT_MOLDRAW2D_EXPORT void setACS1996Options(MolDrawOptions &opts,
                                              double meanBondLen = 1.0);
RDKIT_MOLDRAW2D_EXPORT double meanBondLength(const ROMol &mol, int confId = -1);
}  // namespace MolDraw2DUtils

}  // namespace RDKit
#endif  // MOLDRAW2DUTILS_H
