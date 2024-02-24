//
//  Copyright (C) 2024 Novartis and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <string>
#include <vector>

namespace RDKit {
namespace MinimalLib {

struct DrawingDetails {
  int width = -1;
  int height = -1;
  int offsetx = 0;
  int offsety = 0;
  int panelWidth = -1;
  int panelHeight = -1;
  bool noFreetype = false;
  bool kekulize = true;
  bool addChiralHs = true;
  bool wedgeBonds = true;
  bool forceCoords = false;
  bool wavyBonds = false;
  std::string legend;
  std::vector<int> atomIds;
  std::vector<int> bondIds;
};

struct MolDrawingDetails : public DrawingDetails {
  std::map<int, DrawColour> atomMap;
  std::map<int, DrawColour> bondMap;
  std::map<int, double> radiiMap;
};

struct RxnDrawingDetails : public DrawingDetails {
  bool highlightByReactant = false;
  std::vector<DrawColour> highlightColorsReactants;
};

extern std::string process_mol_details(const std::string &details,
                                       MolDrawingDetails &molDrawingDetails);
[[deprecated(
    "please use the overload taking MolDrawingDetails& as parameter")]] extern std::
    string
    process_mol_details(const std::string &details, int &width, int &height,
                        int &offsetx, int &offsety, std::string &legend,
                        std::vector<int> &atomIds, std::vector<int> &bondIds,
                        std::map<int, DrawColour> &atomMap,
                        std::map<int, DrawColour> &bondMap,
                        std::map<int, double> &radiiMap, bool &kekulize,
                        bool &addChiralHs, bool &wedgeBonds, bool &forceCoords,
                        bool &wavyBonds);

extern std::string process_rxn_details(const std::string &details,
                                       RxnDrawingDetails &rxnDrawingDetails);
[[deprecated(
    "please use the overload taking RxnDrawingDetails& as parameter")]] extern std::
    string
    process_rxn_details(const std::string &details, int &width, int &height,
                        int &offsetx, int &offsety, std::string &legend,
                        std::vector<int> &atomIds, std::vector<int> &bondIds,
                        bool &kekulize, bool &highlightByReactant,
                        std::vector<DrawColour> &highlightColorsReactants);

}  // namespace MinimalLib
}  // namespace RDKit
