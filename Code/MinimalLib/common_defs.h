//
//  Copyright (C) 2024 Novartis and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Chirality.h>
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
  bool useMolBlockWedging = false;
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

class DrawerFromDetails {
 public:
  virtual ~DrawerFromDetails() {}
  std::string draw_mol(const ROMol &mol) {
    MolDrawingDetails molDrawingDetails;
    molDrawingDetails.width = d_width;
    molDrawingDetails.height = d_height;
    if (!d_details.empty()) {
      auto problems = process_mol_details(d_details, molDrawingDetails);
      if (!problems.empty()) {
        return problems;
      }
    }
    initDrawer(molDrawingDetails);
    const ROMol *molPtr = &mol;
    std::unique_ptr<ROMol> origWedgingMol;
    if (molDrawingDetails.useMolBlockWedging) {
      origWedgingMol.reset(new ROMol(mol));
      Chirality::reapplyMolBlockWedging(*origWedgingMol);
      molPtr = origWedgingMol.get();
      drawer().drawOptions().useMolBlockWedging = false;
    }
    drawer().setOffset(molDrawingDetails.offsetx, molDrawingDetails.offsety);

    MolDraw2DUtils::prepareAndDrawMolecule(
        drawer(), *molPtr, molDrawingDetails.legend, &molDrawingDetails.atomIds,
        &molDrawingDetails.bondIds,
        molDrawingDetails.atomMap.empty() ? nullptr
                                          : &molDrawingDetails.atomMap,
        molDrawingDetails.bondMap.empty() ? nullptr
                                          : &molDrawingDetails.bondMap,
        molDrawingDetails.radiiMap.empty() ? nullptr
                                           : &molDrawingDetails.radiiMap,
        -1, molDrawingDetails.kekulize, molDrawingDetails.addChiralHs,
        molDrawingDetails.wedgeBonds, molDrawingDetails.forceCoords,
        molDrawingDetails.wavyBonds);

    return finalizeDrawing();
  }
  std::string draw_rxn(const ChemicalReaction &rxn) {
    RxnDrawingDetails rxnDrawingDetails;
    rxnDrawingDetails.width = d_width;
    rxnDrawingDetails.height = d_height;
    if (!d_details.empty()) {
      auto problems = process_rxn_details(d_details, rxnDrawingDetails);
      if (!problems.empty()) {
        return problems;
      }
    }
    initDrawer(rxnDrawingDetails);
    if (!rxnDrawingDetails.kekulize) {
      drawer().drawOptions().prepareMolsBeforeDrawing = false;
    }
    drawer().drawReaction(
        rxn, rxnDrawingDetails.highlightByReactant,
        !rxnDrawingDetails.highlightByReactant ||
                rxnDrawingDetails.highlightColorsReactants.empty()
            ? nullptr
            : &rxnDrawingDetails.highlightColorsReactants);
    return finalizeDrawing();
  }

 protected:
  void init(int w = -1, int h = -1,
            const std::string &details = std::string()) {
    d_width = w;
    d_height = h;
    d_details = details;
  }
  void updateDrawerParamsFromJSON() {
    if (!d_details.empty()) {
      MolDraw2DUtils::updateDrawerParamsFromJSON(drawer(), d_details);
    }
  }

 private:
  virtual MolDraw2D &drawer() const = 0;
  virtual void initDrawer(const DrawingDetails &drawingDetails) = 0;
  virtual std::string finalizeDrawing() = 0;
  int d_width;
  int d_height;
  std::string d_details;
};

}  // namespace MinimalLib
}  // namespace RDKit
