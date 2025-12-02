#define USE_BETTER_ENUMS
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/JSONHelpers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

#define PT_OPT_GET(opt) opts.opt = pt.get(#opt, opts.opt)

namespace RDKit {
namespace MolDraw2DUtils {

void updateMolDrawOptionsFromJSON(MolDrawOptions &opts, const char *json) {
  PRECONDITION(json, "no parameter string");
  updateMolDrawOptionsFromJSON(opts, std::string(json));
}

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const char *json) {
  updateMolDrawOptionsFromJSON(drawer.drawOptions(), json);
}

void get_rgba(const boost::property_tree::ptree &node, DrawColour &colour) {
  boost::property_tree::ptree::const_iterator itm = node.begin();
  colour.r = itm->second.get_value<float>();
  ++itm;
  colour.g = itm->second.get_value<float>();
  ++itm;
  colour.b = itm->second.get_value<float>();
  ++itm;
  if (itm != node.end()) {
    colour.a = itm->second.get_value<float>();
    ++itm;
  }
}

void get_colour_option(const boost::property_tree::ptree &pt, const char *pnm,
                       DrawColour &colour) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt.find(pnm) == pt.not_found()) {
    return;
  }

  const auto &node = pt.get_child(pnm);
  get_rgba(node, colour);
}

void get_colour_palette_option(const boost::property_tree::ptree &pt,
                               const char *pnm, ColourPalette &palette) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  static const std::map<std::string, void (*)(ColourPalette &)>
      PRESET_PALETTE_MAP{
          {"default", assignDefaultPalette},
          {"avalon", assignAvalonPalette},
          {"cdk", assignCDKPalette},
          {"darkmode", assignDarkModePalette},
          {"bw", assignBWPalette},
      };
  auto atomColourPaletteIt = pt.find(pnm);
  if (atomColourPaletteIt == pt.not_found()) {
    return;
  }
  // Does the "atomColourPalette" key correspond to a terminal value?
  if (atomColourPaletteIt->second.empty()) {
    const auto &v = atomColourPaletteIt->second.data();
    if (v.empty()) {
      return;
    }
    auto paletteName = boost::lexical_cast<std::string>(v);
    boost::algorithm::to_lower(paletteName);
    auto preSetPaletteIt = PRESET_PALETTE_MAP.find(paletteName);
    if (preSetPaletteIt == PRESET_PALETTE_MAP.end()) {
      return;
    }
    // Populate the palette calling the corresponding function
    preSetPaletteIt->second(palette);
  } else {
    for (const auto &atomicNumNode : atomColourPaletteIt->second) {
      int atomicNum = boost::lexical_cast<int>(atomicNumNode.first);
      DrawColour colour;
      get_rgba(atomicNumNode.second, colour);
      palette[atomicNum] = colour;
    }
  }
}

void get_highlight_style_option(const boost::property_tree::ptree &pt,
                                const char *pnm,
                                MultiColourHighlightStyle &mchs) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt.find(pnm) == pt.not_found()) {
    return;
  }
  const auto &node = pt.get_child(pnm);
  auto styleStr = node.get_value<std::string>();
  if (styleStr == "Lasso") {
    mchs = MultiColourHighlightStyle::LASSO;
  } else if (styleStr == "CircleAndLine") {
    mchs = MultiColourHighlightStyle::CIRCLEANDLINE;
  }
}

void updateMolDrawOptionsFromJSON(MolDrawOptions &opts,
                                  const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(atomLabelDeuteriumTritium);
  PT_OPT_GET(dummiesAreAttachments);
  PT_OPT_GET(circleAtoms);
  PT_OPT_GET(splitBonds);
  PT_OPT_GET(continuousHighlight);
  PT_OPT_GET(fillHighlights);
  PT_OPT_GET(highlightRadius);
  PT_OPT_GET(flagCloseContactsDist);
  PT_OPT_GET(includeAtomTags);
  PT_OPT_GET(clearBackground);
  PT_OPT_GET(legendFontSize);
  PT_OPT_GET(legendFraction);
  PT_OPT_GET(maxFontSize);
  PT_OPT_GET(minFontSize);
  PT_OPT_GET(fixedFontSize);
  PT_OPT_GET(baseFontSize);
  PT_OPT_GET(annotationFontScale);
  PT_OPT_GET(fontFile);
  PT_OPT_GET(multipleBondOffset);
  PT_OPT_GET(padding);
  PT_OPT_GET(componentPadding);
  PT_OPT_GET(additionalAtomLabelPadding);
  PT_OPT_GET(noAtomLabels);
  PT_OPT_GET(bondLineWidth);
  PT_OPT_GET(scaleBondWidth);
  PT_OPT_GET(scaleHighlightBondWidth);
  PT_OPT_GET(highlightBondWidthMultiplier);
  PT_OPT_GET(prepareMolsBeforeDrawing);
  PT_OPT_GET(fixedScale);
  PT_OPT_GET(fixedBondLength);
  PT_OPT_GET(rotate);
  PT_OPT_GET(addAtomIndices);
  PT_OPT_GET(addBondIndices);
  PT_OPT_GET(isotopeLabels);
  PT_OPT_GET(dummyIsotopeLabels);
  PT_OPT_GET(addStereoAnnotation);
  PT_OPT_GET(showAllCIPCodes);
  PT_OPT_GET(atomHighlightsAreCircles);
  PT_OPT_GET(centreMoleculesBeforeDrawing);
  PT_OPT_GET(explicitMethyl);
  PT_OPT_GET(includeMetadata);
  PT_OPT_GET(includeRadicals);
  PT_OPT_GET(comicMode);
  PT_OPT_GET(variableBondWidthMultiplier);
  PT_OPT_GET(variableAtomRadius);
  PT_OPT_GET(includeChiralFlagLabel);
  PT_OPT_GET(simplifiedStereoGroupLabel);
  PT_OPT_GET(unspecifiedStereoIsUnknown);
  PT_OPT_GET(singleColourWedgeBonds);
  PT_OPT_GET(useMolBlockWedging);
  PT_OPT_GET(scalingFactor);
  PT_OPT_GET(drawMolsSameScale);
  PT_OPT_GET(useComplexQueryAtomSymbols);
  PT_OPT_GET(bracketsAroundAtomLists);
  PT_OPT_GET(standardColoursForHighlightedAtoms);

  get_colour_option(pt, "highlightColour", opts.highlightColour);
  get_colour_option(pt, "backgroundColour", opts.backgroundColour);
  get_colour_option(pt, "queryColour", opts.queryColour);
  get_colour_option(pt, "legendColour", opts.legendColour);
  get_colour_option(pt, "symbolColour", opts.symbolColour);
  get_colour_option(pt, "annotationColour", opts.annotationColour);
  get_colour_option(pt, "atomNoteColour", opts.atomNoteColour);
  get_colour_option(pt, "bondNoteColour", opts.bondNoteColour);
  get_colour_option(pt, "variableAttachmentColour",
                    opts.variableAttachmentColour);
  get_colour_palette_option(pt, "atomColourPalette", opts.atomColourPalette);
  if (pt.find("atomLabels") != pt.not_found()) {
    for (const auto &item : pt.get_child("atomLabels")) {
      opts.atomLabels[boost::lexical_cast<int>(item.first)] =
          item.second.get_value<std::string>();
    }
  }
  get_highlight_style_option(pt, "multiColourHighlightStyle",
                             opts.multiColourHighlightStyle);
  const auto drawingExtentsIncludeIt = pt.find("drawingExtentsInclude");
  if (drawingExtentsIncludeIt != pt.not_found()) {
    bool haveDrawElementFlags = false;
    auto drawingExtentsInclude = flagsFromJson<DrawElement>(drawingExtentsIncludeIt->second, &haveDrawElementFlags);
    if (haveDrawElementFlags) {
      opts.drawingExtentsInclude = drawingExtentsInclude;
    }
  }
}

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(
    MolDraw2D &drawer, const std::string &json) {
  updateMolDrawOptionsFromJSON(drawer.drawOptions(), json);
}

}  // namespace MolDraw2DUtils
}  // namespace RDKit
