//
//  Copyright (c) 2024, Glysade Inc
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "fragment.h"
#include "bond.h"
#include "node.h"

namespace RDKit {
namespace ChemDraw {
namespace {
const char *sequenceTypeToName(CDXSeqType seqtype) {
  switch (seqtype) {
    case kCDXSeqType_Unknown:
      return "Unknown";
    case kCDXSeqType_Peptide:
      return "Peptide (Helm)";  // HELM peptides
    case kCDXSeqType_Peptide1:
      return "Peptide1 (Single Letter Amino Acid)";  // Single letter amino
                                                     // acids (Legacy biopolymer
                                                     // support)
    case kCDXSeqType_Peptide3:
      return "Peptide3 (Three letter amino acid)";  // Three letter amino acids
                                                    // (Legacy biopolymer
                                                    // support)
    case kCDXSeqType_DNA:
      return "DNA";
    case kCDXSeqType_RNA:
      return "RNA";
    case kCDXSeqType_Biopolymer:
      return "Biopolymer";
    default:
      return "";
  }
}
}  // namespace
bool parse_fragment(RWMol &mol, CDXFragment &fragment, PageData &pagedata,
                    int &missing_frag_id, int external_attachment) {
  int frag_id = fragment.GetObjectID();
  if (fragment.m_sequenceType != kCDXSeqType_Unknown) {
    BOOST_LOG(rdWarningLog)
        << "Unhandled chemdraw sequence type "
        << sequenceTypeToName(fragment.m_sequenceType) << std::endl;
    return false;
  }
  if (frag_id == -1) {
    // ChemDraw simply assigns a new one
    BOOST_LOG(rdWarningLog)
        << "Invalid or missing fragment id from CDXML fragment, assigning new one..."
        << std::endl;
    frag_id = missing_frag_id;
    missing_frag_id--;
  }
  mol.setProp(CDX_FRAG_ID, frag_id);

  // for atom in frag
  std::map<std::pair<int, StereoGroupType>, StereoGroupInfo> sgroups;

  // nodetypes =
  // https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/properties/Node_Type.htm
  bool skip_fragment =
      false;  // is there an irrecoverable error for this fragment

  for (auto child : fragment.ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)child.second->GetTag();
#ifdef DEBUG
    std::cerr << "Data Type: " << id << std::endl;
#endif
    switch (id) {
      case kCDXObj_Node: {
        CDXNode &node = (CDXNode &)(*child.second);
        if (!parse_node(mol, frag_id, node, pagedata, sgroups, missing_frag_id,
                        external_attachment)) {
          skip_fragment = true;
        }
        break;
      }
      case kCDXObj_Bond: {
        CDXBond &bond = (CDXBond &)(*child.second);
        if (!parse_bond(mol, frag_id, bond, pagedata)) {
          skip_fragment = true;
          break;
        }
      }
      case kCDXProp_EndObject: break;
      case kCDXProp_CreationUserName: break;
      case kCDXProp_CreationDate: break;
      case kCDXProp_CreationProgram: break;
      case kCDXProp_ModificationUserName: break;
      case kCDXProp_ModificationDate: break;
      case kCDXProp_ModificationProgram: break;
      case kCDXProp_Unused1: break;
      case kCDXProp_Name: break;
      case kCDXProp_Comment: break;
      case kCDXProp_ZOrder: break;
      case kCDXProp_RegistryNumber: break;
      case kCDXProp_RegistryAuthority: break;
      case kCDXProp_Unused2: break;
      case kCDXProp_RepresentsProperty: break;
      case kCDXProp_IgnoreWarnings: break;
      case kCDXProp_ChemicalWarning: break;
      case kCDXProp_Visible: break;
      case kCDXProp_Transparent: break;
      case kCDXProp_SupersededBy: break;
      case kCDXProp_StructurePerspective: break;
      case kCDXProp_FontTable: break;
      case kCDXProp_2DPosition: break;
      case kCDXProp_3DPosition: break;
      case kCDXProp_2DExtent: break;
      case kCDXProp_3DExtent: break;
      case kCDXProp_BoundingBox: break;
      case kCDXProp_RotationAngle: break;
      case kCDXProp_BoundsInParent: break;
      case kCDXProp_3DHead: break;
      case kCDXProp_3DTail: break;
      case kCDXProp_TopLeft: break;
      case kCDXProp_TopRight: break;
      case kCDXProp_BottomRight: break;
      case kCDXProp_BottomLeft: break;
      case kCDXProp_3DCenter: break;
      case kCDXProp_3DMajorAxisEnd: break;
      case kCDXProp_3DMinorAxisEnd: break;
      case kCDXProp_ColorTable: break;
      case kCDXProp_ForegroundColor: break;
      case kCDXProp_BackgroundColor: break;
      case kCDXProp_FadePercent: break;
      case kCDXProp_Unused8: break;
      case kCDXProp_Unused9: break;
      case kCDXProp_ForegroundAlpha: break;
      case kCDXProp_BackgroundAlpha: break;
      case kCDXProp_HighlightColor: break;
      case kCDXProp_Node_Type: break;
      case kCDXProp_Node_LabelDisplay: break;
      case kCDXProp_Node_Element: break;
      case kCDXProp_Atom_ElementList: break;
      case kCDXProp_Atom_Formula: break;
      case kCDXProp_Atom_Isotope: break;
      case kCDXProp_Atom_Charge: break;
      case kCDXProp_Atom_Radical: break;
      case kCDXProp_Atom_RestrictFreeSites: break;
      case kCDXProp_Atom_RestrictImplicitHydrogens: break;
      case kCDXProp_Atom_RestrictRingBondCount: break;
      case kCDXProp_Atom_RestrictUnsaturatedBonds: break;
      case kCDXProp_Atom_RestrictRxnChange: break;
      case kCDXProp_Atom_RestrictRxnStereo: break;
      case kCDXProp_Atom_AbnormalValence: break;
      case kCDXProp_Unused3: break;
      case kCDXProp_Atom_NumHydrogens: break;
      case kCDXProp_Unused4: break;
      case kCDXProp_Unused5: break;
      case kCDXProp_Atom_HDot: break;
      case kCDXProp_Atom_HDash: break;
      case kCDXProp_Atom_Geometry: break;
      case kCDXProp_Atom_BondOrdering: break;
      case kCDXProp_Node_Attachments: break;
      case kCDXProp_Atom_GenericNickname: break;
      case kCDXProp_Atom_AltGroupID: break;
      case kCDXProp_Atom_RestrictSubstituentsUpTo: break;
      case kCDXProp_Atom_RestrictSubstituentsExactly: break;
      case kCDXProp_Atom_CIPStereochemistry: break;
      case kCDXProp_Atom_Translation: break;
      case kCDXProp_Atom_AtomNumber: break;
      case kCDXProp_Atom_ShowQuery: break;
      case kCDXProp_Atom_ShowStereo: break;
      case kCDXProp_Atom_ShowAtomNumber: break;
      case kCDXProp_Atom_LinkCountLow: break;
      case kCDXProp_Atom_LinkCountHigh: break;
      case kCDXProp_Atom_IsotopicAbundance: break;
      case kCDXProp_Atom_ExternalConnectionType: break;
      case kCDXProp_Atom_GenericList: break;
      case kCDXProp_Atom_ShowTerminalCarbonLabels: break;
      case kCDXProp_Atom_ShowNonTerminalCarbonLabels: break;
      case kCDXProp_Atom_HideImplicitHydrogens: break;
      case kCDXProp_Atom_ShowEnhancedStereo: break;
      case kCDXProp_Atom_EnhancedStereoType: break;
      case kCDXProp_Atom_EnhancedStereoGroupNum: break;
      case kCDXProp_Node_NeedsClean: break;
      case kCDXProp_Atom_ResidueID: break;
      case kCDXProp_Atom_ShowResidueID: break;
      case kCDXProp_Atom_ExternalConnectionNum: break;
      case kCDXProp_Atom_ShowAtomID: break;
      case kCDXProp_Atom_AtomID: break;
      case kCDXProp_Node_HydrogenBondAttachmentAtoms: break;
      case kCDXProp_Node_HydrogenBonds: break;
      case kCDXProp_Mole_Racemic: break;
      case kCDXProp_Mole_Absolute: break;
      case kCDXProp_Mole_Relative: break;
      case kCDXProp_Mole_Formula: break;
      case kCDXProp_Mole_Weight: break;
      case kCDXProp_Frag_ConnectionOrder: break;
      case kCDXProp_Frag_SequenceType: break;
      case kCDXProp_Frag_IsFromGuidedStereo: break;
      case kCDXProp_Frag_IsComplement: break;
      case kCDXProp_Bond_Order: break;
      case kCDXProp_Bond_Display: break;
      case kCDXProp_Bond_Display2: break;
      case kCDXProp_Bond_DoublePosition: break;
      case kCDXProp_Bond_Begin: break;
      case kCDXProp_Bond_End: break;
      case kCDXProp_Bond_RestrictTopology: break;
      case kCDXProp_Bond_RestrictRxnParticipation: break;
      case kCDXProp_Bond_BeginAttach: break;
      case kCDXProp_Bond_EndAttach: break;
      case kCDXProp_Bond_CIPStereochemistry: break;
      case kCDXProp_Bond_BondOrdering: break;
      case kCDXProp_Bond_ShowQuery: break;
      case kCDXProp_Bond_ShowStereo: break;
      case kCDXProp_Bond_CrossingBonds: break;
      case kCDXProp_Bond_ShowRxn: break;
      case kCDXProp_Bond_Connectivity: break;
      case kCDXProp_Bond_BeginExternalNum: break;
      case kCDXProp_Bond_EndExternalNum: break;
      case kCDXProp_Bond_Connectivity_Routed: break;
      case kCDXProp_Text: break;
      case kCDXProp_Justification: break;
      case kCDXProp_LineHeight: break;
      case kCDXProp_WordWrapWidth: break;
      case kCDXProp_LineStarts: break;
      case kCDXProp_LabelAlignment: break;
      case kCDXProp_LabelLineHeight: break;
      case kCDXProp_CaptionLineHeight: break;
      case kCDXProp_InterpretChemically: break;
      case kCDXProp_UTF8Text: break;
      case kCDXProp_MacPrintInfo: break;
      case kCDXProp_WinPrintInfo: break;
      case kCDXProp_PrintMargins: break;
      case kCDXProp_ChainAngle: break;
      case kCDXProp_BondSpacing: break;
      case kCDXProp_BondLength: break;
      case kCDXProp_BoldWidth: break;
      case kCDXProp_LineWidth: break;
      case kCDXProp_MarginWidth: break;
      case kCDXProp_HashSpacing: break;
      case kCDXProp_LabelStyle: break;
      case kCDXProp_CaptionStyle: break;
      case kCDXProp_CaptionJustification: break;
      case kCDXProp_FractionalWidths: break;
      case kCDXProp_Magnification: break;
      case kCDXProp_WidthPages: break;
      case kCDXProp_HeightPages: break;
      case kCDXProp_DrawingSpaceType: break;
      case kCDXProp_Width: break;
      case kCDXProp_Height: break;
      case kCDXProp_PageOverlap: break;
      case kCDXProp_Header: break;
      case kCDXProp_HeaderPosition: break;
      case kCDXProp_Footer: break;
      case kCDXProp_FooterPosition: break;
      case kCDXProp_PrintTrimMarks: break;
      case kCDXProp_LabelStyleFont: break;
      case kCDXProp_CaptionStyleFont: break;
      case kCDXProp_LabelStyleSize: break;
      case kCDXProp_CaptionStyleSize: break;
      case kCDXProp_LabelStyleFace: break;
      case kCDXProp_CaptionStyleFace: break;
      case kCDXProp_LabelStyleColor: break;
      case kCDXProp_CaptionStyleColor: break;
      case kCDXProp_BondSpacingAbs: break;
      case kCDXProp_LabelJustification: break;
      case kCDXProp_FixInplaceExtent: break;
      case kCDXProp_Side: break;
      case kCDXProp_FixInplaceGap: break;
      case kCDXProp_CartridgeData: break;
      case kCDXProp_AminoAcidTermini: break;
      case kCDXProp_ShowSequenceTermini: break;
      case kCDXProp_ShowSequenceBonds: break;
      case kCDXProp_ResidueWrapCount: break;
      case kCDXProp_ResidueBlockCount: break;
      case kCDXProp_Unused10: break;
      case kCDXProp_Unused11: break;
      case kCDXProp_BondSpacingType: break;
      case kCDXProp_LabelStyleFontName: break;
      case kCDXProp_CaptionStyleFontName: break;
      case kCDXProp_ShowSequenceUnlinkedBranches: break;
      case kCDXProp_MonomerRenderingStyle: break;
      case kCDXProp_Window_IsZoomed: break;
      case kCDXProp_Window_Position: break;
      case kCDXProp_Window_Size: break;
      case kCDXProp_Graphic_Type: break;
      case kCDXProp_Line_Type: break;
      case kCDXProp_Arrow_Type: break;
      case kCDXProp_Rectangle_Type: break;
      case kCDXProp_Oval_Type: break;
      case kCDXProp_Orbital_Type: break;
      case kCDXProp_Bracket_Type: break;
      case kCDXProp_Symbol_Type: break;
      case kCDXProp_Curve_Type: break;
      case kCDXProp_Arrowhead_Size: break;
      case kCDXProp_Arc_AngularSize: break;
      case kCDXProp_Bracket_LipSize: break;
      case kCDXProp_Curve_Points: break;
      case kCDXProp_Bracket_Usage: break;
      case kCDXProp_Polymer_RepeatPattern: break;
      case kCDXProp_Polymer_FlipType: break;
      case kCDXProp_BracketedObjects: break;
      case kCDXProp_Bracket_RepeatCount: break;
      case kCDXProp_Bracket_ComponentOrder: break;
      case kCDXProp_Bracket_SRULabel: break;
      case kCDXProp_Bracket_GraphicID: break;
      case kCDXProp_Bracket_BondID: break;
      case kCDXProp_Bracket_InnerAtomID: break;
      case kCDXProp_Curve_Points3D: break;
      case kCDXProp_Arrowhead_Type: break;
      case kCDXProp_Arrowhead_CenterSize: break;
      case kCDXProp_Arrowhead_Width: break;
      case kCDXProp_ShadowSize: break;
      case kCDXProp_Arrow_ShaftSpacing: break;
      case kCDXProp_Arrow_EquilibriumRatio: break;
      case kCDXProp_Arrowhead_Head: break;
      case kCDXProp_Arrowhead_Tail: break;
      case kCDXProp_Fill_Type: break;
      case kCDXProp_Curve_Spacing: break;
      case kCDXProp_Closed: break;
      case kCDXProp_Arrow_Dipole: break;
      case kCDXProp_Arrow_NoGo: break;
      case kCDXProp_CornerRadius: break;
      case kCDXProp_Frame_Type: break;
      case kCDXProp_Arrow_SourceID: break;
      case kCDXProp_Arrow_TargetID: break;
      case kCDXProp_Arrow_IsSmart_Deleted: break;
      case kCDXProp_Picture_Edition: break;
      case kCDXProp_Picture_EditionAlias: break;
      case kCDXProp_MacPICT: break;
      case kCDXProp_WindowsMetafile: break;
      case kCDXProp_OLEObject: break;
      case kCDXProp_EnhancedMetafile: break;
      case kCDXProp_Compressed_MacPICT: break;
      case kCDXProp_Compressed_WindowsMetafile: break;
      case kCDXProp_Compressed_OLEObject: break;
      case kCDXProp_Compressed_EnhancedMetafile: break;
      case kCDXProp_Uncompressed_MacPICT_Size: break;
      case kCDXProp_Uncompressed_WindowsMetafile_Size: break;
      case kCDXProp_Uncompressed_OLEObject_Size: break;
      case kCDXProp_Uncompressed_EnhancedMetafile_Size: break;
      case kCDXProp_GIF: break;
      case kCDXProp_TIFF: break;
      case kCDXProp_PNG: break;
      case kCDXProp_JPEG: break;
      case kCDXProp_BMP: break;
      case kCDXProp_PDF: break;
      case kCDXProp_Spectrum_XSpacing: break;
      case kCDXProp_Spectrum_XLow: break;
      case kCDXProp_Spectrum_XType: break;
      case kCDXProp_Spectrum_YType: break;
      case kCDXProp_Spectrum_XAxisLabel: break;
      case kCDXProp_Spectrum_YAxisLabel: break;
      case kCDXProp_Spectrum_DataPoint: break;
      case kCDXProp_Spectrum_Class: break;
      case kCDXProp_Spectrum_YLow: break;
      case kCDXProp_Spectrum_YScale: break;
      case kCDXProp_TLC_OriginFraction: break;
      case kCDXProp_TLC_SolventFrontFraction: break;
      case kCDXProp_TLC_ShowOrigin: break;
      case kCDXProp_TLC_ShowSolventFront: break;
      case kCDXProp_ShowBorders: break;
      case kCDXProp_TLC_ShowSideTicks: break;
      case kCDXProp_TLC_Rf: break;
      case kCDXProp_TLC_Tail: break;
      case kCDXProp_TLC_ShowRf: break;
      case kCDXProp_GEP_ShowScale: break;
      case kCDXProp_GEP_ScaleUnit: break;
      case kCDXProp_GEP_StartRange: break;
      case kCDXProp_GEP_EndRange: break;
      case kCDXProp_GEP_ShowValue: break;
      case kCDXProp_GEP_Value: break;
      case kCDXProp_GEP_LaneLabelsAngle: break;
      case kCDXProp_GEP_AxisWidth: break;
      case kCDXProp_BioShape_Type: break;
      case kCDXProp_1SubstrateEnzyme_ReceptorSize: break;
      case kCDXProp_Receptor_NeckWidth: break;
      case kCDXProp_HelixProtein_CylinderWidth: break;
      case kCDXProp_HelixProtein_CylinderHeight: break;
      case kCDXProp_HelixProtein_CylinderDistance: break;
      case kCDXProp_HelixProtein_PipeWidth: break;
      case kCDXProp_HelixProtein_Extra: break;
      case kCDXProp_Membrane_ElementSize: break;
      case kCDXProp_Membrane_StartAngle: break;
      case kCDXProp_Membrane_EndAngle: break;
      case kCDXProp_DNA_WaveLength: break;
      case kCDXProp_DNA_WaveWidth: break;
      case kCDXProp_DNA_Offset: break;
      case kCDXProp_DNA_WaveHeight: break;
      case kCDXProp_Gprotein_UpperHeight: break;
      case kCDXProp_NamedAlternativeGroup_TextFrame: break;
      case kCDXProp_NamedAlternativeGroup_GroupFrame: break;
      case kCDXProp_NamedAlternativeGroup_Valence: break;
      case kCDXProp_GeometricFeature: break;
      case kCDXProp_RelationValue: break;
      case kCDXProp_BasisObjects: break;
      case kCDXProp_ConstraintType: break;
      case kCDXProp_ConstraintMin: break;
      case kCDXProp_ConstraintMax: break;
      case kCDXProp_IgnoreUnconnectedAtoms: break;
      case kCDXProp_DihedralIsChiral: break;
      case kCDXProp_PointIsDirected: break;
      case kCDXProp_ChemicalPropertyType: break;
      case kCDXProp_ChemicalPropertyDisplayID: break;
      case kCDXProp_ChemicalPropertyIsActive: break;
      case kCDXProp_ChemicalPropertyUnknown: break;
      case kCDXProp_ChemicalPropertyName: break;
      case kCDXProp_ChemicalPropertyFormula: break;
      case kCDXProp_ChemicalPropertyExactMass: break;
      case kCDXProp_ChemicalPropertyMolWeight: break;
      case kCDXProp_ChemicalPropertyMOverZ: break;
      case kCDXProp_ChemicalPropertyAnalysis: break;
      case kCDXProp_ChemicalPropertyBoilingPoint: break;
      case kCDXProp_ChemicalPropertyMeltingPoint: break;
      case kCDXProp_ChemicalPropertyCriticalTemp: break;
      case kCDXProp_ChemicalPropertyCriticalPressure: break;
      case kCDXProp_ChemicalPropertyCriticalVolume: break;
      case kCDXProp_ChemicalPropertyGibbsEnergy: break;
      case kCDXProp_ChemicalPropertyLogP: break;
      case kCDXProp_ChemicalPropertyMR: break;
      case kCDXProp_ChemicalPropertyHenrysLaw: break;
      case kCDXProp_ChemicalPropertyHeatOfForm: break;
      case kCDXProp_ChemicalPropertytPSA: break;
      case kCDXProp_ChemicalPropertyCLogP: break;
      case kCDXProp_ChemicalPropertyCMR: break;
      case kCDXProp_ChemicalPropertyLogS: break;
      case kCDXProp_ChemicalPropertyPKa: break;
      case kCDXProp_ChemicalPropertyID: break;
      case kCDXProp_ChemicalPropertyFragmentLabel: break;
      case kCDXProp_ChemicalPropertyTypeIUPACAtomNumber: break;
      case kCDXProp_ChemicalPropertyIsChemicallySignificant: break;
      case kCDXProp_ChemicalPropertyExternalBonds: break;
      case kCDXProp_ReactionStep_Atom_Map: break;
      case kCDXProp_ReactionStep_Reactants: break;
      case kCDXProp_ReactionStep_Products: break;
      case kCDXProp_ReactionStep_Plusses: break;
      case kCDXProp_ReactionStep_Arrows: break;
      case kCDXProp_ReactionStep_ObjectsAboveArrow: break;
      case kCDXProp_ReactionStep_ObjectsBelowArrow: break;
      case kCDXProp_ReactionStep_Atom_Map_Manual: break;
      case kCDXProp_ReactionStep_Atom_Map_Auto: break;
      case kCDXProp_RxnAutonumber_Style: break;
      case kCDXProp_RxnAutonumber_Conditions: break;
      case kCDXProp_RxnAutonumber_Start: break;
      case kCDXProp_RxnAutonumber_Format: break;
      case kCDXProp_ObjectTag_Type: break;
      case kCDXProp_Unused6: break;
      case kCDXProp_Unused7: break;
      case kCDXProp_ObjectTag_Tracking: break;
      case kCDXProp_ObjectTag_Persistent: break;
      case kCDXProp_ObjectTag_Value: break;
      case kCDXProp_Positioning: break;
      case kCDXProp_PositioningAngle: break;
      case kCDXProp_PositioningOffset: break;
      case kCDXProp_Sequence_Identifier: break;
      case kCDXProp_CrossReference_Container: break;
      case kCDXProp_CrossReference_Document: break;
      case kCDXProp_CrossReference_Identifier: break;
      case kCDXProp_CrossReference_Sequence: break;
      case kCDXProp_Template_PaneHeight: break;
      case kCDXProp_Template_NumRows: break;
      case kCDXProp_Template_NumColumns: break;
      case kCDXProp_Group_Integral: break;
      case kCDXProp_SG_DataType: break;
      case kCDXProp_SG_PropertyType: break;
      case kCDXProp_SG_DataValue: break;
      case kCDXProp_SG_ComponentIsReactant: break;
      case kCDXProp_SG_ComponentIsHeader: break;
      case kCDXProp_IsHidden: break;
      case kCDXProp_IsReadOnly: break;
      case kCDXProp_IsEdited: break;
      case kCDXProp_SG_ComponentReferenceID: break;
      case kCDXProp_PlasmidMap_NumberBasePairs: break;
      case kCDXProp_PlasmidMap_MarkerStart: break;
      case kCDXProp_PlasmidMap_MarkerOffset: break;
      case kCDXProp_PlasmidMap_MarkerAngle: break;
      case kCDXProp_PlasmidMap_RegionStart: break;
      case kCDXProp_PlasmidMap_RegionEnd: break;
      case kCDXProp_PlasmidMap_RegionOffset: break;
      case kCDXProp_PlasmidMap_RingRadius: break;
      case kCDXProp_RLogic_Group: break;
      case kCDXProp_RLogic_Occurrence: break;
      case kCDXProp_RLogic_RestH: break;
      case kCDXProp_RLogic_IfThenGroup: break;
      case kCDXProp_Annotation_Keyword: break;
      case kCDXProp_Annotation_Content: break;
      case kCDXProp_SplitterPositions: break;
      case kCDXProp_PageDefinition: break;
      case kCDXProp_Property_Rule: break;
      case kCDXProp_Property_DataType: break;
      case kCDXProp_Property_Value: break;
      case kCDXUser_TemporaryBegin: break;
      case kCDXUser_TemporaryEnd: break;
      case kCDXObj_Document: break;
      case kCDXObj_Page: break;
      case kCDXObj_Group: break;
      case kCDXObj_Fragment: break;
      case kCDXObj_Text: break;
      case kCDXObj_Graphic: break;
      case kCDXObj_Curve: break;
      case kCDXObj_EmbeddedObject: break;
      case kCDXObj_NamedAlternativeGroup: break;
      case kCDXObj_TemplateGrid: break;
      case kCDXObj_RegistryNumber: break;
      case kCDXObj_ReactionScheme: break;
      case kCDXObj_ReactionStep: break;
      case kCDXObj_ObjectDefinition: break;
      case kCDXObj_Spectrum: break;
      case kCDXObj_ObjectTag: break;
      case kCDXObj_OleClientItem: break;
      case kCDXObj_Sequence: break;
      case kCDXObj_CrossReference: break;
      case kCDXObj_Splitter: break;
      case kCDXObj_Table: break;
      case kCDXObj_BracketedGroup: break;
      case kCDXObj_BracketAttachment: break;
      case kCDXObj_CrossingBond: break;
      case kCDXObj_Border: break;
      case kCDXObj_Geometry: break;
      case kCDXObj_Constraint: break;
      case kCDXObj_TLCPlate: break;
      case kCDXObj_TLCLane: break;
      case kCDXObj_TLCSpot: break;
      case kCDXObj_ChemicalProperty: break;
      case kCDXObj_Arrow: break;
      case kCDXObj_StoichiometryGrid: break;
      case kCDXObj_SGComponent: break;
      case kCDXObj_SGDatum: break;
      case kCDXObj_BioShape: break;
      case kCDXObj_PlasmidMap: break;
      case kCDXObj_PlasmidMarker: break;
      case kCDXObj_PlasmidRegion: break;
      case kCDXObj_RLogic: break;
      case kCDXObj_RLogicItem: break;
      case kCDXObj_Annotation: break;
      case kCDXObj_GEPPlate: break;
      case kCDXObj_GEPBand: break;
      case kCDXObj_Marker: break;
      case kCDXObj_GEPLane: break;
      case kCDXObj_DocumentProperties: break;
      case kCDXObj_Property: break;
      case kCDXObj_ColoredMolecularArea: break;
      case kCDXObj_UnknownObject: break;       
    }
  }

  // Add the stereo groups
  if (!sgroups.empty()) {
    std::vector<StereoGroup> stereo_groups;
    for (auto &sgroup : sgroups) {
      unsigned gId = 0;
      if (sgroup.second.grouptype != StereoGroupType::STEREO_ABSOLUTE &&
          sgroup.second.sgroup > 0) {
        gId = sgroup.second.sgroup;
      }
      std::vector<Bond *> newBonds;
      stereo_groups.emplace_back(sgroup.second.grouptype, sgroup.second.atoms,
                                 newBonds, gId);
    }
    mol.setStereoGroups(std::move(stereo_groups));
  }

  return !skip_fragment;
}
}
}  // namespace RDKit
