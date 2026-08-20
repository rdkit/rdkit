//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//
#include "MacroMolTemplateLibrary.h"

#include <array>
#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {
namespace {

// Generated from the validated built-in template graphs; no parser is used at runtime.
struct DirectAtomDef { unsigned atomicNum; int formalCharge; unsigned isotope; Atom::ChiralType chiralTag; unsigned explicitHs; bool noImplicit; const char *pdbName; };
struct DirectBondDef { unsigned beginAtomIdx; unsigned endAtomIdx; Bond::BondType type; bool aromatic; Bond::BondStereo stereo; std::vector<unsigned int> stereoAtoms; };
struct BuiltinTemplateDef { const char *symbol; const char *name; MonomerClass monomerClass; const char *originalData; std::vector<DirectAtomDef> atoms; std::vector<DirectBondDef> bonds; std::vector<unsigned int> mainAtomIdxs; std::vector<MacroMolLeavingGroup> leavingGroups; };


const std::array<BuiltinTemplateDef, 22> builtinTemplates{{
    {"A", "ALA", MonomerClass::AminoAcid, "C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CB :1.pdbName. CA :2.pdbName. N :3.pdbName. H :4.pdbName. C :5.pdbName. O :6.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 4, 5}, {
      {{3}, 2, 3, 1},
      {{6}, 4, 6, 2},
     }},
    {"R", "ARG", MonomerClass::AminoAcid, "N=C(N)NCCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD :5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N :9.pdbName. H :10.pdbName. C :11.pdbName. O :12.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NH2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NH1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11}, {
      {{9}, 8, 9, 1},
      {{12}, 10, 12, 2},
     }},
    {"N", "ASN", MonomerClass::AminoAcid, "NC(=O)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. ND2:1.pdbName. CG :2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "ND2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"D", "ASP", MonomerClass::AminoAcid, "O=C([C@H](CC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. OD1:6.pdbName. OD2:7.pdbName. N :8.pdbName. H :9.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OD1"},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7}, {
      {{8}, 7, 8, 1},
      {{9}, 1, 9, 2},
      {{6}, 4, 6, 3},
     }},
    {"C", "CYS", MonomerClass::AminoAcid, "O=C([C@H](CS[H])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG :5.pdbName. HG :6.pdbName. N :7.pdbName. H :8.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {16, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "SG "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "HG "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6}, {
      {{7}, 6, 7, 1},
      {{8}, 1, 8, 2},
      {{5}, 4, 5, 3},
     }},
    {"Q", "GLN", MonomerClass::AminoAcid, "NC(=O)CC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. NE2:1.pdbName. CD :2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. N :7.pdbName. H :8.pdbName. C :9.pdbName. O :10.pdbName. OXT|",
     {
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 8, 9}, {
      {{7}, 6, 7, 1},
      {{10}, 8, 10, 2},
     }},
    {"E", "GLU", MonomerClass::AminoAcid, "O=C([C@H](CCC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N :9.pdbName. H :10.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OE1"},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OE2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 8}, {
      {{9}, 8, 9, 1},
      {{10}, 1, 10, 2},
      {{7}, 5, 7, 3},
     }},
    {"G", "GLY", MonomerClass::AminoAcid, "O=C(CN[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. N :4.pdbName. H :5.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3}, {
      {{4}, 3, 4, 1},
      {{5}, 1, 5, 2},
     }},
    {"H", "HIS", MonomerClass::AminoAcid, "O=C([C@H](Cc1cnc[nH]1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. ND1:9.pdbName. N :10.pdbName. H :11.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 1, false, "ND1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {
      {{10}, 9, 10, 1},
      {{11}, 1, 11, 2},
     }},
    {"I", "ILE", MonomerClass::AminoAcid, "CC[C@H](C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"L", "LEU", MonomerClass::AminoAcid, "CC(C)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"K", "LYS", MonomerClass::AminoAcid, "O=C([C@H](CCCCN[H])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. HZ1:9.pdbName. N :10.pdbName. H :11.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NZ "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "HZ1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 9}, {
      {{10}, 9, 10, 1},
      {{11}, 1, 11, 2},
      {{8}, 7, 8, 3},
     }},
    {"M", "MET", MonomerClass::AminoAcid, "CSCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CE :1.pdbName. SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N :6.pdbName. H :7.pdbName. C :8.pdbName. O :9.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {16, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "SD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 7, 8}, {
      {{6}, 5, 6, 1},
      {{9}, 7, 9, 2},
     }},
    {"F", "PHE", MonomerClass::AminoAcid, "O=C([C@H](Cc1ccccc1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. CD2:10.pdbName. N :11.pdbName. H :12.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {
      {{11}, 10, 11, 1},
      {{12}, 1, 12, 2},
     }},
    {"P", "PRO", MonomerClass::AminoAcid, "O=C([C@@H]1CCCN1[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD :6.pdbName. N :7.pdbName. H :8.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6}, {
      {{7}, 6, 7, 1},
      {{8}, 1, 8, 2},
     }},
    {"S", "SER", MonomerClass::AminoAcid, "O=C([C@H](CO)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N :6.pdbName. H :7.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OG "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5}, {
      {{6}, 5, 6, 1},
      {{7}, 1, 7, 2},
     }},
    {"T", "THR", MonomerClass::AminoAcid, "C[C@@H](O)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N :5.pdbName. H :6.pdbName. C :7.pdbName. O :8.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(1), 1, true, "CB "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6, 7}, {
      {{5}, 4, 5, 1},
      {{8}, 6, 8, 2},
     }},
    {"W", "TRP", MonomerClass::AminoAcid, "O=C([C@H](Cc1c[nH]c2ccccc12)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. CD2:13.pdbName. N :14.pdbName. H :15.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 1, false, "NE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CH2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ3"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE3"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 14, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 15, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {12, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {12, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, {
      {{14}, 13, 14, 1},
      {{15}, 1, 15, 2},
     }},
    {"Y", "TYR", MonomerClass::AminoAcid, "O=C([C@H](Cc1ccc(O)cc1)N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. CE2:10.pdbName. CD2:11.pdbName. N :12.pdbName. H :13.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CZ "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "OH "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {7, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {7, 9, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
      {2, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 4, static_cast<Bond::BondType>(12), true, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, {
      {{12}, 11, 12, 1},
      {{13}, 1, 13, 2},
     }},
    {"V", "VAL", MonomerClass::AminoAcid, "CC(C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CG1:1.pdbName. CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N :5.pdbName. H :6.pdbName. C :7.pdbName. O :8.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG1"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 6, 7}, {
      {{5}, 4, 5, 1},
      {{8}, 6, 8, 2},
     }},
    {"U", "SEC", MonomerClass::AminoAcid, "O=C([C@H](C[SeH])N[H])[OH] |atomProp:0.pdbName. O :1.pdbName. C :2.pdbName. CA :3.pdbName. CB :4.pdbName. SE :5.pdbName. N :6.pdbName. H :7.pdbName. OXT|",
     {
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {34, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "SE "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 7, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5}, {
      {{6}, 5, 6, 1},
      {{7}, 1, 7, 2},
     }},
    {"O", "PYL", MonomerClass::AminoAcid, "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 :8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG :12.pdbName. CB :13.pdbName. CA :14.pdbName. N :15.pdbName. H :16.pdbName. C :17.pdbName. O :18.pdbName. OXT|",
     {
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CG2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE2"},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N2 "},
      {6, 0, 0, static_cast<Atom::ChiralType>(1), 1, true, "CA2"},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C2 "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O2 "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "NZ "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CE "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CD "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CG "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "CB "},
      {6, 0, 0, static_cast<Atom::ChiralType>(2), 1, true, "CA "},
      {7, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "N "},
      {1, 0, 0, static_cast<Atom::ChiralType>(0), 0, true, "H "},
      {6, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "C "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 0, false, "O "},
      {8, 0, 0, static_cast<Atom::ChiralType>(0), 1, true, "OXT"},
     }, {
      {0, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {1, 2, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {2, 3, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {3, 4, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {4, 5, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 6, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 7, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {6, 8, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {8, 9, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {9, 10, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {10, 11, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {11, 12, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {12, 13, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 14, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {14, 15, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {13, 16, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {16, 17, static_cast<Bond::BondType>(2), false, static_cast<Bond::BondStereo>(0), {}},
      {16, 18, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
      {5, 1, static_cast<Bond::BondType>(1), false, static_cast<Bond::BondStereo>(0), {}},
     }, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17}, {
      {{15}, 14, 15, 1},
      {{18}, 16, 18, 2},
     }},
}};

std::unique_ptr<MacroMolTemplate> makeBuiltinTemplate(
    const BuiltinTemplateDef &def) {
  RWMol mol;
  for (const auto &atomDef : def.atoms) {
    auto *atom = new Atom(atomDef.atomicNum);
    atom->setFormalCharge(atomDef.formalCharge);
    atom->setIsotope(atomDef.isotope);
    atom->setChiralTag(atomDef.chiralTag);
    atom->setNumExplicitHs(atomDef.explicitHs);
    atom->setNoImplicit(atomDef.noImplicit);
    if (atomDef.pdbName) atom->setProp("pdbName", std::string(atomDef.pdbName));
    mol.addAtom(atom, true);
  }
  for (const auto &bondDef : def.bonds) {
    auto *bond = new Bond(bondDef.type);
    bond->setIsAromatic(bondDef.aromatic);
    if (bondDef.stereoAtoms.size() == 2) bond->setStereoAtoms(bondDef.stereoAtoms[0], bondDef.stereoAtoms[1]);
    bond->setStereo(bondDef.stereo);
    bond->setBeginAtomIdx(bondDef.beginAtomIdx);
    bond->setEndAtomIdx(bondDef.endAtomIdx);
    mol.addBond(bond, true);
  }
  MacroMolTemplateBuilder builder(mol, def.monomerClass, def.name, def.symbol, def.originalData);
  builder.setMainGroup(def.mainAtomIdxs);
  for (const auto &leavingGroup : def.leavingGroups) builder.addLeavingGroup(leavingGroup);
  return std::move(builder).build();
}

const std::array<std::unique_ptr<const MacroMolTemplate>, builtinTemplates.size()> &getBuiltinTemplates() {
  static const auto templates = [] {
    std::array<std::unique_ptr<const MacroMolTemplate>, builtinTemplates.size()> result{};
    for (std::size_t i = 0; i < builtinTemplates.size(); ++i) result[i] = makeBuiltinTemplate(builtinTemplates[i]);
    return result;
  }();
  return templates;
}

}  // namespace

void addBuiltinMacroMolTemplates(MacroMolTemplateLibrary &templates) {
  for (const auto &builtinTemplate : getBuiltinTemplates()) templates.addTemplate(std::make_unique<MacroMolTemplate>(*builtinTemplate));
}

MacroMolTemplateLibrary &getGlobalMacroMolTemplateLibrary() {
  static MacroMolTemplateLibrary templates;
  static std::once_flag initialized;
  std::call_once(initialized, []() { addBuiltinMacroMolTemplates(templates); });
  return templates;
}

}  // namespace RDKit
