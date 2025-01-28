//
//  Copyright (C) 2018-2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLINTERCHANGEDETAILS_H_FEB2018
#define RD_MOLINTERCHANGEDETAILS_H_FEB2018
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/StereoGroup.h>
namespace RDKit {
namespace MolInterchange {
constexpr int currentMolJSONVersion = 10;
constexpr int currentRDKitJSONVersion = 11;
constexpr int currentRDKitRepresentationVersion = 2;
constexpr int currentChargeRepresentationVersion = 10;
constexpr int currentQueryRepresentationVersion = 10;

static const std::map<std::string, Atom::ChiralType> chilookup = {
    {"unspecified", Atom::CHI_UNSPECIFIED},
    {"cw", Atom::CHI_TETRAHEDRAL_CW},
    {"ccw", Atom::CHI_TETRAHEDRAL_CCW},
    {"other", Atom::CHI_OTHER}};
static const std::map<Atom::ChiralType, std::string> inv_chilookup = {
    {Atom::CHI_UNSPECIFIED, "unspecified"},
    {Atom::CHI_TETRAHEDRAL_CW, "cw"},
    {Atom::CHI_TETRAHEDRAL_CCW, "ccw"},
    {Atom::CHI_OTHER, "other"}};

static const std::map<unsigned int, Bond::BondType> bolookup = {
    {0, Bond::ZERO},   {1, Bond::SINGLE},    {2, Bond::DOUBLE},
    {3, Bond::TRIPLE}, {4, Bond::QUADRUPLE}, {17, Bond::DATIVE}};
static const std::map<Bond::BondType, unsigned int> inv_bolookup = {
    {Bond::ZERO, 0},   {Bond::SINGLE, 1},    {Bond::DOUBLE, 2},
    {Bond::TRIPLE, 3}, {Bond::QUADRUPLE, 4}, {Bond::DATIVE, 17}};

static const std::map<std::string, Bond::BondStereo> stereoBondlookup = {
    {"unspecified", Bond::STEREONONE},
    {"cis", Bond::STEREOCIS},
    {"trans", Bond::STEREOTRANS},
    {"either", Bond::STEREOANY}};
static const std::map<Bond::BondStereo, std::string> inv_stereoBondlookup = {
    {Bond::STEREONONE, "unspecified"}, {Bond::STEREOCIS, "cis"},
    {Bond::STEREOTRANS, "trans"},      {Bond::STEREOZ, "cis"},
    {Bond::STEREOE, "trans"},          {Bond::STEREOANY, "either"}};

static const std::map<std::string, StereoGroupType> stereoGrouplookup = {
    {"abs", StereoGroupType::STEREO_ABSOLUTE},
    {"and", StereoGroupType::STEREO_AND},
    {"or", StereoGroupType::STEREO_OR}};
static const std::map<StereoGroupType, std::string> inv_stereoGrouplookup = {
    {StereoGroupType::STEREO_ABSOLUTE, "abs"},
    {StereoGroupType::STEREO_AND, "and"},
    {StereoGroupType::STEREO_OR, "or"}};

}  // namespace MolInterchange
}  // namespace RDKit
#endif
