//
//  Copyright (C) 2018 Greg Landrum
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
namespace RDKit {
namespace MolInterchange {
static const int currentMolJSONVersion = 10;
static const int currentRDKitRepresentationVersion = 1;
static const int currentChargeRepresentationVersion = 10;

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
    {0, Bond::ZERO}, {1, Bond::SINGLE}, {2, Bond::DOUBLE}, {3, Bond::TRIPLE}};
static const std::map<Bond::BondType, unsigned int> inv_bolookup = {
    {Bond::ZERO, 0}, {Bond::SINGLE, 1}, {Bond::DOUBLE, 2}, {Bond::TRIPLE, 3}};

static const std::map<std::string, Bond::BondStereo> stereolookup = {
    {"unspecified", Bond::STEREONONE},
    {"cis", Bond::STEREOCIS},
    {"trans", Bond::STEREOTRANS},
    {"either", Bond::STEREOANY}};
static const std::map<Bond::BondStereo, std::string> inv_stereolookup = {
    {Bond::STEREONONE, "unspecified"}, {Bond::STEREOCIS, "cis"},
    {Bond::STEREOTRANS, "trans"},      {Bond::STEREOZ, "cis"},
    {Bond::STEREOE, "trans"},          {Bond::STEREOANY, "either"}};
}  // namespace MolInterchange
}  // namespace RDKit
#endif
