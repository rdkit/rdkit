//
//  Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file atomic_data.h

  \brief No user-serviceable parts inside

  This stuff is used by the PeriodicTable interface

*/
#pragma once
#include <RDGeneral/export.h>
#ifndef RD_ATOMIC_DATA_H
#define RD_ATOMIC_DATA_H

#include <vector>
#include <cstdint>
#include <string>
#include <map>

namespace RDKit {

// Names of elements ordered by atomic number
constexpr inline std::uint8_t nElements = 119;
constexpr inline std::array<const char *, nElements> elementNames = {
    "Dummy",         "Hydrogen",    "Helium",       "Lithium",
    "Beryllium",     "Boron",       "Carbon",       "Nitrogen",
    "Oxygen",        "Fluorine",    "Neon",         "Sodium",
    "Magnesium",     "Aluminium",   "Silicon",      "Phosphorus",
    "Sulfur",        "Chlorine",    "Argon",        "Potassium",
    "Calcium",       "Scandium",    "Titanium",     "Vanadium",
    "Chromium",      "Manganese",   "Iron",         "Cobalt",
    "Nickel",        "Copper",      "Zinc",         "Gallium",
    "Germanium",     "Arsenic",     "Selenium",     "Bromine",
    "Krypton",       "Rubidium",    "Strontium",    "Yttrium",
    "Zirconium",     "Niobium",     "Molybdenum",   "Technetium",
    "Ruthenium",     "Rhodium",     "Palladium",    "Silver",
    "Cadmium",       "Indium",      "Tin",          "Antimony",
    "Tellurium",     "Iodine",      "Xenon",        "Caesium",
    "Barium",        "Lanthanum",   "Cerium",       "Praseodymium",
    "Neodymium",     "Promethium",  "Samarium",     "Europium",
    "Gadolinium",    "Terbium",     "Dysprosium",   "Holmium",
    "Erbium",        "Thulium",     "Ytterbium",    "Lutetium",
    "Hafnium",       "Tantalum",    "Tungsten",     "Rhenium",
    "Osmium",        "Iridium",     "Platinum",     "Gold",
    "Mercury",       "Thallium",    "Lead",         "Bismuth",
    "Polonium",      "Astatine",    "Radon",        "Francium",
    "Radium",        "Actinium",    "Thorium",      "Protactinium",
    "Uranium",       "Neptunium",   "Plutonium",    "Americium",
    "Curium",        "Berkelium",   "Californium",  "Einsteinium",
    "Fermium",       "Mendelevium", "Nobelium",     "Lawrencium",
    "Rutherfordium", "Dubnium",     "Seaborgium",   "Bohrium",
    "Hassium",       "Meitnerium",  "Darmstadtium", "Roentgenium",
    "Copernicium",   "Nihonium",    "Flerovium",    "Moscovium",
    "Livermorium",   "Tennessine",  "Oganesson"};

class RDKIT_GRAPHMOL_EXPORT atomicData {
 public:
  atomicData(int atnum, const char *symb, unsigned int row, double rCov,
             double rB0, double rVdw, double mass, int nVal, int commonIsotope,
             double commonIsotopeMass, const std::vector<int> &valences)
      : anum(atnum),
        symb(symb),
        name(elementNames[atnum]),
        rCov(rCov),
        rB0(rB0),
        rVdw(rVdw),
        valence(valences),
        mass(mass),
        nVal(nVal),
        commonIsotope(commonIsotope),
        commonIsotopeMass(commonIsotopeMass),
        row(row) {};
  ~atomicData() = default;

  int AtomicNum() const { return anum; }

  int DefaultValence() const { return valence.front(); }

  int NumValence() const { return static_cast<int>(valence.size()); }

  const std::vector<int> &ValenceList() const { return valence; }

  double Mass() const { return mass; }

  std::string Symbol() const { return symb; }

  std::string Name() const { return name; }

  double Rcov() const { return rCov; }

  double Rb0() const { return rB0; }

  double Rvdw() const { return rVdw; }

  int NumOuterShellElec() const { return nVal; }

  int MostCommonIsotope() const { return commonIsotope; }

  double MostCommonIsotopeMass() const { return commonIsotopeMass; }

  unsigned int Row() const { return row; }
  // maps isotope number -> mass
  std::map<unsigned int, std::pair<double, double>>
      d_isotopeInfoMap;  // available isotopes
 private:
  int anum;                // atomic number
  std::string symb;        // atomic symbol
  std::string name;        // atomic name
  double rCov, rB0, rVdw;  // radii
  std::vector<int>
      valence;  // list of all valences, the first one is the default
  // valence, -1 at the end signifies that any upper valence
  // is tolerated
  double mass;               // atomic mass
  int nVal;                  // number of outer shell electrons
  int commonIsotope;         // most common isotope
  double commonIsotopeMass;  // most common isotopic mass
  unsigned int row;          // row in the periodic table
};

struct RDKIT_GRAPHMOL_EXPORT isotopeInfo {
  unsigned int atomicNumber;
  std::string elementSymbol;
  unsigned int isotope;
  double mass;
  double abundance;
};

extern const std::vector<atomicData> atomicDataVals;
extern const std::vector<isotopeInfo> isotopeDataVals;

}  // namespace RDKit
#endif
