//
//  Copyright (C) 2001-2011 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_PERIODIC_TABLE_H
#define _RD_PERIODIC_TABLE_H

#include <map>
#include <vector>
#include <RDGeneral/types.h>
#include "atomic_data.h"

namespace RDKit {

//! singleton class for retrieving information about atoms
/*!
  Use the singleton like this:

  \verbatim
  const PeriodicTable *tbl = PeriodicTable::getTable();
  tbl->getAtomicWeight(6); // get atomic weight for Carbon
  tbl->getAtomicWeight("C"); // get atomic weight for Carbon
  \endverbatim

*/
class RDKIT_GRAPHMOL_EXPORT PeriodicTable {
 public:
  //! returns a pointer to the singleton PeriodicTable
  /*
      \return a pointer to the singleton ParamCollection

      <b>Notes:</b>
        - do <b>not</b> delete the pointer returned here
        - if the singleton PeriodicTable has already been instantiated and
          the singleton will be returned, otherwise the singleton will
          be constructed.

   */
  static PeriodicTable *getTable();

  ~PeriodicTable() {
    byanum.clear();
    byname.clear();
  }

  //! returns the atomic weight
  double getAtomicWeight(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    double mass = byanum[atomicNumber].Mass();
    return mass;
  }
  //! \overload
  double getAtomicWeight(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol), "Element not found");
    int anum = byname.find(elementSymbol)->second;
    double mass = byanum[anum].Mass();
    return mass;
  }
  //! \overload
  double getAtomicWeight(const char *elementSymbol) const {
    return getAtomicWeight(std::string(elementSymbol));
  }

  //! returns the atomic number
  int getAtomicNumber(const char *elementSymbol) const {
    std::string symb(elementSymbol);

    return getAtomicNumber(symb);
  }
  //! overload
  int getAtomicNumber(const std::string &elementSymbol) const {
    // this little optimization actually makes a measurable difference
    // in molecule-construction time
    int anum = -1;
    if (elementSymbol == "C") {
      anum = 6;
    } else if (elementSymbol == "N") {
      anum = 7;
    } else if (elementSymbol == "O") {
      anum = 8;
    } else {
      STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
      if (iter != byname.end()) {
        anum = iter->second;
      }
    }
    POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
    return anum;
  }

  //! returns the atomic symbol
  std::string getElementSymbol(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Symbol();
  }

  //! returns the full element name
  std::string getElementName(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Name();
  }

  //! returns the atom's van der Waals radius
  double getRvdw(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Rvdw();
  }
  //! \overload
  double getRvdw(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getRvdw(byname.find(elementSymbol)->second);
  }
  //! \overload
  double getRvdw(const char *elementSymbol) const {
    return getRvdw(std::string(elementSymbol));
  }

  //! returns the atom's covalent radius
  double getRcovalent(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Rcov();
  }
  //! \overload
  double getRcovalent(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getRcovalent(byname.find(elementSymbol)->second);
  }
  //! \overload
  double getRcovalent(const char *elementSymbol) const {
    return getRcovalent(std::string(elementSymbol));
  }

  //! returns the atom's bond radius
  double getRb0(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Rb0();
  }
  //! \overload
  double getRb0(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getRb0(byname.find(elementSymbol)->second);
  }
  //! \overload
  double getRb0(const char *elementSymbol) const {
    return getRb0(std::string(elementSymbol));
  }

  //! returns the atom's default valence
  int getDefaultValence(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].DefaultValence();
  }
  //! \overload
  int getDefaultValence(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getDefaultValence(byname.find(elementSymbol)->second);
  }
  //! \overload
  int getDefaultValence(const char *elementSymbol) const {
    return getDefaultValence(std::string(elementSymbol));
  }

  //! returns a vector of all stable valences. For atoms where
  //! we really don't have any idea what a reasonable maximum
  //! valence is (like transition metals), the vector ends with -1
  const INT_VECT &getValenceList(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].ValenceList();
  }
  //! \overload
  const INT_VECT &getValenceList(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getValenceList(byname.find(elementSymbol)->second);
  }
  //! \overload
  const INT_VECT &getValenceList(const char *elementSymbol) const {
    return getValenceList(std::string(elementSymbol));
  }

  //! returns the number of outer shell electrons
  int getNouterElecs(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].NumOuterShellElec();
  }
  //! \overload
  int getNouterElecs(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getNouterElecs(byname.find(elementSymbol)->second);
  }
  //! \overload
  int getNouterElecs(const char *elementSymbol) const {
    return getNouterElecs(std::string(elementSymbol));
  }

  //! returns the number of the most common isotope
  int getMostCommonIsotope(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].MostCommonIsotope();
  }
  //! \overload
  int getMostCommonIsotope(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getMostCommonIsotope(byname.find(elementSymbol)->second);
  }
  //! \overload
  int getMostCommonIsotope(const char *elementSymbol) const {
    return getMostCommonIsotope(std::string(elementSymbol));
  }

  //! returns the mass of the most common isotope
  double getMostCommonIsotopeMass(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].MostCommonIsotopeMass();
  }
  //! \overload
  double getMostCommonIsotopeMass(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getMostCommonIsotopeMass(byname.find(elementSymbol)->second);
  }
  //! \overload
  double getMostCommonIsotopeMass(const char *elementSymbol) const {
    return getMostCommonIsotopeMass(std::string(elementSymbol));
  }

  //! returns the mass of a particular isotope; zero if that
  //! isotope is unknown.
  double getMassForIsotope(UINT atomicNumber, UINT isotope) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    const std::map<unsigned int, std::pair<double, double>> &m =
        byanum[atomicNumber].d_isotopeInfoMap;
    std::map<unsigned int, std::pair<double, double>>::const_iterator item =
        m.find(isotope);
    if (item == m.end()) {
      return 0.0;
    } else {
      return item->second.first;
    }
  }
  //! \overload
  double getMassForIsotope(const std::string &elementSymbol,
                           UINT isotope) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getMassForIsotope(byname.find(elementSymbol)->second, isotope);
  }
  //! \overload
  double getMassForIsotope(const char *elementSymbol, UINT isotope) const {
    return getMassForIsotope(std::string(elementSymbol), isotope);
  }
  //! returns the abundance of a particular isotope; zero if that
  //! isotope is unknown.
  double getAbundanceForIsotope(UINT atomicNumber, UINT isotope) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    const std::map<unsigned int, std::pair<double, double>> &m =
        byanum[atomicNumber].d_isotopeInfoMap;
    std::map<unsigned int, std::pair<double, double>>::const_iterator item =
        m.find(isotope);
    if (item == m.end()) {
      return 0.0;
    } else {
      return item->second.second;
    }
  }
  //! \overload
  double getAbundanceForIsotope(const std::string &elementSymbol,
                                UINT isotope) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getAbundanceForIsotope(byname.find(elementSymbol)->second, isotope);
  }
  //! \overload
  double getAbundanceForIsotope(const char *elementSymbol, UINT isotope) const {
    return getAbundanceForIsotope(std::string(elementSymbol), isotope);
  }

  //! convenience function to determine which atom is more electronegative
  /*!

     check if atom with atomic number \c anum1 is more
     electronegative than the one with \c anum2
     this is rather lame but here is how we do it
       - the atom with the higher number of outer shell electrons
         is considered more electronegative
       - if the # of outer shell elecs are the same
         the atom with the lower atomic weight is more electronegative

  */
  bool moreElectroNegative(UINT anum1, UINT anum2) const {
    PRECONDITION(anum1 < byanum.size(), "Atomic number not found");
    PRECONDITION(anum2 < byanum.size(), "Atomic number not found");
    // FIX: the atomic_data needs to have real electronegativity values
    UINT ne1 = getNouterElecs(anum1);
    UINT ne2 = getNouterElecs(anum2);
    if (ne1 > ne2) {
      return true;
    }
    if (ne1 == ne2) {
      if (anum1 < anum2) {
        return true;
      }
    }
    return false;
  }

  //! returns the maximum recognized atomic number
  UINT getMaxAtomicNumber() const { return byanum.size() - 1; }
  //! returns the row of the periodic table
  UINT getRow(UINT atomicNumber) const {
    PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    return byanum[atomicNumber].Row();
  }
  //! \overload
  UINT getRow(const std::string &elementSymbol) const {
    PRECONDITION(byname.count(elementSymbol),
                 "Element '" + elementSymbol + "' not found");
    return getRow(byname.find(elementSymbol)->second);
  }
  //! \overload
  UINT getRow(const char *elementSymbol) const {
    return getRow(std::string(elementSymbol));
  }

 private:
  PeriodicTable();
  PeriodicTable &operator=(const PeriodicTable &);
  static void initInstance();

  static class std::unique_ptr<PeriodicTable> ds_instance;

  std::vector<atomicData> byanum;
  STR_UINT_MAP byname;
};
};  // namespace RDKit

#endif
