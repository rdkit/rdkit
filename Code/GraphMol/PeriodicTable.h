//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
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
  class PeriodicTable {

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
    };

    
    //! returns the atomic weight
    double getAtomicWeight( UINT atomicNumber ) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      double mass = byanum[atomicNumber].Mass();
      return mass;
    }
    //! \overload
    double getAtomicWeight( const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element not found");
      int anum = byname.find(elementSymbol)->second;
      double mass = byanum[anum].Mass();
      return mass;
    }
    //! \overload
    double getAtomicWeight( char * elementSymbol ) const {
      return getAtomicWeight(std::string(elementSymbol));
    }

    //! returns the atomic number
    int getAtomicNumber( char *elementSymbol ) const {
      std::string symb(elementSymbol);
      
      return getAtomicNumber(symb);
    }
    //! overload
    int getAtomicNumber( const std::string &elementSymbol ) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      int anum = byname.find(elementSymbol)->second;
      return anum;
    }

    //! returns the atomic symbol
    std::string getElementSymbol(UINT atomicNumber) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].Symbol();
    }

    //! returns the atom's van der Waals radius
    double getRvdw(UINT atomicNumber) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].Rvdw();
    }
    //! \overload
    double getRvdw(const std::string &elementSymbol ) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getRvdw(byname.find(elementSymbol)->second);
    }
    //! \overload
    double getRvdw(char *elementSymbol ) const {
      return getRvdw(std::string(elementSymbol));
    }

    //! returns the atom's covalent radius
    double getRcovalent(UINT atomicNumber) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].Rcov();
    }
    //! \overload
    double getRcovalent(const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getRcovalent(byname.find(elementSymbol)->second);
    }
    //! \overload
    double getRcovalent(char *elementSymbol ) const {
      return getRcovalent(std::string(elementSymbol));
    }

    //! returns the atom's bond radius
    double getRb0(UINT atomicNumber) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].Rb0();
    }
    //! \overload
    double getRb0(const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getRb0(byname.find(elementSymbol)->second);
    }
    //! \overload
    double getRb0(char *elementSymbol ) const {
      return getRb0(std::string(elementSymbol));
    }

    //! returns the atom's default valence 
    int getDefaultValence(UINT atomicNumber) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].DefaultValence();
    }
    //! \overload
    int getDefaultValence(const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getDefaultValence(byname.find(elementSymbol)->second);
    }
    //! \overload
    int getDefaultValence(char *elementSymbol ) const {
      return getDefaultValence(std::string(elementSymbol));
    }

    //! returns a vector of all stable valences. For atoms where
    //! we really don't have any idea what a reasonable maximum
    //! valence is (like transition metals), the vector ends with -1
    const INT_VECT &getValenceList( UINT atomicNumber ) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].ValenceList();
    }
    //! \overload
    const INT_VECT &getValenceList( const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getValenceList(byname.find(elementSymbol)->second);
    }
    //! \overload
    const INT_VECT &getValenceList(char *elementSymbol ) const {
      return getValenceList(std::string(elementSymbol));
    }

    //! returns the number of outer shell electrons
    int getNouterElecs( UINT atomicNumber ) const {
      PRECONDITION(atomicNumber<byanum.size(),"Atomic number not found");
      return byanum[atomicNumber].NumOuterShellElec();
    }
    //! \overload
    int getNouterElecs( const std::string &elementSymbol) const {
      PRECONDITION(byname.count(elementSymbol),"Element '" + elementSymbol +"' not found");
      return getNouterElecs(byname.find(elementSymbol)->second);
    }
    //! \overload
    int getNouterElecs(char *elementSymbol ) const {
      return getNouterElecs(std::string(elementSymbol));
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
      PRECONDITION(anum1<byanum.size(),"Atomic number not found");
      PRECONDITION(anum2<byanum.size(),"Atomic number not found");
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


  private:

    PeriodicTable();
    PeriodicTable &operator =( const PeriodicTable & );

    static class PeriodicTable *ds_instance;

    std::vector<atomicData> byanum;
    STR_UINT_MAP byname;
  };

};

#endif
