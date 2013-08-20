//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MonomerInfo.h

  \brief Defines Monomer information classes

*/  
#ifndef _RD_MONOMERINFO_H
#define _RD_MONOMERINFO_H

#include <string>
#include <boost/shared_ptr.hpp>

namespace RDKit{

  //! The abstract base class for atom-level monomer info
  class AtomMonomerInfo {
  public:
    typedef enum {
      UNKNOWN=0,
      PEPTIDERESIDUE,
      OTHER
    } AtomMonomerType;

    virtual ~AtomMonomerInfo() {};

    AtomMonomerInfo() : d_monomerType(UNKNOWN), d_name("") {};
    AtomMonomerInfo(AtomMonomerType typ,const std::string &nm="") : d_monomerType(typ), d_name(nm)  {};
    AtomMonomerInfo(const AtomMonomerInfo &other) : d_monomerType(other.d_monomerType),d_name(other.d_name)  {};
    
    const std::string & getName() const { return d_name; };
    void setName(const std::string &nm) { d_name=nm; };
    AtomMonomerType getMonomerType() const { return d_monomerType; };
    void setMonomerType(AtomMonomerType typ) { d_monomerType=typ; };

    virtual AtomMonomerInfo *copy() const {
      return new AtomMonomerInfo(*this);
    }
  private:
    AtomMonomerType d_monomerType;
    std::string d_name;
  };

  //! Captures atom-level information about peptide residues
  class AtomPeptideResidueInfo : public AtomMonomerInfo {
  public:
    AtomPeptideResidueInfo() : AtomMonomerInfo(PEPTIDERESIDUE) {};
    AtomPeptideResidueInfo(const AtomPeptideResidueInfo &other) : AtomMonomerInfo(other),
                                                                  d_serialNumber(other.d_serialNumber),
                                                                  d_altLoc(other.d_altLoc),
                                                                  d_residueName(other.d_residueName),
                                                                  d_chainId(other.d_chainId),
                                                                  d_insertionCode(other.d_insertionCode),
                                                                  d_occupancy(other.d_occupancy),
                                                                  d_tempFactor(other.d_tempFactor),
                                                                  d_element(other.d_element),
                                                                  d_charge(other.d_charge) {};

    AtomPeptideResidueInfo(std::string atomName,
                           unsigned int serialNumber=0,
                           std::string altLoc="",
                           std::string residueName="",
                           std::string chainId="",
                           std::string insertionCode="",
                           double occupancy=1.0,
                           double tempFactor=0.0,
                           std::string element="",
                           std::string charge="") : AtomMonomerInfo(PEPTIDERESIDUE,atomName),
                                                    d_serialNumber(serialNumber),
                                                    d_altLoc(altLoc),
                                                    d_residueName(residueName),
                                                    d_chainId(chainId),
                                                    d_insertionCode(insertionCode),
                                                    d_occupancy(occupancy),
                                                    d_tempFactor(tempFactor),
                                                    d_element(element),
                                                    d_charge(charge) {};
    
                           
    
    
    unsigned int getSerialNumber() const { return d_serialNumber; };
    void setSerialNumber(unsigned int val) { d_serialNumber=val; };
    const std::string &getAltLoc() const { return d_altLoc; };
    void setAltLoc(const std::string &val) { d_altLoc=val; };
    const std::string &getResidueName() const { return d_residueName; };
    void setResidueName(const std::string &val) { d_residueName=val; };
    const std::string &getChainId() const { return d_chainId; };
    void setChainId(const std::string &val) { d_chainId=val; };
    const std::string &getInsertionCode() const { return d_insertionCode; };
    void setInsertionCode(const std::string &val) { d_insertionCode=val; };
    double getOccupancy() const { return d_occupancy; };
    void setOccupancy(double val) { d_occupancy=val; };
    double getTempFactor() const { return d_tempFactor; };
    void setTempFactor(double val) { d_tempFactor=val; };
    const std::string &getElement() const { return d_element; };
    void setElement(const std::string &val) { d_element=val; };
    const std::string &getCharge() const { return d_charge; };
    void setCharge(const std::string &val) { d_charge=val; };

    AtomMonomerInfo *copy() const {
      return static_cast<AtomMonomerInfo *>(new AtomPeptideResidueInfo(*this));
    }
    
  private:
    // the fields here are from the PDB definition 
    // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM) [9 Aug, 2013]
    unsigned int d_serialNumber;
    std::string d_altLoc;
    std::string d_residueName;
    std::string d_chainId;
    std::string d_insertionCode;
    double d_occupancy;
    double d_tempFactor;
    std::string d_element;
    std::string d_charge;
  };
    
};
//! allows AtomPeptideResidueInfo objects to be dumped to streams
std::ostream & operator<<(std::ostream& target, const RDKit::AtomPeptideResidueInfo &apri);

#endif
