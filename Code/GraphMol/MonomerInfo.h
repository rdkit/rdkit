//
//  Copyright (C) 2013-2024 Greg Landrum and other RDKit contributors
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
#include <RDGeneral/export.h>
#ifndef RD_MONOMERINFO_H
#define RD_MONOMERINFO_H

#include <string>
#include <utility>
#include <boost/shared_ptr.hpp>

namespace RDKit {

//! The abstract base class for atom-level monomer info
class RDKIT_GRAPHMOL_EXPORT AtomMonomerInfo {
 public:
  typedef enum { UNKNOWN = 0, PDBRESIDUE, OTHER, AMINO_ACID, NUCLEIC_ACID } AtomMonomerType;

  virtual ~AtomMonomerInfo() {}

  AtomMonomerInfo() = default;

  /**
   * \param typ the 'type' of monomer this represents:
   *               * UNKNOWN default for base class
   *               * PDBRESIDUE set when AtomPDBResidueInfo is instantiated
   *               * OTHER
   *               * AMINO_ACID
   *               * NUCLEIC_ACID
   * \param nm the atom name (e.g. "CA" using PDB atom naming)
   * \param residueName the residue name (e.g. "ALA")
   * \param resNum the residue number (e.g. 12)
   * \param chainId the chain identifier (e.g. "A", "Heavy", "PEPTIDE1")
   * \param monomerClass further classification of the monomer (e.g, "LGRP", "LINK")
   */
  AtomMonomerInfo(AtomMonomerType typ, std::string nm = "", std::string residueName = "",
                  int resNum = 0, std::string chainId = "", std::string monomerClass = "")
      : d_monomerType(typ), d_name(std::move(nm)), d_residueNumber(resNum),
        d_chainId(std::move(chainId)), d_monomerClass(std::move(monomerClass)),
        d_residueName(std::move(residueName)) {}
  AtomMonomerInfo(const AtomMonomerInfo &other) = default;

  const std::string &getName() const { return d_name; }
  void setName(const std::string &nm) { d_name = nm; }
  AtomMonomerType getMonomerType() const { return d_monomerType; }
  void setMonomerType(AtomMonomerType typ) { d_monomerType = typ; }
  const std::string &getChainId() const { return d_chainId; }
  void setChainId(const std::string &val) { d_chainId = val; }
  int getResidueNumber() const { return d_residueNumber; }
  void setResidueNumber(int val) { d_residueNumber = val; }
  const std::string &getResidueName() const { return d_residueName; }
  void setResidueName(const std::string &val) { d_residueName = val; }
  const std::string &getMonomerClass() const { return d_monomerClass; }
  void setMonomerClass(const std::string &val) { d_monomerClass = val; }

  virtual AtomMonomerInfo *copy() const { return new AtomMonomerInfo(*this); }

 private:
  AtomMonomerType d_monomerType{UNKNOWN};
  std::string d_name{""};
  int d_residueNumber = 0;
  std::string d_chainId = "";
  std::string d_monomerClass = "";
  std::string d_residueName = "";
};

//! Captures atom-level information about peptide residues
class RDKIT_GRAPHMOL_EXPORT AtomPDBResidueInfo : public AtomMonomerInfo {
 public:
  AtomPDBResidueInfo() : AtomMonomerInfo(PDBRESIDUE) {}
  AtomPDBResidueInfo(const AtomPDBResidueInfo &other) = default;

  AtomPDBResidueInfo(const std::string &atomName, int serialNumber = 0,
                     std::string altLoc = "", std::string residueName = "",
                     int residueNumber = 0, std::string chainId = "",
                     std::string insertionCode = "", double occupancy = 1.0,
                     double tempFactor = 0.0, bool isHeteroAtom = false,
                     unsigned int secondaryStructure = 0,
                     unsigned int segmentNumber = 0,
                     std::string monomerClass = "")
      : AtomMonomerInfo(PDBRESIDUE, atomName, residueName, residueNumber, chainId,
                        monomerClass),
        d_serialNumber(serialNumber),
        d_altLoc(std::move(altLoc)),
        d_insertionCode(std::move(insertionCode)),
        d_occupancy(occupancy),
        d_tempFactor(tempFactor),
        df_heteroAtom(isHeteroAtom),
        d_secondaryStructure(secondaryStructure),
        d_segmentNumber(segmentNumber) {}

  int getSerialNumber() const { return d_serialNumber; }
  void setSerialNumber(int val) { d_serialNumber = val; }
  const std::string &getAltLoc() const { return d_altLoc; }
  void setAltLoc(const std::string &val) { d_altLoc = val; }
  const std::string &getInsertionCode() const { return d_insertionCode; }
  void setInsertionCode(const std::string &val) { d_insertionCode = val; }
  double getOccupancy() const { return d_occupancy; }
  void setOccupancy(double val) { d_occupancy = val; }
  double getTempFactor() const { return d_tempFactor; }
  void setTempFactor(double val) { d_tempFactor = val; }
  bool getIsHeteroAtom() const { return df_heteroAtom; }
  void setIsHeteroAtom(bool val) { df_heteroAtom = val; }
  unsigned int getSecondaryStructure() const { return d_secondaryStructure; }
  void setSecondaryStructure(unsigned int val) { d_secondaryStructure = val; }
  unsigned int getSegmentNumber() const { return d_segmentNumber; }
  void setSegmentNumber(unsigned int val) { d_segmentNumber = val; }

  AtomMonomerInfo *copy() const override {
    return static_cast<AtomMonomerInfo *>(new AtomPDBResidueInfo(*this));
  }

 private:
  // the fields here are from the PDB definition
  // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM) [9 Aug, 2013]
  // element and charge are not present since the atom itself stores that
  // information
  unsigned int d_serialNumber = 0;
  std::string d_altLoc = "";
  std::string d_insertionCode = "";
  double d_occupancy = 1.0;
  double d_tempFactor = 0.0;
  // additional, non-PDB fields:
  bool df_heteroAtom = false;  // is this from a HETATM record?
  unsigned int d_secondaryStructure = 0;
  unsigned int d_segmentNumber = 0;
};
};  // namespace RDKit
//! allows AtomPDBResidueInfo objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(
    std::ostream &target, const RDKit::AtomPDBResidueInfo &apri);

#endif
