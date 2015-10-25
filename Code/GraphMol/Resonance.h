//
//  Copyright (C) 2015 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RESONANCE_H__
#define _RESONANCE_H__

#include <vector>
#include <stack>
#include <map>
#include <boost/unordered_set.hpp>

namespace RDKit {
  class ROMol;
  class Atom;
  class Bond;
  class BondElectrons;
  class AtomElectrons;
  class ConjElectrons;
  class CEVect2;
  typedef std::map<unsigned int, BondElectrons *> ConjBondMap;
  typedef std::map<unsigned int, AtomElectrons *> ConjAtomMap;
  typedef std::vector<ConjElectrons *> CEVect;
  typedef std::vector<CEVect2 *> CEVect3;
  typedef std::vector<boost::uint8_t> ConjFP;
  typedef boost::unordered_map<std::size_t, ConjElectrons *> CEMap;
  class ResonanceMolSupplier
  {
    public:
      typedef enum {
        /*! include resonance structures whose octets are less complete
         *  than the the most octet-complete structure */
        ALLOW_INCOMPLETE_OCTETS = (1 << 0),
        /*! include resonance structures featuring charge separation also
        *   when uncharged resonance structures exist */
        ALLOW_CHARGE_SEPARATION = (1 << 1),
        /*! enumerate all possible degenerate Kekule resonance structures
         *  (the default is to include just one) */
        KEKULE_ALL = (1 << 2),
        /*! if the UNCONSTRAINED_CATIONS flag is not set, positively
         *  charged atoms left and right of N with an incomplete octet are
         *  acceptable only if the conjugated group has a positive total
         *  formal charge */
        UNCONSTRAINED_CATIONS = (1 << 3),
        /*! if the UNCONSTRAINED_ANIONS flag is not set, negatively
         *  charged atoms left of N are acceptable only if the conjugated
         *  group has a negative total formal charge */
        UNCONSTRAINED_ANIONS = (1 << 4)
      } ResonanceFlags;
      /*! 
       *   \param mol        - the starter molecule
       *   \param flags      - flags which influence criteria to generate
       *                       resonance structures
       *   \param maxStructs - maximum number of complete resonance
       *                       structures generated
       */
      ResonanceMolSupplier(ROMol &mol, unsigned int flags = 0,
        unsigned int maxStructs = 1000);
      ~ResonanceMolSupplier();
      /*! Returns a reference to the Kekulized form of the ROMol the
       * ResonanceMolSupplier was initialized with */
      const ROMol &mol() const
      {
        return *d_mol;
      }
      /*! Returns the flags the ResonanceMolSupplier was initialized with
       */
      const unsigned int flags() const
      {
        return d_flags;
      }
      /*! Given a bond index, it returns the index of the conjugated
       *  group the bond belongs to */
      unsigned int getBondConjGrpIdx(unsigned int bi) const;
      /*! Given an atom index, it returns the index of the conjugated
       *  group the atom belongs to */
      unsigned int getAtomConjGrpIdx(unsigned int ai) const;
      /*! Returns the number of resonance structures in the
       *  ResonanceMolSupplier */
      unsigned int length() const
      {
        return d_length;
      };
      /*! Resets the ResonanceMolSupplier index */
      void reset();
      /*! Returns a pointer to the next resonance structure as a ROMol,
       * or NULL if there are no more resonance structures left.
       * The caller is responsible for freeing memory associated to
       * the pointer */
      ROMol *next();
      /*! Returns true if there are no more resonance structures left */
      bool atEnd() const;
      /*! Sets the ResonanceMolSupplier index to idx */
      void moveTo(unsigned int idx);
      /*! Returns a pointer to the resonance structure with index idx as
       *  a ROMol. The index generates complete resonance structures by
       *  combining ConjElectrons objects for the respective conjugated
       *  groups in a breadth-first fashion, in order to return the most
       *  stable complete resonance structures first.
       * The caller is responsible for freeing memory associated to
       * the pointer */
      ROMol *operator[](unsigned int idx) const;
    private:
      typedef struct CEPerm {
        unsigned int idx;
        std::vector<unsigned int> v;
      } CEPerm;
      unsigned int d_nConjGrp;
      unsigned int d_length;
      unsigned int d_flags;
      unsigned int d_maxStructs;
      unsigned int d_idx;
      CEVect3 d_ceVect3;
      void buildCEMap(CEMap &ceMap, unsigned int conjGrpIdx);
      const ROMol *d_mol;
      std::vector<unsigned int> d_bondConjGrpIdx;
      std::vector<unsigned int> d_atomConjGrpIdx;
      std::vector<unsigned int> d_enumIdx;
      // disable copy constructor and assignment operator
      ResonanceMolSupplier(const ResonanceMolSupplier&);
      ResonanceMolSupplier &operator=(const ResonanceMolSupplier&);
      void assignConjGrpIdx();
      void resizeCeVect();
      void trimCeVect2();
      void prepEnumIdxVect();
      void idxToCEPerm(unsigned int idx,
        std::vector<unsigned int> &c) const;
      void storeCEMap(CEMap &ceMap, unsigned int conjGrpIdx);
      void enumerateNbArrangements(CEMap &ceMap, CEMap &ceMapTmp);
      void pruneStructures(CEMap &ceMap);
      void assignBondsFormalChargesHelper(ROMol &mol,
        std::vector<unsigned int> &c) const;
      ROMol *assignBondsFormalCharges(std::vector<unsigned int> &c) const;
      static bool cePermCompare(const CEPerm *a, const CEPerm *b);
  };
}
#endif
