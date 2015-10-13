// $Id$
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
  typedef std::map<unsigned int, BondElectrons *> ConjBondMap;
  typedef std::map<unsigned int, AtomElectrons *> ConjAtomMap;
  typedef std::vector<std::vector<ConjElectrons *> *> CEVectVect;
  typedef std::vector<boost::uint8_t> ConjFP;
  typedef boost::unordered_map<std::size_t, ConjElectrons *> CEMap;
  class ResonanceMolSupplier
  {
    public:
      typedef enum {
        /*! include resonance structures whose octets are less complete
         *  than the the most octet-complete structure */
        ALLOW_INCOMPLETE_OCTETS = (1 << 0),
        /*! enumerate all possible degenerate Kekule resonance structures
         *  (the default is to include just one) */
        KEKULE_ALL = (1 << 1),
        /*! if the UNCONSTRAINED_CATIONS flag is not set, positively
         *  charged atoms left and right of N with an incomplete octet are
         *  acceptable only if the conjugated group has a positive total
         *  formal charge */
        UNCONSTRAINED_CATIONS = (1 << 2),
        /*! if the UNCONSTRAINED_ANIONS flag is not set, negatively
         *  charged atoms left of N are acceptable only if the conjugated
         *  group has a negative total formal charge */
        UNCONSTRAINED_ANIONS = (1 << 3),
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
      CEVectVect d_ceVectVect;
      CEMap d_ceMapTmp;
      CEMap d_ceMap;
      void buildCEMap(unsigned int conjGrpIdx);
      const ROMol *d_mol;
      std::vector<unsigned int> d_bondConjGrpIdx;
      std::vector<unsigned int> d_atomConjGrpIdx;
      std::vector<unsigned int> d_enumIdx;
      // disable copy constructor and assignment operator
      ResonanceMolSupplier(const ResonanceMolSupplier&);
      ResonanceMolSupplier &operator=(const ResonanceMolSupplier&);
      void assignConjGrpIdx();
      void resizeCeVect();
      void trimCeVectVect();
      void prepEnumIdxVect();
      void idxToCEPerm(unsigned int idx,
        std::vector<unsigned int> &c) const;
      void storeCEMap(unsigned int conjGrpIdx);
      void enumerateNbArrangements();
      void pruneStructures();
      void assignBondsFormalChargesHelper(ROMol &mol,
        std::vector<unsigned int> &c) const;
      ROMol *assignBondsFormalCharges(std::vector<unsigned int> &c) const;
      static bool resonanceStructureCompare
        (const ConjElectrons *a, const ConjElectrons *b);
      static bool cePermCompare(const CEPerm *a, const CEPerm *b);
      static bool vectSizeCompare(const std::vector<ConjElectrons *> *a,
        const std::vector<ConjElectrons *> *b);
  };

  class ConjElectrons {
    public:
      typedef enum {
        HAVE_CATION_RIGHT_OF_N = (1 << 0)
      } ConjElectronsFlags;
      typedef enum {
        FP_BONDS = (1 << 0),
        FP_ATOMS = (1 << 1)
      } FPFlags;
      ConjElectrons(ResonanceMolSupplier *parent,
        unsigned int groupIdx);
      ConjElectrons(const ConjElectrons &ce);
      ~ConjElectrons();
      unsigned int groupIdx() const
      {
        return d_groupIdx;
      };
      unsigned int currElectrons() const
      {
        return d_currElectrons;
      };
      unsigned int totalElectrons() const
      {
        return d_totalElectrons;
      };
      void decrCurrElectrons(unsigned int d);
      AtomElectrons *getAtomElectronsWithIdx(unsigned int ai);
      BondElectrons *getBondElectronsWithIdx(unsigned int bi);
      int getAtomConjGrpIdx(unsigned int ai) const;
      int getBondConjGrpIdx(unsigned int bi) const;
      void pushToBeginStack(unsigned int ai);
      bool popFromBeginStack(unsigned int &ai);
      bool isBeginStackEmpty() {
        return d_beginAIStack.empty();
      };
      unsigned int absFormalCharges() const {
        return d_absFormalCharges;
      };
      unsigned int wtdFormalCharges() const {
        return d_wtdFormalCharges;
      };
      unsigned int lowestFcIndex() const;
      unsigned int lowestMultipleBondIndex() const;
      int allowedChgLeftOfN() const
      {
        return d_allowedChgLeftOfN;
      };
      void decrAllowedChgLeftOfN(int d)
      {
        d_allowedChgLeftOfN -= d;
      };
      int totalFormalCharge() const
      {
        return d_totalFormalCharge;
      };
      unsigned int fcSameSignDist() const
      {
        return d_fcSameSignDist;
      };
      unsigned int fcOppSignDist() const
      {
        return d_fcOppSignDist;
      };
      unsigned int nbMissing() const
      {
        return d_nbMissing;
      }
      bool hasCationRightOfN() const
      {
        return static_cast<bool>(d_flags & HAVE_CATION_RIGHT_OF_N);
      };
      void enumerateNonBonded(CEMap &ceMap);
      void initCeFromMol();
      void assignNonBonded();
      void assignFormalCharge();
      bool assignFormalChargesAndStore
        (CEMap &ceMap, unsigned int fpFlags);
      void assignBondsFormalChargesToMol(ROMol &mol);
      bool checkCharges();
      void computeMetrics();
      bool storeFP(CEMap &ceMap, unsigned int flags);
      ResonanceMolSupplier *parent() const
      {
        return d_parent;
      };
    private:
      unsigned int d_groupIdx;
      unsigned int d_totalElectrons;
      unsigned int d_currElectrons;
      unsigned int d_numFormalCharges;
      int d_totalFormalCharge;
      int d_allowedChgLeftOfN;
      unsigned int d_absFormalCharges;
      unsigned int d_wtdFormalCharges;
      unsigned int d_fcSameSignDist;
      unsigned int d_fcOppSignDist;
      unsigned int d_nbMissing;
      boost::uint8_t d_flags;
      ConjBondMap d_conjBondMap;
      ConjAtomMap d_conjAtomMap;
      std::stack<unsigned int> d_beginAIStack;
      ResonanceMolSupplier *d_parent;
      ConjElectrons &operator=(const ConjElectrons&);
      unsigned int countTotalElectrons();
      void computeDistFormalCharges();
      void checkOctets();
  };
  class AtomElectrons {
    public:
      typedef enum {
        LAST_BOND = (1 << 0),
        DEFINITIVE = (1 << 1),
        STACKED = (1 << 2)
      } AtomElectronsFlags;
      typedef enum {
        NEED_CHARGE_BIT = 1
      } AllowedBondFlag;
      AtomElectrons(ConjElectrons *parent, const Atom *a);
      AtomElectrons(ConjElectrons *parent, const AtomElectrons &ae);
      ~AtomElectrons() {};
      boost::uint8_t findAllowedBonds(unsigned int bi);
      bool hasOctet() const
      {
        return ((d_nb + d_tv * 2) == 8);
      };
      bool isLastBond() const
      {
        return (d_flags & LAST_BOND);
      };
      void setLastBond() {
        d_flags |= LAST_BOND;
      };
      bool isDefinitive() const
      {
        return (d_flags & DEFINITIVE);
      };
      void setDefinitive()
      {
        d_flags |= DEFINITIVE;
      };
      bool isStacked() const
      {
        return (d_flags & STACKED);
      };
      void setStacked() {
        d_flags |= STACKED;
      };
      void clearStacked() {
        d_flags &= ~STACKED;
      };
      unsigned int conjGrpIdx() const
      {
        return d_parent->parent()->getAtomConjGrpIdx(d_atom->getIdx());
      };
      void finalizeAtom();
      unsigned int nb() const
      {
        return d_nb;
      };
      unsigned int tv() const
      {
        return d_tv;
      };
      unsigned int oe() const
      {
        return PeriodicTable::getTable()->getNouterElecs
          (d_atom->getAtomicNum());
      };
      int fc() const
      {
        return d_fc;
      };
      void tvIncr(unsigned int i) {
        d_tv += i;
      };
      unsigned int neededNbForOctet() const
      {
        return (8 - (2 * d_tv + d_nb));
      }
      const Atom *atom()
      {
        return d_atom;
      }
      void initTvNbFcFromAtom();
      void assignNonBonded(unsigned int nb) {
        d_nb = nb;
      }
      void assignFormalCharge() {
        d_fc = oe() - (d_nb + d_tv);
      }
      bool isNbrCharged(unsigned int bo, unsigned int atomicNum = 0);
    private:
      boost::uint8_t d_nb;
      boost::uint8_t d_tv;
      boost::int8_t d_fc;
      boost::uint8_t d_flags;
      const Atom *d_atom;
      ConjElectrons *d_parent;
      AtomElectrons &operator=(const AtomElectrons&);
      boost::uint8_t canAddBondWithOrder(unsigned int bi,
        unsigned int bo);
      void allConjBondsDefinitiveBut(unsigned int bi);
  };
  class BondElectrons {
    public:
      typedef enum {
        DEFINITIVE = (1 << 0)
      } BondElectronsFlags;
      BondElectrons(ConjElectrons *parent, const Bond *b);
      BondElectrons(ConjElectrons *parent, const BondElectrons &be);
      ~BondElectrons() {};
      bool isDefinitive() const
      {
        return (d_flags & DEFINITIVE);
      };
      void setDefinitive() {
        d_flags |= DEFINITIVE;
      };
      int conjGrpIdx() const
      {
        return d_parent->getBondConjGrpIdx(d_bond->getIdx());
      };
      void setOrder(unsigned int bo);
      unsigned int order() const
      {
        return d_bo;
      };
      unsigned int orderFromBond();
      void initOrderFromBond()
      {
        d_bo = orderFromBond();
      };
      const Bond *bond() {
        return d_bond;
      };
    private:
      boost::uint8_t d_bo;
      boost::uint8_t d_flags;
      const Bond *d_bond;
      ConjElectrons *d_parent;
      BondElectrons &operator=(const BondElectrons&);
  };
}
#endif
