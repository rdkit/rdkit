//
//  Copyright (C) 2015 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RESONANCE_H__
#define _RESONANCE_H__

#include <vector>
#include <stack>
#include <map>
#include <unordered_map>

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
typedef std::vector<std::uint8_t> ConjFP;
typedef std::unordered_map<std::size_t, ConjElectrons *> CEMap;

/*!
 * Create a derived class from this abstract base class
 * and implement the operator()() method.
 * The operator()() method is called at each iteration of the
 * algorithm, and provides a mechanism to monitor or stop
 * its progress.
 * To have your callback called, pass an instance of your
 * derived class to ResonanceMolSupplier::setProgressCallback().
 */
class RDKIT_GRAPHMOL_EXPORT ResonanceMolSupplierCallback {
  friend class ResonanceMolSupplier;

 public:
  ResonanceMolSupplierCallback() {}
  virtual ~ResonanceMolSupplierCallback() {}
  /*! Get the number of conjugated groups this molecule has. */
  unsigned int getNumConjGrps() const { return d_progress.size(); }
  /*! Get the maximum number of resonance structures the supplier
      is allowed to generate.
      \return The number of structures.
   */
  size_t getMaxStructures() const { return d_maxStructs; }
  /*! Get the number of resonance structures generated so far
      for a certain conjugated group.
      \param conjGrpIdx the conjugated group index.
      \return the number of resonance structures generated so far
              for this group.
   */
  size_t getNumStructures(unsigned int conjGrpIdx) const;
  /*! Get the number of resonance structures with a different score
      (i.e., non-degenerate) generated so far for a certain conjugated
      group.
      \param conjGrpIdx the conjugated group index.
      \return the number of non-degenrate resonance structures generated
              so far for this group.
   */
  size_t getNumDiverseStructures(unsigned int conjGrpIdx) const;
  /*! Pure virtual; this must be implemented in the derived class.
      \return true if the resonance structure generation should continue;
              false if the resonance structure generation should stop.
   */
  virtual bool operator()() = 0;

 private:
  struct ResonanceProgress {
    size_t d_totalStructs;
    size_t d_diverseStructs;
  };
  size_t d_maxStructs{0};
  std::vector<ResonanceProgress> d_progress;
};

class RDKIT_GRAPHMOL_EXPORT ResonanceMolSupplier {
 public:
  typedef enum {
    /*! include resonance structures whose octets are less complete
     *  than the most octet-complete structure */
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
     *  formal charge.
     *  UNCONSTRAINED_CATIONS implies ALLOW_INCOMPLETE_OCTETS
     *  and ALLOW_CHARGE_SEPARATION.
     */
    UNCONSTRAINED_CATIONS = (1 << 3),
    /*! if the UNCONSTRAINED_ANIONS flag is not set, negatively
     *  charged atoms left of N are acceptable only if the conjugated
     *  group has a negative total formal charge.
     *  UNCONSTRAINED_ANIONS implies ALLOW_CHARGE_SEPARATION.
     */
    UNCONSTRAINED_ANIONS = (1 << 4)
  } ResonanceFlags;
  /*!
   *   \param mol        - the starter molecule
   *   \param flags      - flags which influence criteria to generate
   *                       resonance structures
   *   \param maxStructs - maximum number of complete resonance
   *                       structures generated
   *   \param numThreads - the number of threads used to carry out the
   *                       resonance structure enumeration (defaults
   *                       to 1; 0 selects the number of concurrent
   *                       threads supported by the hardware; negative
   *                       values are added to the number of
   *                       concurrent threads supported by the
   *                       hardware)
   */
  ResonanceMolSupplier(ROMol &mol, unsigned int flags = 0,
                       unsigned int maxStructs = 1000);
  ~ResonanceMolSupplier();
  /*! Returns a reference to the Kekulized form of the ROMol the
   *  ResonanceMolSupplier was initialized with */
  const ROMol &mol() const { return *d_mol; }
  /*! Returns the flags the ResonanceMolSupplier was initialized with
   */
  unsigned int flags() const { return d_flags; }
  /*! Returns the number of individual conjugated groups
      in the molecule */
  unsigned int getNumConjGrps() const { return d_nConjGrp; }
  /*! Given a bond index, it returns the index of the conjugated
   *  group the bond belongs to, or -1 if it is not conjugated */
  int getBondConjGrpIdx(unsigned int bi) const;
  /*! Given an atom index, it returns the index of the conjugated
   *  group the atom belongs to, or -1 if it is not conjugated */
  int getAtomConjGrpIdx(unsigned int ai) const;
  /*! Sets the number of threads to be used to enumerate resonance
   *  structures (defaults to 1; 0 selects the number of concurrent
   *  threads supported by the hardware; negative values are added
   *  to the number of concurrent threads supported by the hardware)
   */
  void setNumThreads(int numThreads = 1);
  /*! Pass a pointer to an instance of a class derived from
   *  ResonanceMolSupplierCallback, which must implement the callback()
   *  method. The ResonanceMolSupplier takes ownership of the pointer,
   *  which must not be deleted.
   */
  void setProgressCallback(ResonanceMolSupplierCallback *callback);
  /*! Get the pointer to the ResonanceMolSupplierCallback subclass
   *  instance (do not delete it, has the ResonanceMolSupplier takes
   *  ownership of the pointer), or nullptr if none was set.
   */
  ResonanceMolSupplierCallback *getProgressCallback() const {
    return d_callback.get();
  }
  /*! Returns true if the resonance structure generation was canceled
   *  before completion due to the ResMolSupplierProgressCallback
   *  having returned false.
   */
  bool wasCanceled() const { return d_wasCanceled; }
  /*! Ask ResonanceMolSupplier to enumerate resonance structures
   *  (automatically done as soon as any attempt to access them is
   *  made) */
  void enumerate();
  /*! Returns true if resonance structure enumeration has already
   *  happened */
  bool getIsEnumerated() { return d_isEnumerated; }
  /*! Returns the number of resonance structures in the
   *  ResonanceMolSupplier */
  unsigned int length();
  /*! Resets the ResonanceMolSupplier index */
  void reset();
  /*! Returns true if there are no more resonance structures left */
  bool atEnd();
  /*! Returns a pointer to the next resonance structure as a ROMol,
   * or NULL if there are no more resonance structures left.
   * The caller is responsible for freeing memory associated to
   * the pointer */
  ROMol *next();
  /*! Sets the ResonanceMolSupplier index to idx */
  void moveTo(unsigned int idx);
  /*! Returns a pointer to the resonance structure with index idx as
   *  a ROMol. The index generates complete resonance structures by
   *  combining ConjElectrons objects for the respective conjugated
   *  groups in a breadth-first fashion, in order to return the most
   *  stable complete resonance structures first.
   * The caller is responsible for freeing memory associated to
   * the pointer */
  ROMol *operator[](unsigned int idx);

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
  unsigned int d_numThreads;
  bool d_isEnumerated;
  bool d_wasCanceled;
  CEVect3 d_ceVect3;
  const ROMol *d_mol;
  std::vector<int> d_bondConjGrpIdx;
  std::vector<int> d_atomConjGrpIdx;
  std::vector<size_t> d_enumIdx;
  std::unique_ptr<ResonanceMolSupplierCallback> d_callback;
  // disable copy constructor and assignment operator
  ResonanceMolSupplier(const ResonanceMolSupplier &);
  ResonanceMolSupplier &operator=(const ResonanceMolSupplier &);
  void mainLoop(unsigned int ti, unsigned int nt);
  void assignConjGrpIdx();
  void resizeCeVect();
  void trimCeVect2();
  void prepEnumIdxVect();
  void idxToCEPerm(size_t idx, std::vector<unsigned int> &c) const;
  void setResonanceMolSupplierLength();
  void buildCEMap(CEMap &ceMap, unsigned int conjGrpIdx);
  void storeCEMap(const CEMap &ceMap, unsigned int conjGrpIdx);
  void assignBondsFormalChargesHelper(ROMol &mol,
                                      std::vector<unsigned int> &c) const;
  ROMol *assignBondsFormalCharges(std::vector<unsigned int> &c) const;
  static bool cePermCompare(const CEPerm *a, const CEPerm *b);
};
}  // namespace RDKit
#endif
