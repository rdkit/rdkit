/*
 *  sketcherMinimizerMolecule.h
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */
#include "CoordgenConfig.hpp"
#include <iostream>
#include <utility>
#include <vector>

#ifndef sketcherMINIMIZERMOLECULE_H
#define sketcherMINIMIZERMOLECULE_H

class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerPointF;

class sketcherMinimizerFragment;

/* class to define a molecule */
class EXPORT_COORDGEN sketcherMinimizerMolecule
{
  public:
    sketcherMinimizerMolecule();
    ~sketcherMinimizerMolecule();

    // create a new atom and add it to the molecule
    sketcherMinimizerAtom* addNewAtom();

    // create a new bond and add it to the molecule
    sketcherMinimizerBond* addNewBond(sketcherMinimizerAtom* at1,
                                      sketcherMinimizerAtom* at2);

    //    void fromChmMol(ChmMol& mol);

    std::vector<sketcherMinimizerAtom*>& getAtoms() { return _atoms; }
    std::vector<sketcherMinimizerBond*>& getBonds() { return _bonds; }
    std::vector<sketcherMinimizerRing*>& getRings() { return _rings; }
    std::vector<sketcherMinimizerFragment*>& getFragments()
    {
        return _fragments;
    }
    void setFragments(std::vector<sketcherMinimizerFragment*> fragments)
    {
        _fragments = std::move(fragments);
    }

    /* set this molecule to require force-field minimization */
    void requireMinimization();

    /* return true if this molecule requires a force field-based minimization
     (i.e. has clashes that cannot be solved otherwise) */
    bool minimizationIsRequired();
    std::vector<sketcherMinimizerAtom*> _atoms;
    std::vector<sketcherMinimizerBond*> _bonds;
    std::vector<sketcherMinimizerRing*> _rings;
    std::vector<sketcherMinimizerBond*> m_proximityRelations;

    std::vector<sketcherMinimizerFragment*> _fragments;

    /* set the given fragment as the main fragment of the kinematic chain */
    void setMainFragment(sketcherMinimizerFragment* fragment)
    {
        m_mainFragment = fragment;
    }
    sketcherMinimizerFragment* getMainFragment() { return m_mainFragment; }

    bool fixed;
    bool hasFixedFragments, hasConstrainedFragments;
    bool needToAlignNonRingAtoms;
    bool needToAlignWholeMolecule;
    bool isPlaced; // used by arrangeMultipleMolecules

    /* return the total charge of the molecule */
    int totalCharge();

    /* set the top left and bottom right points of the molecule's bounding
     * rectangle */
    void boundingBox(sketcherMinimizerPointF& min,
                     sketcherMinimizerPointF& max);

    /* return the coordinates of the center of the molecule */
    sketcherMinimizerPointF center();

    /* recalculate structure elements (e.g. rings) */
    static void forceUpdateStruct(std::vector<sketcherMinimizerAtom*>& atoms,
                                  std::vector<sketcherMinimizerBond*>& bonds,
                                  std::vector<sketcherMinimizerRing*>& rings);

    /* calculate neighbor info of each atom */
    static void
    assignBondsAndNeighbors(std::vector<sketcherMinimizerAtom*>& atoms,
                            std::vector<sketcherMinimizerBond*>& bonds);

    /* run a SSSR algorithm */
    static void findRings(std::vector<sketcherMinimizerBond*>& bonds,
                          std::vector<sketcherMinimizerRing*>& rings);

    /* convenience function for the SSSR algorithm */
    static sketcherMinimizerRing* closeRing(sketcherMinimizerBond* bond);

    /* convenience function for the SSSR algorithm */
    static void addRing(sketcherMinimizerRing* ring,
                        std::vector<sketcherMinimizerRing*>& rings);

  private:
    sketcherMinimizerFragment* m_mainFragment;
    bool m_requireMinimization;
};

#endif // sketcherMINIMIZERMOLECULE_H
