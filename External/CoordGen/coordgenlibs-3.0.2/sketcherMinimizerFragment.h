/*
 *
 *
 *  Created by Nicola Zonta on 03/05/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERFRAGMENT
#define sketcherMINIMIZERFRAGMENT

#include "sketcherMinimizerMaths.h"

#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "CoordgenConfig.hpp"

class sketcherMinimizerAtom;
class sketcherMinimizerBond;
class sketcherMinimizerRing;
class sketcherMinimizerFragment;

/*
 abstract class for fragment degree of freedom
 */
class CoordgenFragmentDOF
{
  public:
    CoordgenFragmentDOF(sketcherMinimizerFragment* fragment);
    virtual ~CoordgenFragmentDOF();

    /* set the current value as the optimal value for this DOF */
    void storeCurrentValueAsOptimal();

    /* load the optimal value */
    void setToOptimalValue();

    /* cycle through the various possible values of this DOF */
    void changeState();

    /* add the given atom to the DoF */
    void addAtom(sketcherMinimizerAtom* atom);

    /* get the penalty associated with the current state */
    virtual float getCurrentPenalty() const;

    /* return the number of states */
    virtual int numberOfStates() const = 0;

    /* return the tier of this DoF. Lower tier DoFs are considered first in the
     * minimization */
    virtual int tier() const = 0;

    /* apply the current DoF value to the atoms coordinates */
    virtual void apply() const = 0;

    /* return the fragment this DoF refers to */
    sketcherMinimizerFragment* getFragment() const;

    /* return the current state */
    short unsigned int getCurrentState();

    /* set the given state as current */
    void setState(short unsigned int state);
    short unsigned int m_currentState, m_optimalState;

  protected:
    std::vector<sketcherMinimizerAtom*> m_atoms;
    sketcherMinimizerFragment* m_fragment;
};

/*
 change the angle of the fragment wrt its parent
 */
class CoordgenRotateFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenRotateFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;
};

/*
 flip the fragment along the bond to its parent
 */
class CoordgenFlipFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenFlipFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;
};

/*
 scale the coordinates of the given atoms
 */
class CoordgenScaleAtomsDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenScaleAtomsDOF(sketcherMinimizerAtom* pivotAtom);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;

  private:
    sketcherMinimizerAtom* m_pivotAtom;
};

/*
 scale the coordinates of the whole fragment
 */
class CoordgenScaleFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenScaleFragmentDOF(sketcherMinimizerFragment* fragment);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;
};

/*
 extend or shorten the bond to the parent fragment
 */
class CoordgenChangeParentBondLengthFragmentDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenChangeParentBondLengthFragmentDOF(
        sketcherMinimizerFragment* fragment);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;
};

/*
 Invert the direction of a bond (e.g. a substituent to a macrocycle can be
 placed towards the inside
 of the ring)
 */
class CoordgenInvertBondDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenInvertBondDOF(sketcherMinimizerAtom* pivotAtom,
                          sketcherMinimizerAtom* boundAtom);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;

  private:
    sketcherMinimizerAtom* m_pivotAtom;
    sketcherMinimizerAtom* m_boundAtom;
};

/* flip a ring fused with another with more than 2 bonds (e.g. ring in a
 * macrocycle) */
class CoordgenFlipRingDOF : public CoordgenFragmentDOF
{
  public:
    CoordgenFlipRingDOF(sketcherMinimizerRing* ring,
                        const std::vector<sketcherMinimizerAtom*>& fusionAtoms);
    int numberOfStates() const override;
    int tier() const override;
    void apply() const override;
    float getCurrentPenalty() const override;

  private:
    sketcherMinimizerAtom* m_pivotAtom1;
    sketcherMinimizerAtom* m_pivotAtom2;
    int m_penalty;
};

/* class that represents a rigid molecular fragment */
class EXPORT_COORDGEN sketcherMinimizerFragment
{
  public:
    sketcherMinimizerFragment();
    ~sketcherMinimizerFragment();

    /* return the total weight of the fragment */
    unsigned int totalWeight() const;

    /* return the number of double bonds in the fragment */
    unsigned int countDoubleBonds() const;

    /* return the number of heavy atoms in the fragment */
    unsigned int countHeavyAtoms() const;

    /* return the number of constrained atoms in the fragment */
    unsigned int countConstrainedAtoms() const;

    /* return the number of fixed atoms in the fragment */
    unsigned int countFixedAtoms() const;

    /* add an atom to this fragment */
    void addAtom(sketcherMinimizerAtom* atom);

    /* add a bond to this fragment */
    void addBond(sketcherMinimizerBond* bond);

    /* add a ring to this fragment */
    void addRing(sketcherMinimizerRing* ring);

    /* add a degree of freedom to this fragment */
    void addDof(CoordgenFragmentDOF* dof);

    /* get all degrees of freedom of this fragment */
    const std::vector<CoordgenFragmentDOF*>& getDofs();

    /* mark the given bond as interfragment */
    void addInterFragmentBond(sketcherMinimizerBond* bond)
    {
        _interFragmentBonds.push_back(bond);
    }

    std::vector<sketcherMinimizerAtom*> getAtoms() const { return m_atoms; }
    std::vector<sketcherMinimizerBond*> getBonds() const { return m_bonds; }
    std::vector<sketcherMinimizerRing*> getRings() const { return m_rings; }

    std::vector<sketcherMinimizerAtom*>& atoms() { return m_atoms; }
    std::vector<sketcherMinimizerBond*>& bonds() { return m_bonds; }
    std::vector<sketcherMinimizerRing*>& rings() { return m_rings; }

    /* return the parent fragment */
    sketcherMinimizerFragment* getParent() const { return m_parent; }

    /* set the given fragment as parent */
    void setParent(sketcherMinimizerFragment* parent) { m_parent = parent; }

    /* add the given degree of freedom to the given atom which will be modified
     * by it */
    void addDofToAtom(sketcherMinimizerAtom* atom, CoordgenFragmentDOF* dof)
    {
        m_dofsForAtom[atom].push_back(dof);
    }

    /* return the degrees of freedom that would modify the given atom */
    std::vector<CoordgenFragmentDOF*>&
    getDofsOfAtom(sketcherMinimizerAtom* atom)
    {
        return m_dofsForAtom[atom];
    }

    /* set the coordinates of each atom to the template coordinates */
    void setAllCoordinatesToTemplate();

    /* save coordinates information */
    void storeCoordinateInformation();

    std::vector<sketcherMinimizerBond*> _interFragmentBonds;
    std::vector<sketcherMinimizerFragment*> _children;
    std::map<sketcherMinimizerAtom*, sketcherMinimizerPointF> _coordinates;
    sketcherMinimizerPointF _bondToParentCoordinatesStart;
    sketcherMinimizerPointF _bondToParentCoordinatesEnd;
    bool fixed, isTemplated, constrained, constrainedFlip;
    bool isChain;
    sketcherMinimizerBond* _bondToParent;
    float longestChainFromHere;
    size_t numberOfChildrenAtoms;
    float numberOfChildrenAtomsRank;

    /* translate and rotate the fragment and set the resulting coordinates to
     * every atom */
    void setCoordinates(const sketcherMinimizerPointF& position, float angle);

    /* get the dof that refers to flipping this fragment */
    CoordgenFragmentDOF* getFlipDof() const { return m_dofs[0]; }

  private:
    sketcherMinimizerFragment* m_parent;
    std::vector<sketcherMinimizerAtom*> m_atoms;
    std::vector<sketcherMinimizerBond*> m_bonds;
    std::vector<sketcherMinimizerRing*> m_rings;
    std::vector<CoordgenFragmentDOF*> m_dofs;
    std::map<sketcherMinimizerAtom*, std::vector<CoordgenFragmentDOF*>>
        m_dofsForAtom;
};

#endif // sketcherMINIMIZERFRAGMENT
