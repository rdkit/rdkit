/*
 *  sketcherMinimizerEZConstrainInteraction.h
 *
 *  Created by Nicola Zonta
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZEREZCONSTRAININTERACTION
#define sketcherMINIMIZEREZCONSTRAININTERACTION

#include "sketcherMinimizerInteraction.h"

/* forcefield constrain to avoid EZ inversion */
class sketcherMinimizerEZConstrainInteraction
    : public sketcherMinimizerInteraction
{
  public:
    sketcherMinimizerEZConstrainInteraction(sketcherMinimizerAtom* at1,
                                            sketcherMinimizerAtom* at2,
                                            sketcherMinimizerAtom* at3,
                                            sketcherMinimizerAtom* at4,
                                            bool isZ)
        : sketcherMinimizerInteraction(at1, at2)
    {
        atom3 = at3;
        atom4 = at4;
        m_isZ = isZ;
        m_forceMovement = false;
    };
    ~sketcherMinimizerEZConstrainInteraction() override = default;

    /* calculate the energy of the interaction */
    void energy(float& e) override
    {
        if (inversion()) {
            e += 5000;
        }
    };

    /* calculate the forces and apply them */
    void score(float& totalE, bool = false) override
    {
        if (!inversion()) {
            return;
        }
        energy(totalE);
        sketcherMinimizerPointF projection1 =
            sketcherMinimizerMaths::projectPointOnLine(
                atom1->coordinates, atom2->coordinates, atom3->coordinates);
        sketcherMinimizerPointF projection2 =
            sketcherMinimizerMaths::projectPointOnLine(
                atom4->coordinates, atom2->coordinates, atom3->coordinates);
        sketcherMinimizerAtom* sideAtom = atom1;
        sketcherMinimizerAtom* doubleBondAtom = atom2;
        sketcherMinimizerPointF projection = projection1;
        if (sketcherMinimizerMaths::squaredDistance(atom1->coordinates,
                                                    projection1) >
            sketcherMinimizerMaths::squaredDistance(atom4->coordinates,
                                                    projection2)) {
            sideAtom = atom4;
            doubleBondAtom = atom3;
            projection = projection2;
        }
        sketcherMinimizerPointF force = projection - sideAtom->coordinates;
        if (m_forceMovement) {
            sideAtom->coordinates += force;
            doubleBondAtom->coordinates -= force;
            sideAtom->force = sketcherMinimizerPointF(0, 0);
            doubleBondAtom->force = sketcherMinimizerPointF(0, 0);
        } else {
            force.normalize();
            force *= 10;
            sideAtom->force += force;
            doubleBondAtom->force -= force;
        }
    };

    /* check if the E/Z configuration is inverted */
    bool inversion()
    {
        return sketcherMinimizerMaths::sameSide(
                   atom1->coordinates, atom4->coordinates, atom2->coordinates,
                   atom3->coordinates) != m_isZ;
    }
    sketcherMinimizerAtom* atom3;
    sketcherMinimizerAtom* atom4;
    float k2;
    bool m_isZ;
    bool m_forceMovement;
};

#endif // sketcherMINIMIZEREZCONSTRAININTERACTION
