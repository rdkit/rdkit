/*
 *  sketcherMinimizerClashInteraction.h
 *
 *  Created by Nicola Zonta on 19/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERCLASHMINIMIZERINTERACTION
#define sketcherMINIMIZERCLASHMINIMIZERINTERACTION

#include "sketcherMinimizerInteraction.h"
#include "sketcherMinimizerMaths.h"
#include <iostream>

/* forcefield clash */
class sketcherMinimizerClashInteraction : public sketcherMinimizerInteraction
{
  public:
    sketcherMinimizerClashInteraction(sketcherMinimizerAtom* at1,
                                      sketcherMinimizerAtom* at2,
                                      sketcherMinimizerAtom* at3)
        : sketcherMinimizerInteraction(at1, at2)
    {
        atom3 = at3;
        restV = 900;

        k2 = 0.1f;
    }
    ~sketcherMinimizerClashInteraction() override = default;

    /* calculate the energy of the clash */
    void energy(float& e) override
    {
        squaredDistance = sketcherMinimizerMaths::squaredDistancePointSegment(
            atom2->coordinates, atom1->coordinates, atom3->coordinates);
        if (squaredDistance > restV) {
            return;
        }

        float dr = restV - squaredDistance;
        if (dr > 0) {
            e += 0.5f * k * k2 * dr;
        }
    };

    /* calculate the forces of the clash and apply them */
    void score(float& totalE, bool skipForce = false) override
    {
        energy(totalE);
        if (skipForce) {
            return;
        }
        if (squaredDistance > restV) {
            return;
        }

        sketcherMinimizerPointF atomP = atom2->coordinates;
        sketcherMinimizerPointF bondP1 = atom1->coordinates;
        sketcherMinimizerPointF bondP2 = atom3->coordinates;

        sketcherMinimizerPointF projection =
            sketcherMinimizerMaths::projectPointOnLine(atomP, bondP1, bondP2);
        sketcherMinimizerPointF f = atomP - projection;
        f.normalize();
        f *= (restV - squaredDistance) * k * k2;
        atom2->force += f;
        atom1->force -= f * 0.5;
        atom3->force -= f * 0.5;
    }
    bool isClash() override { return true; };

    float k2;
    sketcherMinimizerAtom* atom3;

  private:
    float squaredDistance;
};

#endif // sketcherMINIMIZERCLASHINTERACTION
