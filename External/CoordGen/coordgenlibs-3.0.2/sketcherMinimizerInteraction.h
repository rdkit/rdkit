/*
 *  sketcherMinimizerInteraction.h
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERINTERACTION
#define sketcherMINIMIZERINTERACTION

#include "sketcherMinimizerAtom.h"

/* abstract class for force field interactions */
class sketcherMinimizerInteraction
{
  public:
    sketcherMinimizerInteraction(sketcherMinimizerAtom* at1,
                                 sketcherMinimizerAtom* at2)
    {
        atom1 = at1;
        atom2 = at2;
        k = 1.f;
        restV = 50;
        // minimizationPhase = 0;
    };
    virtual ~sketcherMinimizerInteraction() = default;

    /* return energy associated with it */
    virtual void energy(float& e)
    {
        sketcherMinimizerPointF l = atom1->coordinates - atom2->coordinates;
        float dr = sqrt(l.x() * l.x() + l.y() * l.y()) - restV;
        e += 0.5f * k * dr * dr;
    };

    /* calculate and apply forces */
    virtual void score(float& totalE, bool = false)
    {
        sketcherMinimizerPointF l = atom1->coordinates - atom2->coordinates;
        if (l.x() > 0 && l.x() < SKETCHER_EPSILON) {
            l.setX(SKETCHER_EPSILON);
        } else if (l.x() < 0 && l.x() > -SKETCHER_EPSILON) {
            l.setX(-SKETCHER_EPSILON);
        }
        float delta = 0.05f;
        float e1 = 0.f;
        float e2 = 0.f;
        float dx, dy;
        atom1->coordinates.rx() += delta;
        energy(e1);
        atom1->coordinates.rx() -= 2 * delta;
        energy(e2);
        dx = (e2 - e1) / (2 * delta);
        atom1->coordinates.rx() += delta;
        dy = dx * l.y() / l.x();

        totalE += (e2 + e1) * .5f;
        sketcherMinimizerPointF dForce(dx, dy);
        atom1->force += dForce;
        atom2->force -= dForce;
    };
    virtual bool isClash() { return false; };
    float k;
    float restV; // rest value
    sketcherMinimizerAtom* atom1;
    sketcherMinimizerAtom* atom2;
    //    int minimizationPhase;
};

#endif // sketcherMINIMIZERINTERACTION
