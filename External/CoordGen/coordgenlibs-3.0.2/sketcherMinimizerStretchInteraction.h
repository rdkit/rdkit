/*
 *  sketcherMinimizerStretchInteraction.h
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERSTRETCHMINIMIZERINTERACTION
#define sketcherMINIMIZERSTRETCHMINIMIZERINTERACTION

#include "sketcherMinimizerInteraction.h"

/* force field bond stretches */
class sketcherMinimizerStretchInteraction : public sketcherMinimizerInteraction
{
  public:
    sketcherMinimizerStretchInteraction(sketcherMinimizerAtom* at1,
                                        sketcherMinimizerAtom* at2)
        : sketcherMinimizerInteraction(at1, at2)
    {
    }
    ~sketcherMinimizerStretchInteraction() override = default;

    /* calculate forces and apply them */
    void score(float& totalE, bool = false) override
    {
        energy(totalE);
        sketcherMinimizerPointF l = atom1->coordinates - atom2->coordinates;
        float m = l.length();
        float dr = 0;
        if (m < restV - tolerance) dr =  restV - tolerance - m;
        else if (m > restV + tolerance) dr =  restV + tolerance - m;
        else return;

        float shortBondThreshold = restV * 0.4f;
        float penaltyForVeryShortBonds = (shortBondThreshold - m);
        if (penaltyForVeryShortBonds < 0) {
            penaltyForVeryShortBonds = 0;
        }
        if (m > SKETCHER_EPSILON) {
            l /= m;
        }
        l *= (k * dr + penaltyForVeryShortBonds * 10);
        atom1->force += l;
        atom2->force -= l;
    }
    float tolerance = 0;
};

#endif // sketcherMINIMIZERSTRETCHINTERACTION
