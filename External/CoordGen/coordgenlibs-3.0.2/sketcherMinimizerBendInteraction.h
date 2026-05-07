/*
 *  sketcherMinimizerBendInteraction.h
 *
 *  Created by Nicola Zonta on 13/04/2010.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#ifndef sketcherMINIMIZERBENDMINIMIZERINTERACTION
#define sketcherMINIMIZERBENDMINIMIZERINTERACTION

#include "sketcherMinimizerInteraction.h"

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

/* forcefield class to represent angle bends */
class sketcherMinimizerBendInteraction : public sketcherMinimizerInteraction
{
  public:
    sketcherMinimizerBendInteraction(sketcherMinimizerAtom* at1,
                                     sketcherMinimizerAtom* at2,
                                     sketcherMinimizerAtom* at3)
        : sketcherMinimizerInteraction(at1, at2)
    {
        atom3 = at3;
        restV = 120;
        k2 = 0.05f;
        //    multipleSnap = false;
        isRing = false;
    }
    ~sketcherMinimizerBendInteraction() override = default;

    /* calculate energy associated with the current state */
    void energy(float& e) override
    {
        float dA = angle() - restV;
        e += 0.5f * k * k2 * dA * dA * 10;

        //   qDebug () << restV << "  " << angle ()<<endl;
    };

    /* calculate forces of the interaction */
    void score(float& totalE, bool = false) override
    {
        float a = angle();

        if (a < 0) {
            a = -a;
        }
        float target = restV;
        if (target > 180) {
            target = 360 - target; // this is needed when the angle function is
        }
        // based on cos and only works in [0, 180[ .
        // not needed if using atan2
        /*    if (multipleSnap) {
                vector <int > targets;
                targets .push_back(60);
                targets .push_back(90);
                targets .push_back(120);
            //    targets .push_back(150);

                target = targets [0];
                float distance = target - a;
                if (distance < 0) distance =  -distance;
                for (unsigned int i =1; i < targets.size (); i++) {
                    float newtarget = targets [i];
                    float newdistance = newtarget - a;
                    if (newdistance < 0) newdistance = - newdistance;
                    if (newdistance < distance) {
                        target = newtarget;
                        distance = newdistance;
                    }
                }

            }
             */
        float dA = target - a;
        energy(totalE);
        float x1 = atom1->coordinates.x();
        float y1 = atom1->coordinates.y();
        float x2 = atom2->coordinates.x();
        float y2 = atom2->coordinates.y();
        float x3 = atom3->coordinates.x();
        float y3 = atom3->coordinates.y();

        float v1x = x1 - x2;
        float v1y = y1 - y2;
        float v2x = x3 - x2;
        float v2y = y3 - y2;
        float v3x = x3 - x1;
        float v3y = y3 - y1;

        float newk2 = k2;
        // if (minimizationPhase < 1) newk2 *= 5;
        sketcherMinimizerPointF n1(v1y, -v1x);
        sketcherMinimizerPointF n2(v2y, -v2x);

        if ((n1.x() * v3x + n1.y() * v3y) > 0) {
            n1 *= -1; // dot product n1 v3
        }
        if ((n2.x() * v3x + n2.y() * v3y) < 0) {
            n2 *= -1; // dot product n2 v3
        }

        float q1 = sqrt(n1.x() * n1.x() + n1.y() * n1.y());
        if (q1 < SKETCHER_EPSILON) {
            q1 = SKETCHER_EPSILON;
        }

        float q2 = sqrt(n2.x() * n2.x() + n2.y() * n2.y());
        if (q2 < SKETCHER_EPSILON) {
            q2 = SKETCHER_EPSILON;
        }

        n1 /= q1;
        n2 /= q2;
        n1 *= k * newk2 * dA;
        n2 *= k * newk2 * dA;

        atom1->force += n1;
        atom3->force += n2;
        atom2->force -= n1 + n2;
    };

    /* calculate angle between the three atoms */
    float angle()
    {
        float x1 = atom1->coordinates.x();
        float y1 = atom1->coordinates.y();
        float x2 = atom2->coordinates.x();
        float y2 = atom2->coordinates.y();
        float x3 = atom3->coordinates.x();
        float y3 = atom3->coordinates.y();
        float v1x = x1 - x2;
        float v1y = y1 - y2;
        float v2x = x3 - x2;
        float v2y = y3 - y2;

        float d = sqrt(v1x * v1x + v1y * v1y) * sqrt(v2x * v2x + v2y * v2y);
        if (d < SKETCHER_EPSILON) {
            d = SKETCHER_EPSILON;
        }
        float cosine = (v1x * v2x + v1y * v2y) / d;

        if (cosine < -1) {
            cosine = -1;
        } else if (cosine > 1) {
            cosine = 1;
        }
        return float((acos(cosine)) * 180 / M_PI);
    }
    sketcherMinimizerAtom* atom3;
    float k2;
    //    bool multipleSnap; // used in tetracoordinated centers to get 120 - 90
    //    -90 -60 angles
    bool isRing;
};

#endif // sketcherMINIMIZERBENDINTERACTION
