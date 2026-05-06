/*
 *  sketcherMinimizerRing.cpp
 *
 *  Created by Nicola Zonta on 24/05/2011.
 *   Copyright Schrodinger, LLC. All rights reserved.
 *
 */

#include "sketcherMinimizerRing.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMaths.h"

using namespace std;

sketcherMinimizerRing::sketcherMinimizerRing()
    : visited(false), coordinatesGenerated(false)
{

    //    assert (0);
    side = false;
}

sketcherMinimizerRing::~sketcherMinimizerRing() = default;

sketcherMinimizerPointF sketcherMinimizerRing::findCenter()
{
    sketcherMinimizerPointF o(0.f, 0.f);
    for (auto& _atom : _atoms) {
        o += _atom->coordinates;
    }
    o /= _atoms.size();
    return o;
}

bool sketcherMinimizerRing::isBenzene()
{
    if (_atoms.size() != 6) {
        return false;
    }
    for (auto& _atom : _atoms) {
        if (_atom->atomicNumber != 6) {
            return false;
        }
    }
    for (auto a : _atoms) {
        bool found = false;
        for (auto& bond : a->bonds) {
            if (bond->bondOrder == 2) {
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }

    return true;
}

bool sketcherMinimizerRing::isAromatic() // not chemically accurate, but good
                                         // enough for minimizer
{
    size_t bonds = _bonds.size();
    int doubleBonds = 0;
    int NSOCount = 0;
    for (auto& _bond : _bonds) {
        if (_bond->bondOrder == 2) {
            doubleBonds++;
        }
    }
    for (auto& _atom : _atoms) {
        int an = _atom->atomicNumber;
        bool doubleBound = false;
        for (auto& bond : _atom->bonds) {
            if (bond->bondOrder == 2) {
                doubleBound = true;
            }
        }

        if (!doubleBound) {
            if (an == 8 || an == 7 || an == 16) {
                NSOCount++;
            }
        }
    }
    if (bonds == 6 && doubleBonds == 3) {
        return true;
    }
    if (bonds == 5 && doubleBonds == 2 && NSOCount == 1) {
        return true;
    }
    return false;
}

bool sketcherMinimizerRing::containsAtom(const sketcherMinimizerAtom* a) const
{
    for (auto _atom : _atoms) {
        if (_atom == a) {
            return true;
        }
    }
    return false;
}
bool sketcherMinimizerRing::containsBond(sketcherMinimizerBond* b)
{
    for (auto& _bond : _bonds) {
        if (_bond == b) {
            return true;
        }
    }
    return false;
}
bool sketcherMinimizerRing::isFusedWith(sketcherMinimizerRing* ring)
{
    for (auto& i : fusedWith) {
        if (i == ring) {
            return true;
        }
    }
    return false;
}
std::vector<sketcherMinimizerAtom*> sketcherMinimizerRing::getFusionAtomsWith(
    const sketcherMinimizerRing* ring) const
{
    for (unsigned int i = 0; i < fusedWith.size(); i++) {
        if (fusedWith[i] == ring) {
            return fusionAtoms[i];
        }
    }
    std::vector<sketcherMinimizerAtom*> empty;
    return empty;
}

bool sketcherMinimizerRing::sameAs(sketcherMinimizerRing* ring)
{
    if (!(_bonds.size() == ring->_bonds.size())) {
        return false;
    }
    for (auto& _bond : _bonds) {
        if (!ring->containsBond(_bond)) {
            return false;
        }
    }
    return true;
}

bool sketcherMinimizerRing::contains(const sketcherMinimizerPointF& p)
{

    int n = 0;
    for (auto b : _bonds) {
        if ((p.y() < b->startAtom->coordinates.y() &&
             p.y() > b->endAtom->coordinates.y()) ||
            (p.y() > b->startAtom->coordinates.y() &&
             p.y() < b->endAtom->coordinates.y())) {
            sketcherMinimizerPointF v =
                b->endAtom->coordinates - b->startAtom->coordinates;
            if (v.y() > SKETCHER_EPSILON || v.y() < -SKETCHER_EPSILON) {
                v *= (p.y() - b->startAtom->coordinates.y()) / v.y();
                v += b->startAtom->coordinates;
                if (p.x() > v.x()) {
                    n++;
                }
            }
        }
    }
    return (n % 2) != 0;
}

//    std::vector <sketcherMinimizerBond *> _bonds;
