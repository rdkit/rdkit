/*
 Contributors: Nicola Zonta
 Copyright Schrodinger, LLC. All rights reserved
 */

#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMaths.h"
#include "sketcherMinimizerRing.h"

static const float INVERTED_MACROCYCLE_BOND_PENALTY = 100.f;
static const float SCALE_FRAGMENT_PENALTY = 500.f;
static const float SCALE_ATOMS_PENALTY = 50.f;
static const float ROTATE_FRAGMENT_PENALTY = 400.f;
static const float BREAK_CHAIN_PENALTY = 10.f;
static const float CHANGE_PARENT_BOND_PENALTY = 200.f;
static const float FLIP_RING_PENALTY = 200.f;
static const float FLIP_CONSTRAINED_FRAGMENT_PENALTY = 1000.f;

static const int FLIP_FRAGMENT_TIER = 0;
static const int INVERT_BOND_TIER = 1;
static const int ROTATE_FRAGMENT_TIER = 3;
static const int ELONGATE_PARENT_BOND_TIER = 2;
static const int SCALE_ATOMS_TIER = 4;
static const int FLIP_RING_TIER = 1;
static const int SCALE_FRAGMENT_TIER = 5;

CoordgenFragmentDOF::CoordgenFragmentDOF(sketcherMinimizerFragment* fragment)
    : m_currentState(0), m_optimalState(0), m_fragment(fragment)
{
}

CoordgenFragmentDOF::~CoordgenFragmentDOF() = default;

short unsigned int CoordgenFragmentDOF::getCurrentState()
{
    return m_currentState;
}

void CoordgenFragmentDOF::setState(short unsigned int state)
{
    m_currentState = state;
}

void CoordgenFragmentDOF::storeCurrentValueAsOptimal()
{
    m_optimalState = m_currentState;
}

void CoordgenFragmentDOF::setToOptimalValue()
{
    m_currentState = m_optimalState;
}

void CoordgenFragmentDOF::changeState()
{
    m_currentState++;
    m_currentState = m_currentState % numberOfStates();
}

void CoordgenFragmentDOF::addAtom(sketcherMinimizerAtom* atom)
{
    m_atoms.push_back(atom);
    atom->fragment->addDofToAtom(atom, this);
}

float CoordgenFragmentDOF::getCurrentPenalty() const
{
    return 0.f;
}

sketcherMinimizerFragment* CoordgenFragmentDOF::getFragment() const
{
    return m_fragment;
}

CoordgenFlipFragmentDOF::CoordgenFlipFragmentDOF(
    sketcherMinimizerFragment* fragment)
    : CoordgenFragmentDOF(fragment)
{
}

float CoordgenFlipFragmentDOF::getCurrentPenalty() const
{
    float penalty = 0.f;
    if (m_currentState != 0 && m_fragment->constrainedFlip) {
        penalty += FLIP_CONSTRAINED_FRAGMENT_PENALTY;
    }
    if (m_fragment->isChain && m_fragment->getParent() &&
        m_fragment->getParent()->isChain) {
        penalty += BREAK_CHAIN_PENALTY;
    }
    return penalty;
}

int CoordgenFlipFragmentDOF::numberOfStates() const
{
    if (m_fragment->getParent() == nullptr) {
        return 1;
    }
    return 2;
}

int CoordgenFlipFragmentDOF::tier() const
{
    return FLIP_FRAGMENT_TIER;
}

void CoordgenFlipFragmentDOF::apply() const
{
    if (m_currentState != 0) {
        for (auto& atom : m_fragment->_coordinates) {
            atom.first->coordinates.setY(-atom.first->coordinates.y());
        }
    }
}

CoordgenScaleFragmentDOF::CoordgenScaleFragmentDOF(
    sketcherMinimizerFragment* fragment)
    : CoordgenFragmentDOF(fragment)
{
}

int CoordgenScaleFragmentDOF::numberOfStates() const
{
    if (m_fragment->getRings().empty()) {
        return 1;
    }
    return 5;
}

int CoordgenScaleFragmentDOF::tier() const
{
    return SCALE_FRAGMENT_TIER;
}

void CoordgenScaleFragmentDOF::apply() const
{
    if (m_currentState != 0) {
        auto scale = static_cast<float>(pow(1.4, (m_currentState + 1) / 2));
        if (m_currentState % 2 == 0) {
            scale = 1 / scale;
        }
        for (auto& atom : m_fragment->_coordinates) {
            atom.first->setCoordinates(atom.first->getCoordinates() * scale);
        }
    }
}

float CoordgenScaleFragmentDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return SCALE_FRAGMENT_PENALTY * ((m_currentState + 1) / 2);
    }
    return 0.f;
}

CoordgenScaleAtomsDOF::CoordgenScaleAtomsDOF(sketcherMinimizerAtom* pivotAtom)
    : CoordgenFragmentDOF(pivotAtom->getFragment()), m_pivotAtom(pivotAtom)
{
}

int CoordgenScaleAtomsDOF::numberOfStates() const
{
    return 2;
}

int CoordgenScaleAtomsDOF::tier() const
{
    return SCALE_ATOMS_TIER;
}

void CoordgenScaleAtomsDOF::apply() const
{
    if (m_currentState != 0) {
        float scale = 0.4f;
        for (auto atom : m_atoms) {
            auto distance =
                atom->getCoordinates() - m_pivotAtom->getCoordinates();
            atom->setCoordinates(distance * scale +
                                 m_pivotAtom->getCoordinates());
        }
    }
}

float CoordgenScaleAtomsDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return SCALE_ATOMS_PENALTY * m_atoms.size();
    }
    return 0.f;
}

CoordgenChangeParentBondLengthFragmentDOF::
    CoordgenChangeParentBondLengthFragmentDOF(
        sketcherMinimizerFragment* fragment)
    : CoordgenFragmentDOF(fragment)
{
}

int CoordgenChangeParentBondLengthFragmentDOF::numberOfStates() const
{
    return 7;
}

int CoordgenChangeParentBondLengthFragmentDOF::tier() const
{
    return ELONGATE_PARENT_BOND_TIER;
}

void CoordgenChangeParentBondLengthFragmentDOF::apply() const
{
    if (m_currentState != 0) {
        auto scale = static_cast<float>(pow(1.6, (m_currentState + 1) / 2));
        if (m_currentState % 2 == 0) {
            scale = 1 / scale;
        }
        float MoveBy = BONDLENGTH * (scale - 1);
        for (auto& atom : m_fragment->_coordinates) {
            atom.first->coordinates.setX(atom.first->coordinates.x() + MoveBy);
        }
    }
}

float CoordgenChangeParentBondLengthFragmentDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return CHANGE_PARENT_BOND_PENALTY * ((m_currentState + 1) / 2);
    }
    return 0.f;
}

CoordgenRotateFragmentDOF::CoordgenRotateFragmentDOF(
    sketcherMinimizerFragment* fragment)
    : CoordgenFragmentDOF(fragment)
{
}

int CoordgenRotateFragmentDOF::numberOfStates() const
{
    if (m_fragment->getParent() == nullptr) {
        return 1;
    }
    return 5;
}

int CoordgenRotateFragmentDOF::tier() const
{
    return ROTATE_FRAGMENT_TIER;
}

void CoordgenRotateFragmentDOF::apply() const
{
    if (m_currentState != 0) {
        auto angle =
            static_cast<float>(M_PI / 180 * 15 * ((m_currentState + 1) / 2));
        if (m_currentState % 2 == 0) {
            angle = -angle;
        }
        float sine = sin(angle);
        float cosine = cos(angle);
        for (auto& atom : m_fragment->_coordinates) {
            sketcherMinimizerPointF origin(-BONDLENGTH, 0);
            sketcherMinimizerPointF coords =
                atom.first->getCoordinates() - origin;
            coords.rotate(sine, cosine);
            atom.first->setCoordinates(coords + origin);
        }
    }
}

float CoordgenRotateFragmentDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return ROTATE_FRAGMENT_PENALTY * ((m_currentState + 1) / 2);
    }
    return 0.f;
}

CoordgenInvertBondDOF::CoordgenInvertBondDOF(sketcherMinimizerAtom* pivotAtom,
                                             sketcherMinimizerAtom* boundAtom)
    : CoordgenFragmentDOF(pivotAtom->getFragment()), m_pivotAtom(pivotAtom),
      m_boundAtom(boundAtom)
{
    assert(pivotAtom->bondTo(boundAtom) != nullptr);
    addAtom(boundAtom);
}

int CoordgenInvertBondDOF::numberOfStates() const
{
    return 2;
}

int CoordgenInvertBondDOF::tier() const
{
    return INVERT_BOND_TIER;
}

void CoordgenInvertBondDOF::apply() const
{
    if (m_currentState != 0) {
        sketcherMinimizerPointF pivot = m_pivotAtom->getCoordinates();
        sketcherMinimizerPointF bondDir = m_boundAtom->getCoordinates() - pivot;
        sketcherMinimizerPointF normal(bondDir.y(), -bondDir.x());
        sketcherMinimizerPointF point1 = pivot + normal;
        sketcherMinimizerPointF point2 = pivot - normal;

        for (auto& atom : m_atoms) {
            atom->setCoordinates(sketcherMinimizerMaths::mirrorPoint(
                atom->getCoordinates(), point1, point2));
        }
    }
}

float CoordgenInvertBondDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return INVERTED_MACROCYCLE_BOND_PENALTY;
    }
    return 0.f;
}

CoordgenFlipRingDOF::CoordgenFlipRingDOF(
    sketcherMinimizerRing* ring,
    const std::vector<sketcherMinimizerAtom*>& fusionAtoms)
    : CoordgenFragmentDOF((*fusionAtoms.begin())->getFragment()),
      m_pivotAtom1(*(fusionAtoms.begin())),
      m_pivotAtom2(*(fusionAtoms.rbegin())),
      m_penalty(std::abs(int(ring->size() - 2 * fusionAtoms.size() + 2)))
{
    for (auto atom : ring->getAtoms()) {
        addAtom(atom);
    }
}

int CoordgenFlipRingDOF::numberOfStates() const
{
    return 2;
}

int CoordgenFlipRingDOF::tier() const
{
    return FLIP_RING_TIER;
}

void CoordgenFlipRingDOF::apply() const
{
    if (m_currentState != 0) {
        for (auto& atom : m_atoms) {
            atom->setCoordinates(sketcherMinimizerMaths::mirrorPoint(
                atom->getCoordinates(), m_pivotAtom1->getCoordinates(),
                m_pivotAtom2->getCoordinates()));
        }
    }
}

float CoordgenFlipRingDOF::getCurrentPenalty() const
{
    if (m_currentState != 0) {
        return FLIP_RING_PENALTY * m_penalty;
    }
    return 0.f;
}

sketcherMinimizerFragment::sketcherMinimizerFragment()
    : fixed(false), isTemplated(false), constrained(false), isChain(false),
      _bondToParent(nullptr), longestChainFromHere(0.f),
      numberOfChildrenAtoms(0), numberOfChildrenAtomsRank(0.f),
      m_parent(nullptr)
{
    m_dofs.push_back(new CoordgenFlipFragmentDOF(this));
    //    m_dofs.push_back(new CoordgenScaleFragmentDOF(this));
    m_dofs.push_back(new CoordgenChangeParentBondLengthFragmentDOF(this));
    m_dofs.push_back(new CoordgenRotateFragmentDOF(this));
}

sketcherMinimizerFragment::~sketcherMinimizerFragment()
{
    for (auto dof : m_dofs) {
        delete dof;
    }
}

void sketcherMinimizerFragment::addDof(CoordgenFragmentDOF* dof)
{
    m_dofs.push_back(dof);
}

const std::vector<CoordgenFragmentDOF*>& sketcherMinimizerFragment::getDofs()
{
    return m_dofs;
}

unsigned int sketcherMinimizerFragment::totalWeight() const
{
    int n = 0;
    for (auto m_atom : m_atoms) {
        n += m_atom->atomicNumber + m_atom->_implicitHs;
    }
    return n;
}
unsigned int sketcherMinimizerFragment::countDoubleBonds() const
{
    int n = 0;
    for (auto m_bond : m_bonds) {
        if (m_bond->bondOrder == 2) {
            ++n;
        }
    }
    return n;
}
unsigned int sketcherMinimizerFragment::countHeavyAtoms() const
{
    int n = 0;
    for (auto m_atom : m_atoms) {
        if (m_atom->atomicNumber != 6) {
            ++n;
        }
    }
    return n;
}

unsigned int sketcherMinimizerFragment::countConstrainedAtoms() const
{
    int n = 0;
    for (auto atom : m_atoms) {
        if (atom->constrained) {
            ++n;
        }
    }
    return n;
}

unsigned int sketcherMinimizerFragment::countFixedAtoms() const
{
    int n = 0;
    for (auto atom : m_atoms) {
        if (atom->fixed) {
            ++n;
        }
    }
    return n;
}

void sketcherMinimizerFragment::addAtom(sketcherMinimizerAtom* atom)
{
    m_atoms.push_back(atom);
    atom->setFragment(this);
}
void sketcherMinimizerFragment::addBond(sketcherMinimizerBond* bond)
{
    m_bonds.push_back(bond);
}
void sketcherMinimizerFragment::addRing(sketcherMinimizerRing* ring)
{
    m_rings.push_back(ring);
}

void sketcherMinimizerFragment::setAllCoordinatesToTemplate()
{
    for (sketcherMinimizerAtom* atom : m_atoms) {
        atom->setCoordinatesToTemplate();
    }
    if (_bondToParent) {
        _bondToParent->startAtom->setCoordinatesToTemplate();
        _bondToParent->endAtom->setCoordinatesToTemplate();
    }
    for (sketcherMinimizerFragment* child : _children) {
        child->_bondToParent->startAtom->setCoordinatesToTemplate();
        child->_bondToParent->endAtom->setCoordinatesToTemplate();
    }
}

void sketcherMinimizerFragment::storeCoordinateInformation()
{

    _coordinates.clear();
    sketcherMinimizerPointF origin(0.f, 0.f);
    float angle = 0.f;
    if (_bondToParent) {
        origin = _bondToParent->endAtom->getCoordinates();
        angle = atan2(_bondToParent->startAtom->coordinates.y() - origin.y(),
                      -_bondToParent->startAtom->coordinates.x() + origin.x());
    } else {
        if (!constrained && !fixed) {
            origin = m_atoms[0]->getCoordinates();
        }
    }

    float cosine = cos(-angle);
    float sine = sin(-angle);

    for (auto at : m_atoms) {
        sketcherMinimizerPointF c = at->coordinates - origin;
        c.rotate(sine, cosine);
        _coordinates[at] = c;
    }

    for (auto child : _children) {
        sketcherMinimizerAtom* at = child->_bondToParent->endAtom;
        sketcherMinimizerPointF c = at->coordinates - origin;
        c.rotate(sine, cosine);
        _coordinates[at] = c;
    }
}

void sketcherMinimizerFragment::setCoordinates(
    const sketcherMinimizerPointF& position, float angle)
{
    float sine = sin(angle);
    float cosine = cos(angle);
    assert(_coordinates.size() == m_atoms.size() + _children.size());
    for (auto& atom : _coordinates) {
        atom.first->setCoordinates(atom.second);
    }
    for (auto dof : m_dofs) {
        dof->apply();
    }
    for (auto& coords : _coordinates) {
        sketcherMinimizerAtom* atom = coords.first;
        sketcherMinimizerPointF initialCoordinates = atom->getCoordinates();
        initialCoordinates.rotate(sine, cosine);
        atom->setCoordinates(initialCoordinates + position);
    }
}
