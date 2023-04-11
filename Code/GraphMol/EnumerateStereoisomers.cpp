#include "EnumerateStereoisomers.h"

namespace RDKit {
    StereoEnumerationOptions::StereoEnumerationOptions(
        bool try_embedding, bool only_unassigned, bool only_stereo_groups,
        bool unique, unsigned int max_isomers, unsigned int rand) {
            try_embedding = try_embedding;
            only_unassigned = only_unassigned;
            only_stereo_groups = only_stereo_groups;
            unique = unique;
            max_isomers = max_isomers;
            rand = rand;
    };

    _BondFlipper::_BondFlipper(Bond *bond) {
        bond = bond;
    };

    void _BondFlipper::flip(bool flag) {
        if (flag) {
            bond->setStereo(Bond::STEREOCIS);
        } else {
            bond->setStereo(Bond::STEREOTRANS);
        }
    }

    _AtomFlipper::_AtomFlipper(Atom* atom) {
        atom = atom;
    }

    void _AtomFlipper::flip(bool flag) {
        if (flag) {
            atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        } else {
            atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        }
    }
}