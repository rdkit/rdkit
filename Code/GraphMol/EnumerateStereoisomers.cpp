#include "EnumerateStereoisomers.h"

namespace RDKit {
    StereoEnumerationOptions::StereoEnumerationOptions(
        bool try_embedding, bool only_unassigned, bool only_stereo_groups,
        bool unique, unsigned int max_isomers, unsigned int rand) :
            try_embedding(try_embedding),
            only_unassigned(only_unassigned),
            only_stereo_groups(only_stereo_groups),
            unique(unique),
            max_isomers(max_isomers),
            rand(rand) {}

    _BondFlipper::_BondFlipper(Bond *bond) :
        bond(bond) {};

    void _BondFlipper::flip(bool flag) {
        if (flag) {
            bond->setStereo(Bond::STEREOCIS);
        } else {
            bond->setStereo(Bond::STEREOTRANS);
        }
    };

    _AtomFlipper::_AtomFlipper(Atom* atom) {
        atom = atom;
    };

    void _AtomFlipper::flip(bool flag) {
        if (flag) {
            atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        } else {
            atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        }
    };

    _StereoGroupFlipper::_StereoGroupFlipper(RDKit::StereoGroup* group) {
        _original_parities = std::vector<std::tuple<Atom*, RDKit::Atom::ChiralType> >();
        for (Atom* atom : group->getAtoms()) {
            _original_parities.push_back(std::make_tuple(atom, atom->getChiralTag()));
        }
    };

    void _StereoGroupFlipper::flip(bool flag) {
        if (flag) {
            for (auto& atom_parity : _original_parities) {
                std::get<0>(atom_parity)->setChiralTag(std::get<1>(atom_parity));
            } 
        } else {
            for (auto& atom_parity : _original_parities) {
                if (std::get<1>(atom_parity) == Atom::CHI_TETRAHEDRAL_CW) {
                    std::get<0>(atom_parity)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
                } else if (std::get<1>(atom_parity) == Atom::CHI_TETRAHEDRAL_CCW) {
                    std::get<0>(atom_parity)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
                }
            }
        }
    };

    std::vector<_Flipper*> _get_flippers(ROMol* mol, const StereoEnumerationOptions options) {
        std::vector<RDKit::Chirality::StereoInfo> potential_stereo = RDKit::Chirality::findPotentialStereo(*mol);

        std::vector<_Flipper*> flippers = std::vector<_Flipper*>();
        if (!options.only_stereo_groups) {
            for (auto atom : mol->atoms()) {
                if (atom->hasProp("_ChiralityPossible")) {
                    if (!options.only_unassigned || atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
                        flippers.push_back(new _AtomFlipper(atom));
                    }
                }
            }
            
            for (auto bond : mol->bonds()) {
                Bond::BondStereo bstereo = bond->getStereo();
                if (bstereo != Bond::STEREONONE) {
                    if (!options.only_unassigned || bstereo == Bond::STEREOANY) {
                        flippers.push_back(new _BondFlipper(bond));
                    }
                }
            }
        } 

        if (options.only_unassigned) {
            for (auto group : mol->getStereoGroups()) {
                if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
                    flippers.push_back(new _StereoGroupFlipper(&group));
                }
            }
        }
        return flippers;
    };
    
    unsigned int get_stereoisomer_count(ROMol* mol, const StereoEnumerationOptions options=StereoEnumerationOptions()) {
        std::vector<_Flipper*> flippers = _get_flippers(mol, options);
        return std::pow(2, flippers.size());
    }
}