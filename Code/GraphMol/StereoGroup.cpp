#include "StereoGroup.h"


namespace RDKit {

StereoGroup::StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms)
    : d_grouptype(grouptype), d_atoms(atoms) {}

StereoGroupType StereoGroup::getGroupType() const {return d_grouptype;}

const std::vector<Atom *>& StereoGroup::getAtoms() const {return d_atoms;}

void remove_groups_with_atom(const Atom* atom, std::vector<StereoGroup>& groups)
{
    auto contains_atom = [atom](const StereoGroup &group) {
      return std::find(group.getAtoms().cbegin(), group.getAtoms().cend(),
                       atom) != group.getAtoms().cend();
    };
    groups.erase(
        std::remove_if(groups.begin(), groups.end(), contains_atom),
        groups.end());
}

void remove_groups_with_atoms(const std::vector<Atom*>& atoms, std::vector<StereoGroup>& groups)
{
    auto contains_any_atom = [atoms](const StereoGroup &group) {
        for (auto atom: atoms) {
            if (std::find(group.getAtoms().cbegin(), group.getAtoms().cend(), atom) != group.getAtoms().cend()) {
                return true;
            }
        }
        return false;
    };
    groups.erase(
        std::remove_if(groups.begin(), groups.end(), contains_any_atom),
        groups.end());
}

}  // namespace RDKit
