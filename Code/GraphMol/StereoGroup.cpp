#include <algorithm>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "StereoGroup.h"
#include "Atom.h"

namespace RDKit {

namespace {
void checkAndStoreId(boost::dynamic_bitset<> &ids, StereoGroup &sg) {
  const auto groupId = sg.getId();
  if (groupId == 0) {
    return;
  } else if (groupId >= ids.size()) {
    ids.resize(groupId + 1);
  }
  if (ids[groupId]) {
    // This is likely to happen in reactions if two reactants have groups with
    // the same ids. In this case we just reset the id of one of them so that it
    // is reassigned later.
    sg.setId(0);
  } else {
    ids[groupId] = true;
  }
};

void updateId(const boost::dynamic_bitset<> &ids, unsigned &nextId,
              StereoGroup &sg) {
  if (sg.getId() == 0) {
    ++nextId;
    while (nextId < ids.size() && ids[nextId]) {
      ++nextId;
    }
    sg.setId(nextId);
  }
};
}  // namespace

StereoGroup::StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
                         unsigned id)
    : d_grouptype(grouptype), d_atoms(atoms), d_id{id} {}
StereoGroup::StereoGroup(StereoGroupType grouptype,
                         const std::vector<Atom *> &atoms, unsigned id)
    : d_grouptype(grouptype), d_atoms(std::move(atoms)), d_id{id} {}

StereoGroupType StereoGroup::getGroupType() const { return d_grouptype; }

const std::vector<Atom *> &StereoGroup::getAtoms() const { return d_atoms; }

void removeGroupsWithAtom(const Atom *atom, std::vector<StereoGroup> &groups) {
  auto containsAtom = [atom](const StereoGroup &group) {
    return std::find(group.getAtoms().cbegin(), group.getAtoms().cend(),
                     atom) != group.getAtoms().cend();
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAtom),
               groups.end());
}

void removeGroupsWithAtoms(const std::vector<Atom *> &atoms,
                           std::vector<StereoGroup> &groups) {
  auto containsAnyAtom = [&atoms](const StereoGroup &group) {
    for (auto atom : atoms) {
      if (std::find(group.getAtoms().cbegin(), group.getAtoms().cend(), atom) !=
          group.getAtoms().cend()) {
        return true;
      }
    }
    return false;
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAnyAtom),
               groups.end());
}

void assignStereoGroupIds(std::vector<StereoGroup> &groups) {
  if (groups.empty()) {
    return;
  }

  boost::dynamic_bitset<> andIds;
  boost::dynamic_bitset<> orIds;

  for (auto &sg : groups) {
    if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
      checkAndStoreId(andIds, sg);
    } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
      checkAndStoreId(orIds, sg);
    }
  }

  unsigned andId = 0;
  unsigned orId = 0;
  for (auto &sg : groups) {
    if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
      updateId(andIds, andId, sg);
    } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
      updateId(andIds, orId, sg);
    }
  }
}

}  // namespace RDKit

std::ostream &operator<<(std::ostream &target, const RDKit::StereoGroup &stg) {
  switch (stg.getGroupType()) {
    case RDKit::StereoGroupType::STEREO_ABSOLUTE:
      target << "ABS";
      break;
    case RDKit::StereoGroupType::STEREO_OR:
      target << "OR ";
      break;
    case RDKit::StereoGroupType::STEREO_AND:
      target << "AND";
      break;
  }
  target << " id: " << stg.getId();
  target << " atoms: { ";
  for (auto atom : stg.getAtoms()) {
    target << atom->getIdx() << ' ';
  }
  target << '}';

  return target;
}
