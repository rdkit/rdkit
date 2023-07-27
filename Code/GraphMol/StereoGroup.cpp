#include <algorithm>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "StereoGroup.h"
#include "Atom.h"
#include "ROMol.h"

namespace RDKit {

namespace {
void storeIdsInUse(boost::dynamic_bitset<> &ids, StereoGroup &sg) {
  const auto groupId = sg.getWriteId();
  if (groupId == 0) {
    return;
  } else if (groupId >= ids.size()) {
    ids.resize(groupId + 1);
  }
  if (ids[groupId]) {
    // This id is duplicate, let's reset it so we can reassign it later
    BOOST_LOG(rdWarningLog)
        << "StereoGroup ID " << groupId
        << " is used by more than one group, and will be reassined"
        << std::endl;
    sg.setWriteId(0);
  } else {
    ids[groupId] = true;
  }
};

void assignMissingIds(const boost::dynamic_bitset<> &ids, unsigned &nextId,
                      StereoGroup &sg) {
  if (sg.getWriteId() == 0) {
    ++nextId;
    while (nextId < ids.size() && ids[nextId]) {
      ++nextId;
    }
    sg.setWriteId(nextId);
  }
};
}  // namespace

StereoGroup::StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
                         unsigned readId)
    : d_grouptype(grouptype), d_atoms(atoms), d_readId{readId} {}
StereoGroup::StereoGroup(StereoGroupType grouptype,
                         const std::vector<Atom *> &atoms, unsigned readId)
    : d_grouptype(grouptype), d_atoms(std::move(atoms)), d_readId{readId} {}

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
      storeIdsInUse(andIds, sg);
    } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
      storeIdsInUse(orIds, sg);
    }
  }

  unsigned andId = 0;
  unsigned orId = 0;
  for (auto &sg : groups) {
    if (sg.getGroupType() == StereoGroupType::STEREO_AND) {
      assignMissingIds(andIds, andId, sg);
    } else if (sg.getGroupType() == StereoGroupType::STEREO_OR) {
      assignMissingIds(andIds, orId, sg);
    }
  }
}

void forwardStereoGroupIds(ROMol &mol) {
  auto stgs = mol.getStereoGroups();
  for (auto &stg : stgs) {
    stg.setWriteId(stg.getReadId());
  }
  mol.setStereoGroups(stgs);
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
  target << " rId: " << stg.getReadId();
  target << " wRd: " << stg.getWriteId();
  target << " atoms: { ";
  for (auto atom : stg.getAtoms()) {
    target << atom->getIdx() << ' ';
  }
  target << '}';

  return target;
}
