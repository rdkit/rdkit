//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
#include <GraphMol/FileParsers/MolToMacroMol.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RDKit {
namespace {

constexpr const char *MACROMOL_GROUP = "MACROMOL";

struct AttachmentPointSpec {
  unsigned int sourceAtomIdx;
  int attachmentIdx;
};

struct GroupSpec {
  const MacroMolTemplate *templ = nullptr;
  std::vector<unsigned int> sourceAtomIndices;
  std::vector<AttachmentPointSpec> attachmentPoints;

  std::optional<int> getAttachmentIdx(unsigned int atomIdx) const {
    const auto it = std::find_if(
        attachmentPoints.begin(), attachmentPoints.end(),
        [atomIdx](const AttachmentPointSpec &spec) {
          return spec.sourceAtomIdx == atomIdx;
        });
    if (it == attachmentPoints.end()) {
      return std::nullopt;
    }
    return it->attachmentIdx;
  }
};

struct MainGroupQuery {
  std::unique_ptr<RWMol> query;
  std::vector<unsigned int> queryToTemplateAtom;
};

MainGroupQuery makeMainGroupQuery(const MacroMolTemplate &templ) {
  std::vector<unsigned int> mainAtoms(templ.getMainAtomIdxs().begin(),
                                      templ.getMainAtomIdxs().end());
  std::sort(mainAtoms.begin(), mainAtoms.end());
  const std::unordered_set<unsigned int> mainAtomSet(mainAtoms.begin(),
                                                     mainAtoms.end());

  auto query = std::make_unique<RWMol>(templ.getMol());
  query->beginBatchEdit();
  for (auto atom : query->atoms()) {
    if (!mainAtomSet.count(atom->getIdx())) {
      query->removeAtom(atom);
    }
  }
  query->commitBatchEdit();
  query->updatePropertyCache(false);
  return {std::move(query), std::move(mainAtoms)};
}

std::vector<MatchVectType> findTemplateMatches(const ROMol &mol,
                                               const RWMol &query) {
  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.useChirality = true;
  params.useQueryQueryMatches = false;
  params.uniquify = false;
  params.maxMatches = 0;
  auto matches = SubstructMatch(mol, query, params);

  auto matchKey = [](const MatchVectType &match) {
    std::vector<int> sourceAtoms;
    std::vector<std::pair<int, int>> queryToSource;
    for (const auto &pair : match) {
      sourceAtoms.push_back(pair.second);
      queryToSource.push_back(pair);
    }
    std::sort(sourceAtoms.begin(), sourceAtoms.end());
    std::sort(queryToSource.begin(), queryToSource.end());
    return std::make_pair(sourceAtoms, queryToSource);
  };
  std::sort(matches.begin(), matches.end(),
            [&matchKey](const auto &lhs, const auto &rhs) {
              return matchKey(lhs) < matchKey(rhs);
            });
  return matches;
}

bool violatesAttachmentInvariant(const ROMol &mol, const GroupSpec &group) {
  auto atoms = group.sourceAtomIndices;
  std::sort(atoms.begin(), atoms.end());
  for (const auto atomIdx : atoms) {
    unsigned int externalBondCount = 0;
    for (const auto *bond : mol.atomBonds(mol.getAtomWithIdx(atomIdx))) {
      if (!std::binary_search(atoms.begin(), atoms.end(),
                              bond->getOtherAtomIdx(atomIdx))) {
        ++externalBondCount;
      }
    }
    const unsigned int allowed =
        group.getAttachmentIdx(atomIdx).has_value() ? 1u : 0u;
    if (externalBondCount > allowed) {
      return true;
    }
  }
  return false;
}

std::optional<GroupSpec> makeGroupForMatch(
    const MatchVectType &match, const ROMol &mol,
    const MacroMolTemplate &templ, const MainGroupQuery &mainQuery,
    const std::vector<bool> &claimed) {
  GroupSpec group;
  group.templ = &templ;
  std::unordered_map<unsigned int, unsigned int> templateToSource;
  for (const auto &[queryIdx, sourceIdxInt] : match) {
    const auto sourceIdx = static_cast<unsigned int>(sourceIdxInt);
    if (claimed[sourceIdx]) {
      return std::nullopt;
    }
    group.sourceAtomIndices.push_back(sourceIdx);
    templateToSource[mainQuery.queryToTemplateAtom[queryIdx]] = sourceIdx;
  }

  std::unordered_set<unsigned int> seenAttachAtoms;
  for (const auto &leavingGroup : templ.getLeavingGroups()) {
    const auto it = templateToSource.find(leavingGroup.attachAtomIdx);
    if (it == templateToSource.end()) {
      return std::nullopt;
    }
    if (!seenAttachAtoms.insert(it->second).second) {
      throw ValueErrorException(
          "MolToMacroMol does not support multiple attachment points on one "
          "template atom");
    }
    group.attachmentPoints.push_back({it->second, leavingGroup.attachPoint});
  }
  if (violatesAttachmentInvariant(mol, group)) {
    return std::nullopt;
  }
  std::sort(group.sourceAtomIndices.begin(), group.sourceAtomIndices.end());
  return group;
}

std::vector<GroupSpec> findGroups(const ROMol &mol,
                                  const MacroMolTemplateLibrary &templates) {
  std::vector<GroupSpec> groups;
  std::vector<bool> claimed(mol.getNumAtoms(), false);
  auto orderedTemplates = templates.getTemplates();
  std::stable_sort(
      orderedTemplates.begin(), orderedTemplates.end(),
      [](const auto *lhs, const auto *rhs) {
        return lhs->getMainAtomIdxs().size() > rhs->getMainAtomIdxs().size();
      });
  for (const auto *templ : orderedTemplates) {
    const auto mainQuery = makeMainGroupQuery(*templ);
    for (const auto &match : findTemplateMatches(mol, *mainQuery.query)) {
      auto group = makeGroupForMatch(match, mol, *templ, mainQuery, claimed);
      if (!group) {
        continue;
      }
      for (const auto atomIdx : group->sourceAtomIndices) {
        claimed[atomIdx] = true;
      }
      groups.push_back(std::move(*group));
    }
  }
  return groups;
}

void addGroupMetadata(RWMol &result, const GroupSpec &group) {
  SubstanceGroup sgroup(&result, "SUP");
  sgroup.setProp<bool>(MACROMOL_GROUP, true);
  sgroup.setProp<std::string>("CLASS",
                              monomerClassToString(group.templ->getMonomerClass()));
  sgroup.setProp<std::string>("LABEL", group.templ->getSymbol());
  sgroup.setAtoms(group.sourceAtomIndices);
  for (const auto &attach : group.attachmentPoints) {
    sgroup.addAttachPoint(attach.sourceAtomIdx, -1,
                          std::to_string(attach.attachmentIdx));
  }
  addSubstanceGroup(result, std::move(sgroup));
}

std::unique_ptr<RWMol> annotate(const ROMol &mol,
                                const std::vector<GroupSpec> &groups) {
  auto result = std::make_unique<RWMol>(mol);
  for (const auto &group : groups) {
    addGroupMetadata(*result, group);
  }
  return result;
}

struct CollapseGroup {
  std::vector<unsigned int> atoms;
  std::unordered_map<unsigned int, int> attachmentByAtom;
  const MacroMolTemplate *templ = nullptr;
  unsigned int outputAtom = std::numeric_limits<unsigned int>::max();
};

const MacroMolTemplate *findTemplate(const MacroMolTemplateLibrary &templates,
                                     const SubstanceGroup &sgroup) {
  const auto monomerClass = monomerClassFromString(
      sgroup.getProp<std::string>("CLASS"));
  const auto label = sgroup.getProp<std::string>("LABEL");
  auto *templ = templates.getBySymbol(monomerClass, label);
  if (!templ) {
    templ = templates.getByName(monomerClass, label);
  }
  return templ;
}

void copyRegularBond(MacroMol &result, const Bond &source,
                     unsigned int beginIdx, unsigned int endIdx,
                     const std::vector<unsigned int> &outputForAtom) {
  std::unique_ptr<Bond> bond(source.copy());
  bond->setBeginAtomIdx(beginIdx);
  bond->setEndAtomIdx(endIdx);
  for (auto &stereoAtom : bond->getStereoAtoms()) {
    if (stereoAtom >= 0 &&
        static_cast<unsigned int>(stereoAtom) < outputForAtom.size()) {
      stereoAtom = static_cast<int>(
          outputForAtom[static_cast<unsigned int>(stereoAtom)]);
    }
  }
  result.addBond(bond.get(), true);
  bond.release();
}

}  // namespace

std::unique_ptr<RWMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates) {
  return annotate(mol, findGroups(mol, templates));
}

std::unique_ptr<MacroMol> CollapseMacroMol(
    const ROMol &groupedMol, const MacroMolTemplateLibrary &templates) {
  auto result = std::make_unique<MacroMol>();
  std::vector<int> groupForAtom(groupedMol.getNumAtoms(), -1);
  std::vector<CollapseGroup> groups;

  for (const auto &sgroup : getSubstanceGroups(groupedMol)) {
    bool isMacroGroup = false;
    if (!sgroup.getPropIfPresent<bool>(MACROMOL_GROUP, isMacroGroup) ||
        !isMacroGroup) {
      continue;
    }
    CollapseGroup group;
    group.atoms = sgroup.getAtoms();
    group.templ = findTemplate(templates, sgroup);
    if (!group.templ || group.atoms.empty()) {
      throw ValueErrorException("invalid MacroMol SubstanceGroup");
    }
    for (const auto atomIdx : group.atoms) {
      if (atomIdx >= groupedMol.getNumAtoms() || groupForAtom[atomIdx] != -1) {
        throw ValueErrorException("overlapping or invalid MacroMol groups");
      }
      groupForAtom[atomIdx] = static_cast<int>(groups.size());
    }
    for (const auto &attach : sgroup.getAttachPoints()) {
      if (attach.aIdx >= groupedMol.getNumAtoms() ||
          groupForAtom[attach.aIdx] != static_cast<int>(groups.size()) ||
          attach.lvIdx != -1) {
        throw ValueErrorException("invalid MacroMol attachment point");
      }
      int attachmentIdx = 0;
      try {
        attachmentIdx = std::stoi(attach.id);
      } catch (const std::exception &) {
        throw ValueErrorException("invalid MacroMol attachment point ID");
      }
      group.attachmentByAtom.emplace(attach.aIdx, attachmentIdx);
    }
    groups.push_back(std::move(group));
  }

  std::vector<unsigned int> outputForAtom(groupedMol.getNumAtoms());
  for (const auto *atom : groupedMol.atoms()) {
    const auto atomIdx = atom->getIdx();
    const auto groupIdx = groupForAtom[atomIdx];
    if (groupIdx >= 0) {
      auto &group = groups[static_cast<unsigned int>(groupIdx)];
      if (group.outputAtom == std::numeric_limits<unsigned int>::max()) {
        group.outputAtom = result->addMacroAtom(group.templ->getSymbol(),
                                                group.templ->getMonomerClass());
      }
      outputForAtom[atomIdx] = group.outputAtom;
    } else {
      auto copy = std::unique_ptr<Atom>(atom->copy());
      outputForAtom[atomIdx] = result->addAtom(copy.get(), true, true);
      copy.release();
    }
  }

  for (const auto *bond : groupedMol.bonds()) {
    const auto begin = bond->getBeginAtomIdx();
    const auto end = bond->getEndAtomIdx();
    const auto beginGroup = groupForAtom[begin];
    const auto endGroup = groupForAtom[end];
    if (beginGroup >= 0 && beginGroup == endGroup) {
      continue;
    }

    const auto beginAttach = beginGroup >= 0
                                 ? groups[static_cast<unsigned int>(beginGroup)]
                                       .attachmentByAtom.find(begin)
                                 : decltype(groups.front().attachmentByAtom)::const_iterator{};
    const auto endAttach = endGroup >= 0
                               ? groups[static_cast<unsigned int>(endGroup)]
                                     .attachmentByAtom.find(end)
                               : decltype(groups.front().attachmentByAtom)::const_iterator{};
    const bool hasBeginAttach =
        beginGroup >= 0 && beginAttach !=
                              groups[static_cast<unsigned int>(beginGroup)]
                                  .attachmentByAtom.end();
    const bool hasEndAttach =
        endGroup >= 0 && endAttach !=
                            groups[static_cast<unsigned int>(endGroup)]
                                .attachmentByAtom.end();

    if (beginGroup >= 0 || endGroup >= 0) {
      if ((beginGroup >= 0 && !hasBeginAttach) ||
          (endGroup >= 0 && !hasEndAttach)) {
        throw ValueErrorException(
            "MacroMol group has an unannotated crossing bond");
      }
      const auto bondType = bond->getBondType();
      if (beginGroup >= 0 && endGroup >= 0) {
        result->addMacroBond(outputForAtom[begin], outputForAtom[end],
                             beginAttach->second, endAttach->second,
                             bondType);
      } else if (beginGroup >= 0) {
        result->addMacroAtomToAtomBond(outputForAtom[begin], outputForAtom[end],
                                       beginAttach->second, bondType);
      } else {
        result->addAtomToMacroAtomBond(outputForAtom[begin], outputForAtom[end],
                                       endAttach->second, bondType);
      }
    } else {
      copyRegularBond(*result, *bond, outputForAtom[begin],
                      outputForAtom[end], outputForAtom);
    }
  }
  return result;
}

}  // namespace RDKit
