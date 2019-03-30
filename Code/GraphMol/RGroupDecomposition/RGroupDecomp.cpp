//  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "RGroupDecomp.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FMCS/FMCS.h>
#include <boost/scoped_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <set>
#include <utility>
#include <vector>

//#define DEBUG

namespace RDKit {

// Attachment Points
//  labeled cores => isotopes
//  atom mappings
//  atom indices => use -1 - atom index, range is [-1, ...., -num_atoms]
namespace {
const std::string RLABEL = "tempRlabel";
const std::string SIDECHAIN_RLABELS = "sideChainRlabels";
const std::string done = "RLABEL_PROCESSED";

bool setLabel(Atom *atom, int label, std::set<int> &labels, int &maxLabel,
              bool relabel, const std::string &type) {
  if (type == "IsotopeLabels") {
    atom->setIsotope(0);
  }

  if (label) {
    if (labels.find(label) != labels.end()) {
      if (relabel)
        label = maxLabel + 1;
      else
        // XXX FIX me - get label id
        throw ValueErrorException(
            std::string("Duplicate label in input, current type is:") + type);
    }

    atom->setProp<int>(RLABEL, label);
    labels.insert(label);
    maxLabel = label + 1;
    return true;
  }
  return false;
}

bool hasDummy(const RWMol &core) {
  for (RWMol::ConstAtomIterator atIt = core.beginAtoms();
       atIt != core.endAtoms(); ++atIt) {
    if ((*atIt)->getAtomicNum() == 0) return true;
  }
  return false;
}
}  // namespace

bool RGroupDecompositionParameters::prepareCore(RWMol &core,
                                                const RWMol *alignCore) {
  const bool relabel = labels & RelabelDuplicateLabels;
  if (alignCore && (alignment & MCS)) {
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ROMOL_SPTR(new ROMol(core)));
    mols.push_back(ROMOL_SPTR(new ROMol(*alignCore)));
    MCSResult res = findMCS(mols);
    if (res.isCompleted()) {
      RWMol *m = SmartsToMol(res.SmartsString);
      if (m) {
        MatchVectType match1;
        MatchVectType match2;

        bool target_matched1 = SubstructMatch(core, *m, match1);
        bool target_matched2 = SubstructMatch(*alignCore, *m, match2);
        CHECK_INVARIANT(match1.size() == match2.size(),
                        "Matches should be the same size in prepareCore");

        if (target_matched1 && target_matched2) {
          for (size_t i = 0; i < match1.size(); ++i) {
            int queryAtomIdx1 = match1[i].first;
            int coreAtomIdx = match1[i].second;
            int queryAtomIdx2 = match2[i].first;
            int alignCoreAtomIdx = match2[i].second;
            CHECK_INVARIANT(queryAtomIdx1 == queryAtomIdx2,
                            "query atoms aren't the same");
            const Atom *coreAtm = core.getAtomWithIdx(coreAtomIdx);
            const Atom *alignCoreAtm =
                alignCore->getAtomWithIdx(alignCoreAtomIdx);
            int rlabel = alignCoreAtm->getProp<int>(RLABEL);
            coreAtm->setProp(RLABEL, rlabel);
          }
        }
        delete m;
      }
    }
  }

  std::set<int> foundLabels;

  int maxLabel = 0;
  int nextOffset = 0;
  std::map<int, int> atomToLabel;

  for (RWMol::AtomIterator atIt = core.beginAtoms(); atIt != core.endAtoms();
       ++atIt) {
    Atom *atom = *atIt;
    bool found = false;

    if (atom->hasProp(RLABEL)) found = true;

    if (!found && (labels & IsotopeLabels)) {
      if (setLabel(atom, rdcast<int>(atom->getIsotope()), foundLabels, maxLabel,
                   relabel, "IsotopeLabels"))
        found = true;
    }

    if (!found && (labels & AtomMapLabels)) {
      if (setLabel(atom, rdcast<int>(atom->getAtomMapNum()), foundLabels,
                   maxLabel, relabel, "AtomMapLabels"))
        found = true;
    }

    if (!found && (labels & AtomIndexLabels)) {
      if (setLabel(atom, indexOffset - atom->getIdx(), foundLabels, maxLabel,
                   relabel, "IndexLabels"))
        nextOffset++;
      found = true;
    }

    int rlabel;
    if (atom->getPropIfPresent(RLABEL, rlabel)) {
      atomToLabel[atom->getIdx()] = rlabel;
    }
  }

  indexOffset -= nextOffset;

  MolOps::AdjustQueryParameters adjustParams;
  adjustParams.makeDummiesQueries = true;
  adjustParams.adjustDegree = false;
  adjustQueryProperties(core, &adjustParams);

  for (auto &it : atomToLabel)
    core.getAtomWithIdx(it.first)->setProp(RLABEL, it.second);
  return true;
}

namespace {
// RGroup Class to hold the attached bits

struct RGroupData {
  boost::shared_ptr<RWMol> combinedMol;
  std::vector<boost::shared_ptr<ROMol>> mols;  // All the mols in the rgroup
  std::set<std::string> smilesSet;             // used for rgroup equivalence
  std::string
      smiles;  // smiles for all the mols in the rgroup (with attachments)
  std::set<int> attachments;  // attachment points
  bool labelled;

 private:
  RGroupData(const RGroupData &rhs);

 public:
  RGroupData()
      : combinedMol(),
        mols(),
        smilesSet(),
        smiles(),
        attachments(),
        labelled(false) {}

  void add(boost::shared_ptr<ROMol> newMol,
           const std::vector<int> &rlabel_attachments) {
    // some fragments can be add multiple times if they are cyclic
    for (auto &mol : mols) {
      if (newMol.get() == mol.get()) return;
    }

    labelled = false;
    std::copy(rlabel_attachments.begin(), rlabel_attachments.end(),
              std::inserter(attachments, attachments.end()));

    mols.push_back(newMol);
    std::string smi = MolToSmiles(*newMol, true);
    // REVIEW: we probably shouldn't be using a set here... the merging of
    // duplicates is likely not what we want
    smilesSet.insert(smi);
    if (!combinedMol.get()) {
      combinedMol = boost::shared_ptr<RWMol>(new RWMol(*mols[0].get()));
    } else {
      ROMol *m = combineMols(*combinedMol.get(), *newMol.get());
      m->updateProps(*combinedMol.get());
      combinedMol.reset(new RWMol(*m));
      delete m;
    }
    smiles = getSmiles();
    combinedMol->setProp(common_properties::internalRgroupSmiles, smiles);
  }

  std::map<int, int> getNumBondsToRlabels() const {
    std::map<int, int> rlabelsUsedCount;

    for (ROMol::AtomIterator atIt = combinedMol->beginAtoms();
         atIt != combinedMol->endAtoms(); ++atIt) {
      Atom *atom = *atIt;
      int rlabel;
      if (atom->getPropIfPresent<int>(RLABEL, rlabel))
        rlabelsUsedCount[rlabel] += 1;
    }
    return rlabelsUsedCount;
  }

  bool isHydrogen() const {  // is the rgroup all Hs
    for (const auto &mol : mols) {
      for (ROMol::AtomIterator atIt = mol->beginAtoms();
           atIt != mol->endAtoms(); ++atIt) {
        if ((*atIt)->getAtomicNum() > 1) return false;
      }
    }
    return true;
  }

 private:
  std::string getSmiles()
      const {  // compute the canonical smiles for the attachments
    std::string s;
    for (const auto &it : smilesSet) {
      if (s.length()) s += ".";
      s += it;
    }
    return s;
  }
};
}  // namespace

namespace {
typedef boost::shared_ptr<RGroupData> RData;

typedef std::map<int, RData> R_DECOMP;
struct RGroupMatch {
  // RGroupMatch is the decomposition for a single molecule
  size_t core_idx;   // index of the matching core
  R_DECOMP rgroups;  // rlabel->RGroupData mapping

  RGroupMatch(size_t core_index, R_DECOMP input_rgroups)
      : core_idx(core_index), rgroups(std::move(input_rgroups)) {}
};

void ADD_MATCH(R_DECOMP &match, int rlabel) {
  if (match.find(rlabel) == match.end())
    match[rlabel] = boost::make_shared<RGroupData>();
}

struct CartesianProduct {
  std::vector<size_t> permutation;
  std::vector<size_t> sizes;
  size_t maxPermutations;
  size_t permutationCount;
  CartesianProduct(const std::vector<size_t> &inputSizes)
      : permutation(inputSizes.size(), 0),
        sizes(inputSizes),
        permutationCount(0) {
    maxPermutations = 1;
    for (unsigned long size : sizes)
      maxPermutations *= size;  // may overflow....
  }

  bool next() {
    ++permutationCount;
    if (permutationCount == 1) {
      return true;
    }

    return increment(0);
  }

  bool increment(size_t rowToIncrement) {
    if (permutationCount > maxPermutations) return false;

    permutation[rowToIncrement] += 1;
    size_t max_index_of_row = sizes[rowToIncrement] - 1;
    if (permutation[rowToIncrement] > max_index_of_row) {
      permutation[rowToIncrement] = 0;
      return increment(rowToIncrement + 1);
    }
    return true;
  }
};

// stupid total score
double score(const std::vector<size_t> &permutation,
             const std::vector<std::vector<RGroupMatch>> &matches,
             const std::set<int> &labels) {
  double score = 1.;

#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Scoring permutation "
            << " num matches: " << matches.size() << std::endl;
#endif

  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    std::map<std::string, int> matchSet;
    std::map<std::set<int>, int> linkerMatchSet;
    std::map<std::string, int> onlyH;

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg != matches[m][permutation[m]].rgroups.end()) {
#ifdef DEBUG
        std::cerr << "  combined: " << MolToSmiles(*rg->second->combinedMol)
                  << std::endl;
        std::cerr << " RGroup: " << rg->second->smiles << " "
                  << rg->second->smiles.find_first_not_of("0123456789[]*H:.");
#endif
        matchSet[rg->second->smiles] += 1;
        // detect whether or not this is an H
        if (rg->second->smiles.find_first_not_of("0123456789[]*H:.") >=
            rg->second->smiles.length()) {
          onlyH[rg->second->smiles] = 0;
        } else {
          onlyH[rg->second->smiles] = 1;
        }
#ifdef DEBUG
        std::cerr << " " << rg->second->combinedMol->getNumAtoms(false)
                  << " isH: " << onlyH[rg->second->smiles]
                  << " score: " << matchSet[rg->second->smiles] << std::endl;
#endif
        // XXX Use fragment counts to see if we are linking cycles?
        if (rg->second->smiles.find(".") == std::string::npos &&
            rg->second->attachments.size() > 1) {
          linkerMatchSet[rg->second->attachments]++;
#ifdef DEBUG
          std::cerr << " Linker Score: "
                    << linkerMatchSet[rg->second->attachments]++ << std::endl;
#endif
        }
      }
    }

    // get the counts for each rgroup found and sort in reverse order
    std::vector<float> equivalentRGroupCount;

    for (std::map<std::string, int>::const_iterator it = matchSet.begin();
         it != matchSet.end(); ++it) {
#ifdef DEBUG
      std::cerr << " equiv: " << it->first << " " << it->second << " "
                << permutation.size() << std::endl;
#endif

      // if the rgroup is hydrogens, only consider if the group is all
      //  hydrogen, otherwise score based on the non hydrogens
      if (onlyH[it->first]) {
        if (static_cast<size_t>(it->second) == permutation.size()) {
          equivalentRGroupCount.push_back(static_cast<float>(it->second));
        } else {
          // hydrogens in a mixed group don't contribute to the score
          equivalentRGroupCount.push_back(0.0);
        }
      } else {
        equivalentRGroupCount.push_back(static_cast<float>(it->second));
      }
    }
    std::sort(equivalentRGroupCount.begin(), equivalentRGroupCount.end(),
              std::greater<float>());

    double tempScore = 1.;
    // score the sets from the largest to the smallest
    //  each smaller set gets penalized (i+1) below
    //  1.0 is the perfect score
    for (size_t i = 0; i < equivalentRGroupCount.size(); ++i) {
      auto lscore =
          equivalentRGroupCount[i] / ((i + 1) * (double)matches.size());
      if (lscore > 0) tempScore *= lscore * lscore;
#ifdef DEBUG
      std::cerr << "    lscore^2 " << i << ": " << lscore * lscore << std::endl;
#endif
    }

    // overweight linkers with the same attachments points....
    //  because these belong to 2 rgroups we really want these to stay
    //  ** this heuristic really should be taken care of above **
    int maxLinkerMatches = 0;
    for (std::map<std::set<int>, int>::const_iterator it =
             linkerMatchSet.begin();
         it != linkerMatchSet.end(); ++it) {
      if (it->second > 1) {
        if (it->second > maxLinkerMatches) maxLinkerMatches = it->second;
      }
    }
#ifdef DEBUG
    std::cerr << "Max Linker Matches :" << maxLinkerMatches << std::endl;
#endif
    double increment = 1.0;        // no change in score
    double linkerIncrement = 1.0;  // no change in score
    if (maxLinkerMatches) {
      linkerIncrement = (double)(maxLinkerMatches) / (double)matches.size();
    } else {
      increment = tempScore;
    }

    score += increment * linkerIncrement;
#ifdef DEBUG
    std::cerr << "Increment: " << increment
              << " Linker_Increment: " << linkerIncrement << std::endl;
    std::cerr << "increment*linkerIncrement: " << increment * linkerIncrement
              << std::endl;
    std::cerr << "Score = " << score << std::endl;
#endif
  }

  return score;
}
}  // namespace

const unsigned int EMPTY_CORE_LABEL = -100000;

struct RGroupDecompData {
  // matches[mol_idx] == vector of potential matches
  std::map<int, RWMol> cores;
  std::map<std::string, int> newCores;  // new "cores" found along the way
  int newCoreLabel;
  RGroupDecompositionParameters params;

  std::vector<std::vector<RGroupMatch>> matches;
  std::set<int> labels;
  std::vector<size_t> permutation;
  std::map<int, std::vector<int>> userLabels;

  std::vector<int> processedRlabels;
  std::map<int, boost::shared_ptr<RWMol>> labelledCores;

  std::map<int, int> finalRlabelMapping;

  RGroupDecompData(const RWMol &inputCore,
                   RGroupDecompositionParameters inputParams)
      : cores(),
        newCores(),
        newCoreLabel(EMPTY_CORE_LABEL),
        params(std::move(inputParams)) {
    cores[0] = inputCore;
    prepareCores();
  }

  RGroupDecompData(const std::vector<ROMOL_SPTR> &inputCores,
                   RGroupDecompositionParameters inputParams)
      : cores(),
        newCores(),
        newCoreLabel(EMPTY_CORE_LABEL),
        params(std::move(inputParams)) {
    for (size_t i = 0; i < inputCores.size(); ++i) {
      cores[i] = *inputCores[i].get();
    }

    prepareCores();
  }

  void prepareCores() {
    size_t idx = 0;
    for (auto coreIt = cores.begin(); coreIt != cores.end(); ++coreIt, ++idx) {
      RWMol *alignCore = coreIt->first ? &cores[0] : nullptr;
      params.prepareCore(coreIt->second, alignCore);
      labelledCores[coreIt->first] =
          boost::shared_ptr<RWMol>(new RWMol(coreIt->second));
    }
  }

  void setRlabel(Atom *atom, int rlabel) {
    // XXX Fix me - use parameters to decide what to do.  Currenty does
    // everything
    if (params.rgroupLabelling & AtomMap) atom->setAtomMapNum(rlabel);

    if (params.rgroupLabelling & MDLRGroup) {
      std::string dLabel = "R" + std::to_string(rlabel);
      atom->setProp(common_properties::dummyLabel, dLabel);
      setAtomRLabel(atom, rlabel);
    }

    if (params.rgroupLabelling & Isotope) atom->setIsotope(rlabel);
  }

  void prune() {  // prune all but the current "best" permutation of matches
    for (size_t mol_idx = 0; mol_idx < permutation.size(); ++mol_idx) {
      std::vector<RGroupMatch> keepVector;
      keepVector.push_back(matches[mol_idx][permutation[mol_idx]]);
      matches[mol_idx] = keepVector;
    }
    permutation = std::vector<size_t>(matches.size(), 0);
  }

  // Return the RGroups with the current "best" permutation
  //  of matches.
  std::vector<RGroupMatch> GetCurrentBestPermutation() const {
    const bool removeAllHydrogenRGroups = params.removeAllHydrogenRGroups;

    std::vector<RGroupMatch> result;  // std::map<int, RGroup> > result;
    for (size_t i = 0; i < permutation.size(); ++i) {
      PRECONDITION(i < matches.size(), "Best Permutation mol idx out of range");
      PRECONDITION(permutation[i] < matches[i].size(),
                   "Selected match at permutation out of range");
      result.push_back(matches[i][permutation[i]]);
    }

    if (removeAllHydrogenRGroups) {
      // if a label is all hydrogens, remove it
      for (int label : labels) {
        bool allH = true;
        for (auto &i : result) {
          R_DECOMP::const_iterator rgroup = i.rgroups.find(label);
          if (rgroup == i.rgroups.end() || !rgroup->second->isHydrogen()) {
            allH = false;
            break;
          }
        }

        if (allH) {
          for (auto &i : result) {
            i.rgroups.erase(label);
          }
        }
      }
    }
    return result;
  }

  void relabelCore(RWMol &mol, std::map<int, int> &mappings,
                   const std::set<int> &userLabels,
                   const std::set<int> &indexLabels,
                   std::map<int, std::vector<int>> extraAtomRLabels) {
    // Now remap to proper rlabel ids
    //  if labels are positive, they come from User labels
    //  if they are negative, they come from indices and should be
    //  numbered *after* the user labels.
    //
    //  Some indices are attached to multiple bonds,
    //   these rlabels should be incrementally added last
    int count = 0;
    std::map<int, Atom *> atoms;

    // a core only has one labelled index
    //  a secondary structure extraAtomRLabels contains the number
    //  of bonds between this atom and the side chain

    // a sidechain atom has a vector of the attachments back to the
    //  core that takes the place of numBondsToRlabel

    std::map<int, std::vector<int>> bondsToCore;

    for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      Atom *atom = *atIt;

      if (atom->hasProp(RLABEL)) {
        int rlabel = (*atIt)->getProp<int>(RLABEL);  // user label
        PRECONDITION(atoms.find(rlabel) == atoms.end(),
                     "Duplicate labels in rgroup core!");
        atoms[rlabel] = *atIt;
      }
    }

    std::vector<std::pair<Atom *, Atom *>> atomsToAdd;  // adds -R if necessary

    // Deal with user supplied labels
    for (int userLabel : userLabels) {
      auto atm = atoms.find(userLabel);
      if (atm == atoms.end()) continue;  // label not used in the rgroup
      Atom *atom = atm->second;
      mappings[userLabel] = userLabel;
      if (count < userLabel) count = userLabel;
      if (atom->getAtomicNum() == 0) {  // add to existing dummy/rlabel
        setRlabel(atom, userLabel);
      } else {  // adds new rlabel
        auto *newAt = new Atom(0);
        setRlabel(newAt, userLabel);
        atomsToAdd.push_back(std::make_pair(atom, newAt));
      }
    }

    // Deal with non-user supplied labels
    for (int indexLabel : indexLabels) {
      auto atm = atoms.find(indexLabel);
      if (atm == atoms.end()) continue;  // label not used in the rgroup

      Atom *atom = atm->second;
      mappings[indexLabel] = ++count;
      if (atom->getAtomicNum() == 0) {  // add to dummy
        setRlabel(atom, count);
      } else {
        auto *newAt = new Atom(0);
        setRlabel(newAt, count);
        atomsToAdd.push_back(std::make_pair(atom, newAt));
      }
    }

    // Deal with multiple bonds to the same label
    for (auto &extraAtomRLabel : extraAtomRLabels) {
      auto atm = atoms.find(extraAtomRLabel.first);
      if (atm == atoms.end()) continue;  // label not used in the rgroup
      Atom *atom = atm->second;

      for (size_t i = 0; i < extraAtomRLabel.second.size(); ++i) {
        extraAtomRLabel.second[i] = ++count;
        // Is this necessary?
        PRECONDITION(
            atom->getAtomicNum() > 1,
            "Multiple attachements to a dummy (or hydrogen) is weird.");
        auto *newAt = new Atom(0);
        setRlabel(newAt, count);
        atomsToAdd.push_back(std::make_pair(atom, newAt));
      }
    }

    for (auto &i : atomsToAdd) {
      mol.addAtom(i.second, false, true);
      mol.addBond(i.first, i.second, Bond::SINGLE);
    }
    mol.updatePropertyCache(false);  // this was github #1550
  }

  void relabelRGroup(RGroupData &rgroup, const std::map<int, int> &mappings) {
    PRECONDITION(rgroup.combinedMol.get(), "Unprocessed rgroup");

    RWMol &mol = *rgroup.combinedMol.get();
    if (rgroup.combinedMol->hasProp(done)) {
      rgroup.labelled = true;
      return;
    }

    mol.setProp(done, true);

    std::vector<std::pair<Atom *, Atom *>> atomsToAdd;  // adds -R if necessary

    for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      Atom *atom = *atIt;
      if (atom->hasProp(SIDECHAIN_RLABELS)) {
        atom->setIsotope(0);
        const std::vector<int> &rlabels =
            atom->getProp<std::vector<int>>(SIDECHAIN_RLABELS);
        // switch on atom mappings or rlabels....

        for (int rlabel : rlabels) {
          auto label = mappings.find(rlabel);
          PRECONDITION(label != mappings.end(), "Unprocessed mapping");

          if (atom->getAtomicNum() == 0) {
            setRlabel(atom, label->second);
          } else {
            auto *newAt = new Atom(0);
            setRlabel(newAt, label->second);
            atomsToAdd.push_back(std::make_pair(atom, newAt));
          }
        }
      }
    }

    for (auto &i : atomsToAdd) {
      mol.addAtom(i.second, false, true);
      mol.addBond(i.first, i.second, Bond::SINGLE);
    }

    if (params.removeHydrogensPostMatch) {
      bool implicitOnly = false;
      bool updateExplicitCount = false;
      bool sanitize = false;
      MolOps::removeHs(mol, implicitOnly, updateExplicitCount, sanitize);
    }

    mol.updatePropertyCache(false);  // this was github #1550
    rgroup.labelled = true;
  }

  // relabel the core and sidechains using the specified user labels
  //  if matches exist for non labelled atoms, these are added as well
  void relabel() {
    std::vector<RGroupMatch> best = GetCurrentBestPermutation();

    // get the labels used
    std::set<int> userLabels;
    std::set<int> indexLabels;

    // Go through all the RGroups and find out which labels were
    //  actually used.

    // some atoms will have multiple attachment points, i.e. cycles
    //  split these up into new rlabels if necessary
    //  These are detected at match time
    //  This vector will hold the extra (new) labels required
    std::map<int, std::vector<int>> extraAtomRLabels;

    for (auto &it : best) {
      for (auto rit = it.rgroups.begin(); rit != it.rgroups.end(); ++rit) {
        if (rit->first >= 0) userLabels.insert(rit->first);
        if (rit->first < 0) indexLabels.insert(rit->first);

        std::map<int, int> rlabelsUsedInRGroup =
            rit->second->getNumBondsToRlabels();
        for (auto &numBondsUsed : rlabelsUsedInRGroup) {
          // Make space for the extra labels
          if (numBondsUsed.second > 1) {  // multiple
            extraAtomRLabels[numBondsUsed.first].resize(numBondsUsed.second -
                                                        1);
          }
        }
      }
    }

    finalRlabelMapping.clear();
    for (std::map<int, RWMol>::const_iterator coreIt = cores.begin();
         coreIt != cores.end(); ++coreIt) {
      boost::shared_ptr<RWMol> labelledCore(new RWMol(coreIt->second));
      labelledCores[coreIt->first] = labelledCore;

      relabelCore(*labelledCore.get(), finalRlabelMapping, userLabels,
                  indexLabels, extraAtomRLabels);
    }

    for (auto &it : best) {
      for (auto rit = it.rgroups.begin(); rit != it.rgroups.end(); ++rit) {
        relabelRGroup(*rit->second, finalRlabelMapping);
      }
    }
  }

  bool process(bool pruneMatches, bool finalize = false) {
    if (matches.size() == 0) return false;

    // Exhaustive search, get the MxN matrix
    size_t M = matches.size();
    std::vector<size_t> permutations;
    size_t N = 1;

    for (size_t m = 0; m < M; ++m) {
      size_t sz = matches[m].size();
      permutations.push_back(sz);
      N *= sz;
    }

    permutation = std::vector<size_t>(permutations.size(), 0);

    // run through all possible matches and score each
    //  set
    double best_score = 0;
    std::vector<size_t> best_permutation = permutation;
    size_t count = 0;
#ifdef DEBUG
    std::cerr << "Processing" << std::endl;
#endif
    CartesianProduct iterator(permutations);
    while (iterator.next()) {
      if (count > N) throw ValueErrorException("Next did not finish");
#ifdef DEBUG
      std::cerr << "**************************************************"
                << std::endl;
#endif
      double newscore = score(iterator.permutation, matches, labels);
      if (newscore > best_score) {
#ifdef DEBUG
        std::cerr << " ===> current best:" << newscore << ">" << best_score
                  << std::endl;
#endif
        best_score = newscore;
        best_permutation = iterator.permutation;
      }
    }

    permutation = best_permutation;
    if (pruneMatches || finalize) {
      prune();
    }

    if (finalize) {
      relabel();
    }

    return true;
  }
};

RGroupDecomposition::RGroupDecomposition(
    const ROMol &inputCore, const RGroupDecompositionParameters &params)
    : data(new RGroupDecompData(inputCore, params)) {}

RGroupDecomposition::RGroupDecomposition(
    const std::vector<ROMOL_SPTR> &cores,
    const RGroupDecompositionParameters &params)
    : data(new RGroupDecompData(cores, params)) {}

RGroupDecomposition::~RGroupDecomposition() { delete data; }

int RGroupDecomposition::add(const ROMol &inmol) {
  // get the sidechains if possible
  //  Add hs for better symmeterization
  RWMol mol(inmol);
  bool explicitOnly = false;
  bool addCoords = true;
  MolOps::addHs(mol, explicitOnly, addCoords);

  int core_idx = 0;
  const RWMol *core = nullptr;
  std::vector<MatchVectType> tmatches;

  // Find the first matching core.
  for (std::map<int, RWMol>::const_iterator coreIt = data->cores.begin();
       coreIt != data->cores.end(); ++coreIt) {
    {
      const bool uniquify = false;
      const bool recursionPossible = false;
      const bool useChirality = true;
      SubstructMatch(mol, coreIt->second, tmatches, uniquify, recursionPossible,
                     useChirality);
    }

    if (data->params.onlyMatchAtRGroups) {
      // First find all the core atoms that have user
      //  label and but their indices into core_atoms_with_user_labels
      std::set<int> core_atoms_with_user_labels;

      for (auto atom : coreIt->second.atoms()) {
        if (atom->hasProp(RLABEL)) {
          core_atoms_with_user_labels.insert(atom->getIdx());
        }
      }

      std::vector<MatchVectType> tmatches_filtered;
      for (auto &mv : tmatches) {
        bool passes_filter = true;
        boost::dynamic_bitset<> target_match_indices(mol.getNumAtoms());
        for (auto &match : mv) {
          target_match_indices[match.second] = 1;
        }

        for (auto &match : mv) {
          const Atom *atm = mol.getAtomWithIdx(match.second);
          // is this a labelled rgroup or not?
          if (core_atoms_with_user_labels.find(match.first) ==
              core_atoms_with_user_labels.end()) {
            // nope... if any neighbor is not part of the substructure
            //  make sure we are a hydrogen, otherwise, skip the match
            for (const auto &nbri :
                 boost::make_iterator_range(mol.getAtomNeighbors(atm))) {
              const auto &nbr = mol[nbri];
              if (nbr->getAtomicNum() != 1 &&
                  !target_match_indices[nbr->getIdx()]) {
                passes_filter = false;
                break;
              }
            }
          }
          if (!passes_filter) break;
        }

        if (passes_filter) {
          tmatches_filtered.push_back(mv);
        }
      }
      tmatches = tmatches_filtered;
    }

    if (!tmatches.size()) {
      continue;
    } else {
      if (tmatches.size() > 1) {
        if (data->matches.size() == 0) {
          // Greedy strategy just grabs the first match and
          //  takes the best matches from the rest
          if (data->params.matchingStrategy == Greedy) tmatches.resize(1);
        }
      }
      core = &coreIt->second;
      core_idx = coreIt->first;
      break;
    }
  }

  if (core == nullptr) return -1;

  // strategies
  // ==========
  // Exhaustive - saves all matches and optimizes later exhaustive
  //               May never finish due to combinitorial complexity
  // Greedy - matches to *FIRST* available match
  // GreedyChunks - default - process every N chunks

  //  Should probably scan all mols first to find match with
  //  smallest number of matches...
  size_t size = data->matches.size();

  std::vector<RGroupMatch> potentialMatches;

  for (auto &tmatche : tmatches) {
    boost::scoped_ptr<ROMol> tMol;
    {
      const bool replaceDummies = false;
      const bool labelByIndex = true;
      const bool requireDummyMatch = false;
      tMol.reset(replaceCore(mol, *core, tmatche, replaceDummies, labelByIndex,
                             requireDummyMatch));
    }

    if (tMol) {
      R_DECOMP match;
      // rlabel rgroups
      MOL_SPTR_VECT fragments = MolOps::getMolFrags(*tMol, false);
      for (size_t i = 0; i < fragments.size(); ++i) {
        std::vector<int> attachments;
        boost::shared_ptr<ROMol> &newMol = fragments[i];
        newMol->setProp<int>("core", core_idx);
        newMol->setProp<int>("idx", size);
        newMol->setProp<int>("frag_idx", i);

        for (ROMol::AtomIterator atIt = newMol->beginAtoms();
             atIt != newMol->endAtoms(); ++atIt) {
          Atom *tmp = *atIt;
          unsigned int elno = tmp->getAtomicNum();
          if (elno == 0) {
            unsigned int index =
                tmp->getIsotope();  // this is the index into the core
            // it messes up when there are multiple ?
            int rlabel;
            if (core->getAtomWithIdx(index)->getPropIfPresent(RLABEL, rlabel)) {
              std::vector<int> rlabelsOnSideChain;
              tmp->getPropIfPresent(SIDECHAIN_RLABELS, rlabelsOnSideChain);
              rlabelsOnSideChain.push_back(rlabel);
              tmp->setProp(SIDECHAIN_RLABELS, rlabelsOnSideChain);

              data->labels.insert(rlabel);  // keep track of all labels used
              attachments.push_back(rlabel);
            }
          }
        }

        if (attachments.size() > 0) {
          // reject multiple attachments?
          // what to do with labelled cores ?
          std::string newCoreSmi = MolToSmiles(*newMol, true);

          for (size_t attach_idx = 0; attach_idx < attachments.size();
               ++attach_idx) {
            int rlabel = attachments[attach_idx];
            ADD_MATCH(match, rlabel);
            match[rlabel]->add(newMol, attachments);
          }
        } else {
          // special case, only one fragment
          if (fragments.size() == 1) {  // need to make a new core
            // remove the sidechains
            RWMol newCore(mol);

            for (MatchVectType::const_iterator mvit = tmatche.begin();
                 mvit != tmatche.end(); ++mvit) {
              const Atom *coreAtm = core->getAtomWithIdx(mvit->first);
              Atom *newCoreAtm = newCore.getAtomWithIdx(mvit->second);
              int rlabel;
              if (coreAtm->getPropIfPresent(RLABEL, rlabel)) {
                newCoreAtm->setProp<int>(RLABEL, rlabel);
              }
              newCoreAtm->setProp<bool>("keep", true);
            }

            for (int aIdx = newCore.getNumAtoms() - 1; aIdx >= 0; --aIdx) {
              Atom *atom = newCore.getAtomWithIdx(aIdx);
              if (!atom->hasProp("keep")) newCore.removeAtom(atom);
            }
            if (newCore.getNumAtoms()) {
              std::string newCoreSmi = MolToSmiles(newCore, true);
              // add a new core if possible
              auto newcore = data->newCores.find(newCoreSmi);
              int core_idx = 0;
              if (newcore == data->newCores.end()) {
                core_idx = data->newCores[newCoreSmi] = data->newCoreLabel--;
                data->cores[core_idx] = newCore;
                return add(inmol);
              }
            }
          }
        }
      }

      if (match.size()) {
        potentialMatches.push_back(RGroupMatch(core_idx, match));
      }
    }
  }
  if (potentialMatches.size() == 0) {
    BOOST_LOG(rdWarningLog)
        << "No attachment points in side chains" << std::endl;

    return -1;
  }

  size_t N = 1;
  for (auto &matche : data->matches) {
    size_t sz = matche.size();
    N *= sz;
  }
  // oops, exponential is a pain
  if (N * potentialMatches.size() > 100000) {
    data->permutation = std::vector<size_t>(data->matches.size(), 0);
    data->process(true);
  }

  data->matches.push_back(potentialMatches);
  data->permutation = std::vector<size_t>(data->matches.size(), 0);

  if (size) {
    if (data->params.matchingStrategy & Greedy ||
        (data->params.matchingStrategy & GreedyChunks && size > 1 &&
         size % data->params.chunkSize == 0))
      data->process(true);
  }
  return data->matches.size() - 1;
}

bool RGroupDecomposition::process() {
  try {
    const bool prune = true;
    const bool finalize = true;
    return data->process(prune, finalize);
  } catch (...) {
    return false;
  }
}

RGroupRows RGroupDecomposition::getRGroupsAsRows() const {
  std::vector<RGroupMatch> permutation = data->GetCurrentBestPermutation();

  RGroupRows groups;

  int molidx = 0;
  for (auto it = permutation.begin(); it != permutation.end(); ++it, ++molidx) {
    // make a new rgroup entry
    groups.push_back(RGroupRow());
    RGroupRow &out_rgroups = groups.back();
    out_rgroups["Core"] = data->labelledCores[it->core_idx];

    R_DECOMP &in_rgroups = it->rgroups;

    for (R_DECOMP::const_iterator rgroup = in_rgroups.begin();
         rgroup != in_rgroups.end(); ++rgroup) {
      std::map<int, int>::const_iterator realLabel =
          data->finalRlabelMapping.find(rgroup->first);
      PRECONDITION(realLabel != data->finalRlabelMapping.end(),
                   "unprocessed rlabel, please call process() first.");
      out_rgroups[std::string("R") + std::to_string(realLabel->second)] =
          rgroup->second->combinedMol;
    }
  }
  return groups;
}
//! return rgroups in column order group[attachment_point][molidx] = ROMol
RGroupColumns RGroupDecomposition::getRGroupsAsColumns() const {
  std::vector<RGroupMatch> permutation = data->GetCurrentBestPermutation();

  RGroupColumns groups;

  unsigned int molidx = 0;
  for (auto it = permutation.begin(); it != permutation.end(); ++it, ++molidx) {
    R_DECOMP &in_rgroups = it->rgroups;
    groups["Core"].push_back(data->labelledCores[it->core_idx]);

    for (R_DECOMP::const_iterator rgroup = in_rgroups.begin();
         rgroup != in_rgroups.end(); ++rgroup) {
      std::map<int, int>::const_iterator realLabel =
          data->finalRlabelMapping.find(rgroup->first);
      PRECONDITION(realLabel != data->finalRlabelMapping.end(),
                   "unprocessed rlabel, please call process() first.");
      PRECONDITION(rgroup->second->combinedMol->hasProp(done),
                   "Not done! Call process()");

      std::string r = std::string("R") + std::to_string(realLabel->second);
      RGroupColumn &col = groups[r];
      if (molidx && col.size() < (size_t)(molidx - 1)) col.resize(molidx - 1);
      col.push_back(rgroup->second->combinedMol);
    }
  }
  // Now make all columns equal - this adds empty mols...
  for (auto &group : groups) {
    if (group.second.size() != molidx) {
      group.second.resize(molidx);
    }

    for (size_t idx = 0; idx < group.second.size(); ++idx) {
      if (!group.second[idx].get()) {
        group.second[idx] = boost::make_shared<RWMol>();
      }
    }
  }
  return groups;
}

namespace {
std::vector<unsigned int> Decomp(RGroupDecomposition &decomp,
                                 const std::vector<ROMOL_SPTR> &mols) {
  std::vector<unsigned int> unmatched;
  for (size_t i = 0; i < mols.size(); ++i) {
    int v = decomp.add(*mols[i].get());
    if (v == -1) unmatched.push_back(i);
  }
  decomp.process();
  return unmatched;
}
}  // namespace
unsigned int RGroupDecompose(const std::vector<ROMOL_SPTR> &cores,
                             const std::vector<ROMOL_SPTR> &mols,
                             RGroupRows &rows,
                             std::vector<unsigned int> *unmatchedIndices,
                             const RGroupDecompositionParameters &options) {
  RGroupDecomposition decomp(cores, options);
  std::vector<unsigned int> unmatched = Decomp(decomp, mols);
  if (unmatchedIndices) *unmatchedIndices = unmatched;
  rows = decomp.getRGroupsAsRows();
  return mols.size() - unmatched.size();
}

unsigned int RGroupDecompose(const std::vector<ROMOL_SPTR> &cores,
                             const std::vector<ROMOL_SPTR> &mols,
                             RGroupColumns &columns,
                             std::vector<unsigned int> *unmatchedIndices,
                             const RGroupDecompositionParameters &options) {
  RGroupDecomposition decomp(cores, options);
  std::vector<unsigned int> unmatched = Decomp(decomp, mols);
  if (unmatchedIndices) *unmatchedIndices = unmatched;
  columns = decomp.getRGroupsAsColumns();
  return mols.size() - unmatched.size();
}
}  // namespace RDKit
