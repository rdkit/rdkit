#include <RDGeneral/export.h>

#ifndef RDKIT_EDITABLESUBSTRUCTLIBRARY_H
#define RDKIT_EDITABLESUBSTRUCTLIBRARY_H

#include <RDGeneral/Exceptions.h>
#include "SubstructLibrary.h"
#include <vector>
#include <boost/thread.hpp>
#include <boost/bimap.hpp>
#include <boost/range/combine.hpp>
#include <algorithm>

#include <RDGeneral/RDThreads.h>
#ifdef RDK_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

#include <atomic>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "ChunkedHitlist.h"

#include <chrono>
#include <thread>

namespace RDKit {

template <typename Key, typename MolHolder, typename FpHolder>
class EditableSubstructLibrary : boost::noncopyable {
  template <typename T>
  struct KeyType {};

  Key nameToKey(KeyType<std::string>, const std::string &name) { return name; }

  Key nameToKey(KeyType<int>, const std::string &name) {
    return std::atoi(name.c_str());
  }

  void hitlistSubSearcher(const ROMol &in_query, const Bits &bits,
                          ChunkedHitlist<Key> &hitlist,
                          unsigned int threadGroupNo,
                          unsigned int numThreads) const {
    ROMol query(in_query);
    MatchVectType matchVect;
    for (unsigned int idx = threadGroupNo; idx < molHolder->size();
         idx += numThreads) {
      if (!bits.check(idx)) continue;
      const auto &m = molHolder->getMol(idx);
      const ROMol *mol = m.get();
      if (SubstructMatch(*mol, query, matchVect, bits.recursionPossible,
                         bits.useChirality, bits.useQueryQueryMatches)) {
        hitlist.push(indexToId(idx));
        //std::this_thread::sleep_for(std::chrono::milliseconds(500));
        if (hitlist.isFinished()) {
          break;
        }
      }
    }
  }

  void internalHitlistMatches(const ROMol &query,
                              boost::shared_ptr<ChunkedHitlist<Key>> hitlist,
                              bool recursionPossible, bool useChirality,
                              bool useQueryQueryMatches, int numThreads = -1) {
    boost::shared_lock<boost::shared_mutex> lock(access_mutex);
    if (numThreads == -1)
      numThreads = (int)getNumThreadsToUse(numThreads);
    else
      numThreads = std::min(numThreads, (int)getNumThreadsToUse(numThreads));

    if (molHolder->size() < static_cast<unsigned int>(numThreads))
      numThreads = molHolder->size();

    std::vector<std::future<void>> thread_group;

    Bits bits(fpHolder.get(), query, recursionPossible, useChirality,
              useQueryQueryMatches);

    for (int thread_group_idx = 0; thread_group_idx < numThreads;
         ++thread_group_idx) {
      thread_group.emplace_back(
          std::async(std::launch::async, [this, &query, &bits, hitlist,
                                          thread_group_idx, numThreads] {
            hitlistSubSearcher(query, bits, *hitlist, thread_group_idx,
                               numThreads);
          }));
    }
    for (auto &fut : thread_group) {
      fut.get();
    }
    hitlist->setFinished();
    delete bits.queryBits;
  }

 protected:
  boost::shared_mutex access_mutex;
  boost::shared_ptr<MolHolder> molHolder = boost::make_shared<MolHolder>();
  boost::shared_ptr<FpHolder> fpHolder = boost::make_shared<FpHolder>();
  SubstructLibrary substructLibrary{molHolder, fpHolder};
  boost::bimap<unsigned int, Key> ids;
  using idInfo = typename boost::bimap<unsigned int, Key>::value_type;

  template <typename T>
  static void removeItem(std::vector<T> &c, const unsigned int idx,
                         const bool quickRemove = true) {
    if (idx >= c.size()) throw IndexErrorException(idx);
    if (quickRemove) {
      c[idx] = c.back();
      c.pop_back();
    } else {
      c.erase(c.begin() + idx);
    }
  }

  void checkIdNotPresent(const Key &id) const {
    const auto idRight = ids.right;
    if (idRight.find(id) != idRight.end())
      throw std::invalid_argument(
          "Id " + id + " is already present in substructure library");
  }

  void addFingerprints(const std::vector<boost::shared_ptr<ROMol>> &mols) {
    for (auto &mol : mols) {
      const ExplicitBitVect *bitVector = fpHolder->makeFingerprint(*mol);
      fpHolder->addFingerprint(*bitVector);
      delete bitVector;
    }
  }

  void addFingerprints(const std::vector<std::string> &stringFingerprints) {
    for (auto &strFp : stringFingerprints) {
      const auto *bitVector = new ExplicitBitVect(strFp);
      fpHolder->addFingerprint(*bitVector);
      delete bitVector;
    }
  }

 public:
  EditableSubstructLibrary(){};

  ~EditableSubstructLibrary(){};

  unsigned int addMol(const ROMol &mol) {
    const std::vector<boost::shared_ptr<ROMol>> mols{
        boost::make_shared<ROMol>(mol)};
    return addMols(mols);
  }

  unsigned int addMols(const std::vector<boost::shared_ptr<ROMol>> &mols) {
    std::vector<Key> newIds;
    newIds.reserve(mols.size());
    for (auto mol : mols) {
      if (!mol->hasProp("_Name"))
        throw std::invalid_argument(
            "EditableSubstructLibrary: addMol: molecule does not have _Name "
            "property!");
      const auto name = mol->getProp<std::string>("_Name");
      const Key newId = nameToKey(KeyType<Key>(), name);
      newIds.push_back(newId);
    }

    return addMols(newIds, mols);
  }

  int addMol(const Key &id, const ROMol &mol) {
    const std::vector<boost::shared_ptr<ROMol>> mols{
        boost::make_shared<ROMol>(mol)};
    const std::vector<Key> newIds{id};
    return addMols(newIds, mols);
  }

  unsigned int addMols(const std::vector<Key> newIds,
                       const std::vector<boost::shared_ptr<ROMol>> &mols) {
    boost::unique_lock<boost::shared_mutex> lock(access_mutex);
    auto start = ids.size();
    boost::shared_ptr<ROMol> mol;
    Key newId;
    BOOST_FOREACH (boost::tie(newId, mol), boost::combine(newIds, mols)) {
      checkIdNotPresent(newId);
      ids.insert(idInfo(ids.size(), newId));
      molHolder->addMol(*mol);
      assert(molHolder->size() == ids.size());
    }
    addFingerprints(mols);
    auto nAdded = ids.size() - start;
    assert(nAdded == mols.size());
    return static_cast<unsigned int>(nAdded);
  }

  void removeMol(const Key &id) {
    const std::vector<Key> ids{id};
    removeMols(ids);
  }

  void removeMols(const std::vector<Key> idsToRemove) {
    boost::unique_lock<boost::shared_mutex> lock(access_mutex);
    std::vector<unsigned int> indices;
    indices.reserve(idsToRemove.size());
    std::transform(
        idsToRemove.cbegin(), idsToRemove.cend(), std::back_inserter(indices),
        [this](const std::string &id) { return this->idToIndex(id); });
    // descending sort
    std::sort(indices.rbegin(), indices.rend());
    auto idxToIds = ids.left;
    auto mols = &molHolder->mols;
    auto fps = &fpHolder->fps;
    for (const auto idx : indices) {
      if (idx >= mols->size()) {
        // this can only happen if we have duplicate ids
        throw std::invalid_argument("Duplicate id in remove molecule list!");
      }
      auto lastIdx = mols->size() - 1;
      auto lastId = idxToIds.at(lastIdx);
      idxToIds.erase(lastIdx);
      if (idx < mols->size() - 1) {
        idxToIds.replace_data(idxToIds.find(idx), lastId);
      }
      removeItem(*mols, idx);
      removeItem(*fps, idx);
    }
  }

  boost::shared_ptr<ROMol> getMol(const Key &id) {
    boost::shared_lock<boost::shared_mutex> lock(access_mutex);
    int idx = idToIndex(id);
    return molHolder->getMol(idx);
  }

  const std::vector<std::string> getIds() {
    boost::shared_lock<boost::shared_mutex> lock(access_mutex);
    std::vector<std::string> allIds;
    allIds.reserve(ids.size());
    BOOST_FOREACH (idInfo it, ids) { allIds.push_back(it.right); }
    return allIds;
  }

  unsigned int size() { return substructLibrary.size(); }

  std::vector<Key> getMatches(const ROMol &query, bool recursionPossible = true,
                              bool useChirality = true,
                              bool useQueryQueryMatches = false,
                              int numThreads = -1, int maxResults = -1) {
    boost::shared_lock<boost::shared_mutex> lock(access_mutex);
    auto indexMatches = substructLibrary.getMatches(
        query, recursionPossible, useChirality, useQueryQueryMatches,
        numThreads, maxResults);
    std::vector<Key> ids;
    ids.reserve(indexMatches.size());
    BOOST_FOREACH (int idx, indexMatches) { ids.push_back(indexToId(idx)); }
    return ids;
  }

  boost::shared_ptr<ChunkedHitlist<Key>> getHitlistMatches(
      const ROMol &query, bool recursionPossible = true,
      bool useChirality = true, bool useQueryQueryMatches = false,
      unsigned int chunkSize = 100, int numThreads = -1,
      unsigned int maxResults = 0) {
    auto hitlist =
        boost::make_shared<ChunkedHitlist<Key>>(chunkSize, maxResults);

    // std::async will block as the (unassigned) returned future will wait on its destructor,
    // so use thread then detach

    // The hitlist shared pointer is passed by value.  I thought that passing a 
    // const reference would be fine and it doesn't seem to cause a problem in c++,
    // but creates error in Java.  Not sure why that is.
    auto thread = std::thread([this, &query, hitlist, recursionPossible, useChirality,
                useQueryQueryMatches, numThreads] {
      internalHitlistMatches(query, hitlist, recursionPossible, useChirality,
                             useQueryQueryMatches, numThreads);
    });
    thread.detach();
    return hitlist;
  }

  ExplicitBitVect *makeFingerprint(const ROMol &mol) const  {
    return fpHolder->makeFingerprint(mol);
  }

  // These two guys are not protected by the lock

  const std::string &indexToId(const unsigned int index) const {
    return ids.left.at(index);
  }

  unsigned int idToIndex(const Key &id) const { 
    return ids.right.at(id); 
  }

};

}  // namespace RDKit
#endif  // RDKIT_EDITABLESUBSTRUCTLIBRARY_H
