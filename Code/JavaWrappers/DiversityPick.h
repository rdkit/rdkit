#include <RDGeneral/export.h>
#include <list>
#include <map>
#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <RDGeneral/Exceptions.h>

namespace {
class taniFunctor {
 public:
  taniFunctor(const std::vector<ExplicitBitVect> &ebvs, bool useCache)
      : df_useCache(useCache), d_ebvs(ebvs) {}
  double operator()(unsigned int i, unsigned int j) {
    double res;
    if (df_useCache) {
      std::pair<unsigned int, unsigned int> idxPair(i, j);
      if (this->d_cache.count(idxPair) > 0) {
        res = this->d_cache[idxPair];
      } else {
        res = 1. - TanimotoSimilarity(d_ebvs[i], d_ebvs[j]);
        this->d_cache[idxPair] = res;
      }
    } else {
      res = 1. - TanimotoSimilarity(d_ebvs[i], d_ebvs[j]);
    }
    return res;
  }

 private:
  bool df_useCache;
  const std::vector<ExplicitBitVect> &d_ebvs;
  std::map<std::pair<unsigned int, unsigned int>, double> d_cache;
};
}  // namespace

std::vector<int> pickUsingFingerprints(
    const std::vector<ExplicitBitVect> &ebvs, unsigned int nToPick,
    int seed = -1, std::vector<int> firstPicks = std::vector<int>(),
    bool useCache = true) {
  if (nToPick >= ebvs.size())
    throw ValueErrorException("nToPick is larger than the vector size");
  std::vector<int> res;

  RDPickers::MaxMinPicker picker;
  taniFunctor ftor(ebvs, useCache);
  res = picker.lazyPick(ftor, ebvs.size(), nToPick, firstPicks, seed);
  return res;
}
