//
// Created by gareth on 2/10/19.
//
#include <RDGeneral/export.h>

#ifndef RDKIT_CHUNKED_HITLIST_H
#define RDKIT_CHUNKED_HITLIST_H

#include <RDGeneral/Exceptions.h>
#include "SubstructLibrary.h"
#include <queue>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <atomic>

namespace RDKit {

template <typename T>
class ChunkedHitlist: boost::noncopyable {
  std::queue<T> queue;
  std::mutex mutex;
  const unsigned int chunkSize;
  std::condition_variable cv;
  std::atomic<bool> finished{false};
  const unsigned int maxHits;
  unsigned int noHits = 0;

  bool returnNext() const {
    bool f = finished.load();
    return queue.size() > chunkSize || f;
  }

  std::vector<T> getChunk() {
    std::vector<T> rtn;
    if (maxHits > 0 && noHits >= maxHits) {
      return rtn;
    }
    rtn.reserve(chunkSize);
    while (!queue.empty()) {
      T elem = queue.front();
      rtn.push_back(elem);
      queue.pop();
      noHits++;
      if (maxHits > 0 && noHits >= maxHits) {
        finished = true;
        break;
      }
      if (rtn.size() >= chunkSize) {
        break;
      }
    }
    return rtn;
  }

 public:
  ChunkedHitlist() : chunkSize{100}, maxHits{0} {};

  ChunkedHitlist(unsigned int cs, unsigned int mh)
      : chunkSize{cs}, maxHits{mh} {};

  void push(const T &item) {
    std::lock_guard<std::mutex> lock(mutex);
    queue.push(item);
    if (queue.size() >= chunkSize) {
      cv.notify_one();
    }
  }

  const std::vector<T> next() {
    std::unique_lock<std::mutex> lock(mutex);
    if (finished.load() || queue.size() > chunkSize) {
      return getChunk();
    }
    cv.wait(lock, [this] { return finished.load() || queue.size() >= chunkSize; });
    return getChunk();
  }

  void setFinished() {
    finished = true;
    cv.notify_one();
  }

  bool isFinished() const { return finished.load(); }
};

}  // namespace RDKit

#endif
