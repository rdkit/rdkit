#ifdef RDK_THREADSAFE_SSS
#ifndef PRODUCER_CONSUMER_H
#define PRODUCER_CONSUMER_H
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>
namespace RDKit {
template <class E>
class ProducerConsumer {
 private:
  std::mutex q_mutex;
  std::queue<E> q;
  std::condition_variable cv_produce;
  std::condition_variable cv_consume;
  bool done = false;

 public:
  void produce(size_t items, size_t stock) {
    for (size_t i = 0; i < items; ++i) {
      std::unique_lock<std::mutex> lock(q_mutex);
      cv_produce.wait(lock, [&] { return q.size() < stock; });
      //! in general we can produce the data here
      q.push(i);
      cv_consume.notify_all();
    }
  }

  void consume(size_t& sum) {
    while (!done || !q.empty()) {
      std::unique_lock<std::mutex> lock(q_mutex);
      if (cv_consume.wait_for(lock, std::chrono::seconds(1),
                              [&] { return !q.empty(); })) {
        sum += q.front();
        q.pop();
        cv_produce.notify_all();
      }
    }
  }

  void setDone() {
    std::unique_lock<std::mutex> lock(q_mutex);
    done = true;
  }

  bool isDone() {
    std::unique_lock<std::mutex> lock(q_mutex);
    return done;
  }
};
}  // namespace RDKit
#endif
#endif
