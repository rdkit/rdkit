#ifdef RDK_THREADSAFE_SSS
#ifndef CONCURRENT_QUEUE
#define CONCURRENT_QUEUE
#include <condition_variable>
#include <iostream>
#include <queue>
#include <thread>
#include <vector>

namespace RDKit {
template <class E>
class ConcurrentQueue {
 private:
  size_t capacity;
  bool done;
  std::mutex lock;
  std::condition_variable cv_push, cv_pop;
  std::queue<E> q;

 public:
  ConcurrentQueue<E>(size_t capacity) : capacity(capacity), done(false) {}
  void push(const E& element);
  bool pop(E& element);
  bool isEmpty();
  bool getDone();
  void setDone();
  size_t limit();
  size_t size();
};

template <class E>
void ConcurrentQueue<E>::push(const E& element) {
  std::unique_lock<std::mutex> lk(lock);
  while (q.size() >= capacity) {
    cv_push.wait(lk);
  }
  q.push(element);
  cv_pop.notify_one();
}

template <class E>
bool ConcurrentQueue<E>::pop(E& element) {
  std::unique_lock<std::mutex> lk(lock);
  while (q.empty()) {
    if (done) {
      return false;
    }
    cv_pop.wait(lk);
  }
  element = q.front();
  q.pop();
  cv_push.notify_one();
  return true;
}

template <class E>
bool ConcurrentQueue<E>::isEmpty() {
  std::unique_lock<std::mutex> lk(lock);
  return q.empty();
}

template <class E>
void ConcurrentQueue<E>::setDone() {
  std::unique_lock<std::mutex> lk(lock);
  done = true;
  cv_pop.notify_all();
}

template <class E>
bool ConcurrentQueue<E>::getDone() {
  std::unique_lock<std::mutex> lk(lock);
  return done;
}

template <class E>
size_t ConcurrentQueue<E>::limit() {
  return capacity;
}

template <class E>
size_t ConcurrentQueue<E>::size() {
  std::unique_lock<std::mutex> lk(lock);
  return q.size();
}

}  // namespace RDKit
#endif
#endif
