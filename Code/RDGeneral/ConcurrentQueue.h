#ifdef RDK_THREADSAFE_SSS
#ifndef CONCURRENT_QUEUE
#define CONCURRENT_QUEUE
#include <condition_variable>
#include <queue>
#include <thread>

namespace RDKit {
template <class E>
class ConcurrentQueue {
 private:
  size_t capacity;
  bool done;
  std::mutex lock;
  //! we need two condition variables to establish communication
  //! between popping (consumer threads) and pushing (producer threads)
  std::condition_variable goPush, goPop;
  std::queue<E> q;

 public:
  ConcurrentQueue<E>(size_t capacity) : capacity(capacity), done(false) {}
  //! tries to push an element into the queue if it is not full without
  //! modifying the variable element
  void push(const E& element);

  //! tries to pop an element from the queue if it is not empty and not done
  //! the boolean value indicates the whether popping is successful
  bool pop(E& element);

  //! checks whether the ConcurrentQueue is empty
  bool isEmpty();

  //! returns the value of the variable done
  bool getDone();

  //! sets the variable done = true
  void setDone();

  //! returns the capacity of the ConcurrentQueue
  size_t limit();

  //! returns the current size of the queue
  size_t size();
};

template <class E>
void ConcurrentQueue<E>::push(const E& element) {
  std::unique_lock<std::mutex> lk(lock);
  while (q.size() == capacity) {
    goPush.wait(lk);
  }
  q.push(element);
  goPop.notify_one();
}

template <class E>
bool ConcurrentQueue<E>::pop(E& element) {
  std::unique_lock<std::mutex> lk(lock);
  while (q.empty()) {
    if (done) {
      return false;
    }
    goPop.wait(lk);
  }
  element = q.front();
  q.pop();
  goPush.notify_one();
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
  goPop.notify_all();
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
