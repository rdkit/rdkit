//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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
  size_t d_capacity;
  bool d_done;
  mutable std::mutex d_lock;
  //! we need two condition variables to establish communication
  //! between popping (consumer threads) and pushing (producer threads)
  std::condition_variable d_goPush, d_goPop;
  std::queue<E> d_q;

 public:
  ConcurrentQueue<E>(size_t capacity) : d_capacity(capacity), d_done(false) {}
  //! tries to push an element into the queue if it is not full without
  //! modifying the variable element, if the queue is full then pushing an
  //! element will result in blocking
  void push(const E& element);

  //! tries to pop an element from the queue if it is not empty and not done
  //! the boolean value indicates the whether popping is successful, if the
  //! queue is empty and not done then popping an element will result in
  //! blocking
  bool pop(E& element);

  //! checks whether the ConcurrentQueue is empty
  bool isEmpty() const;

  //! returns the value of the variable done
  bool getDone() const;

  //! sets the variable d_done = true
  void setDone();

  //! returns the current size of the queue
  size_t size() const;
};

template <class E>
void ConcurrentQueue<E>::push(const E& element) {
  std::unique_lock<std::mutex> lk(d_lock);
  while (d_q.size() == d_capacity) {
    d_goPush.wait(lk);
  }
  d_q.push(element);
  d_goPop.notify_one();
}

template <class E>
bool ConcurrentQueue<E>::pop(E& element) {
  std::unique_lock<std::mutex> lk(d_lock);
  while (d_q.empty()) {
    if (d_done) {
      return false;
    }
    d_goPop.wait(lk);
  }
  element = d_q.front();
  d_q.pop();
  d_goPush.notify_one();
  return true;
}

template <class E>
bool ConcurrentQueue<E>::isEmpty() const {
  std::unique_lock<std::mutex> lk(d_lock);
  return d_q.empty();
}

template <class E>
bool ConcurrentQueue<E>::getDone() const {
  std::unique_lock<std::mutex> lk(d_lock);
  return d_done;
}

template <class E>
void ConcurrentQueue<E>::setDone() {
  std::unique_lock<std::mutex> lk(d_lock);
  d_done = true;
  d_goPop.notify_all();
}

template <class E>
size_t ConcurrentQueue<E>::size() const {
  std::unique_lock<std::mutex> lk(d_lock);
  return d_q.size();
}
}  // namespace RDKit
#endif
#endif
