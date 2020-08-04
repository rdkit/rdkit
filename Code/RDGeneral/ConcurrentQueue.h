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
#include <thread>
#include <vector>

namespace RDKit {
template <typename E>
class ConcurrentQueue {
 private:
  unsigned int d_capacity;
  bool d_done;
  std::vector<E> d_elements;
  unsigned int d_head, d_tail;
  mutable std::mutex d_lock;
  std::condition_variable d_notEmpty, d_notFull;

 private:
  ConcurrentQueue<E>(const ConcurrentQueue<E>&);
  ConcurrentQueue<E>& operator=(const ConcurrentQueue<E>&);

 public:
  ConcurrentQueue<E>(unsigned int capacity)
      : d_capacity(capacity), d_done(false), d_head(0), d_tail(0) {
    std::vector<E> elements(capacity);
    d_elements = elements;
  }

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

  //! clears the vector
  void clear();
};

template <typename E>
void ConcurrentQueue<E>::push(const E& element) {
  std::unique_lock<std::mutex> lk(d_lock);
  //! concurrent queue is full so we wait until
  //! it is not full
  while (d_head + d_capacity == d_tail) {
    d_notFull.wait(lk);
  }
  bool wasEmpty = (d_head == d_tail);
  d_elements.at(d_tail % d_capacity) = element;
  d_tail++;
  //! if the concurrent queue was empty before
  //! then it is not any more since we have "pushed" an element
  //! thus we notify all the consumer threads
  if (wasEmpty) {
    d_notEmpty.notify_all();
  }
}

template <typename E>
bool ConcurrentQueue<E>::pop(E& element) {
  std::unique_lock<std::mutex> lk(d_lock);
  //! concurrent queue is empty so we wait until
  //! it is not empty
  while (d_head == d_tail) {
    if (d_done) {
      return false;
    }
    d_notEmpty.wait(lk);
  }
  bool wasFull = (d_head + d_capacity == d_tail);
  element = d_elements.at(d_head % d_capacity);
  d_head++;
  //! if the concurrent queue was full before
  //! then it is not any more since we have "popped" an element
  //! thus we notify all producer threads
  if (wasFull) {
    d_notFull.notify_all();
  }
  return true;
}

template <typename E>
bool ConcurrentQueue<E>::isEmpty() const {
  std::unique_lock<std::mutex> lk(d_lock);
  return (d_head == d_tail);
}

template <typename E>
bool ConcurrentQueue<E>::getDone() const {
  std::unique_lock<std::mutex> lk(d_lock);
  return d_done;
}

template <typename E>
void ConcurrentQueue<E>::setDone() {
  std::unique_lock<std::mutex> lk(d_lock);
  d_done = true;
  d_notEmpty.notify_all();
}

template <typename E>
void ConcurrentQueue<E>::clear() {
  std::unique_lock<std::mutex> lk(d_lock);
  d_elements.clear();
}

}  // namespace RDKit
#endif
#endif
