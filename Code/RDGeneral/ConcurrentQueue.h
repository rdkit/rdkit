#ifdef RDK_THREADSAFE_SSS
#ifndef CONCURRENT_QUEUE
#define CONCURRENT_QUEUE
#include <condition_variable>
#include <iostream>
#include <thread>
#include <vector>

namespace RDKit {
template <class E>
class ConcurrentQueue {
 public:
  int capacity, head, tail;
  std::vector<E> elements;
  std::mutex lock;
  std::condition_variable notEmpty, notFull;
  ConcurrentQueue<E>(int capacity) {
    this->capacity = capacity;
    std::vector<E> elements(capacity);
    this->elements = elements;
    this->head = 0;
    this->tail = 0;
  }
  void push(E element);
  E pop();
  bool isEmpty();
};

template <class E>
void ConcurrentQueue<E>::push(E element) {
  std::unique_lock<std::mutex> lk(lock);
  while (this->head + this->capacity == this->tail) {
    notFull.wait(lk);
  }
  bool wasEmpty = (head == tail);
  this->elements.at(this->tail % this->capacity) = element;
  this->tail++;
  //! the queue was empty before but affter adding the element
  //! it is not anymore, thus notify all the consumers to consume
  if (wasEmpty) {
    notEmpty.notify_all();
  }
}

template <class E>
E ConcurrentQueue<E>::pop() {
  std::unique_lock<std::mutex> lk(lock);
  while (this->head == this->tail) {
    notEmpty.wait(lk);
  }
  bool wasFull = (head + capacity == tail);
  E element = this->elements.at(this->head % this->capacity);
  this->head++;
  //! the queue was full before but now it is not since we have consumed the
  //! element thus notify the producer to producer, i.e. put an element in the
  //! queue
  if (wasFull) {
    notFull.notify_all();
  }
  return element;
}

template <class E>
bool ConcurrentQueue<E>::isEmpty() {
  std::unique_lock<std::mutex> lk(lock);
  return (this->head == this->tail);
}

}  // namespace RDKit
#endif
#endif
