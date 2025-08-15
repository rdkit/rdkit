#ifdef RDK_BUILD_THREADSAFE_SSS
#include <catch2/catch_all.hpp>
#include <RDGeneral/RDLog.h>

#include <functional>
#include <iomanip>
#include <sstream>

#include "ConcurrentQueue.h"

using namespace RDKit;

TEST_CASE("testPushAndPop") {
  ConcurrentQueue<int> *q = new ConcurrentQueue<int>(4);
  int e1, e2, e3;
  REQUIRE(q->isEmpty());

  q->push(1);
  q->push(2);
  q->push(3);

  REQUIRE(!q->isEmpty());

  REQUIRE(q->pop(e1));
  REQUIRE(q->pop(e2));
  REQUIRE(q->pop(e3));

  REQUIRE(e1 == 1);
  REQUIRE(e2 == 2);
  REQUIRE(e3 == 3);

  REQUIRE(q->isEmpty());

  delete (q);
}

void produce(ConcurrentQueue<int> &q, const int numToProduce) {
  for (int i = 0; i < numToProduce; ++i) {
    q.push(i);
  }
}

void consume(ConcurrentQueue<int> &q, std::vector<int> &result) {
  int element;
  while (q.pop(element)) {
    result.push_back(element);
  }
}

//! multithreaded testing for ConcurrentQueue
bool testProducerConsumer(const int numProducerThreads,
                          const int numConsumerThreads) {
  ConcurrentQueue<int> q(5);
  REQUIRE(q.isEmpty());

  const int numToProduce = 10;

  std::vector<std::thread> producers(numProducerThreads);
  std::vector<std::thread> consumers(numConsumerThreads);
  std::vector<std::vector<int>> results(numConsumerThreads);

  //! start producer threads
  for (int i = 0; i < numProducerThreads; i++) {
    producers[i] = std::thread(produce, std::ref(q), numToProduce);
  }
  //! start consumer threads
  for (int i = 0; i < numConsumerThreads; i++) {
    consumers[i] = std::thread(consume, std::ref(q), std::ref(results[i]));
  }

  std::for_each(producers.begin(), producers.end(),
                std::mem_fn(&std::thread::join));

  //! the producer is done producing
  q.setDone();

  std::for_each(consumers.begin(), consumers.end(),
                std::mem_fn(&std::thread::join));
  REQUIRE(q.isEmpty());

  std::vector<int> frequency(numToProduce, 0);
  for (auto &result : results) {
    for (auto &element : result) {
      frequency[element] += 1;
    }
  }
  for (auto &freq : frequency) {
    if (freq != numProducerThreads) {
      return false;
    }
  }
  return true;
}

TEST_CASE("testMultipleTimes") {
  const int trials = 10000;
  //! Single Producer, Single Consumer
  for (int i = 0; i < trials; i++) {
    REQUIRE(testProducerConsumer(1, 1));
  }

  //! Single Producer, Multiple Consumer
  for (int i = 0; i < trials; i++) {
    REQUIRE(testProducerConsumer(1, 5));
  }

  //! Multiple Producer, Single Consumer
  for (int i = 0; i < trials; i++) {
    REQUIRE(testProducerConsumer(5, 1));
  }

  //! Multiple Producer, Multiple Consumer
  for (int i = 0; i < trials; i++) {
    REQUIRE(testProducerConsumer(2, 4));
  }
}
#endif
