#ifdef RDK_BUILD_THREADSAFE_SSS
#include <catch2/catch_all.hpp>
#include <RDGeneral/RDLog.h>

#include <functional>
#include <iomanip>
#include <sstream>

#include "ConcurrentQueue.h"
#include "ControlCHandler.h"

using namespace RDKit;

namespace {
volatile std::sig_atomic_t previousSigIntCount = 0;

void countPreviousSigInt(int signalNumber) {
  if (signalNumber == SIGINT) {
    previousSigIntCount = previousSigIntCount + 1;
  }
}

void exerciseControlCHandler(std::atomic<bool> &start) {
  while (!start.load(std::memory_order_acquire)) {
    std::this_thread::yield();
  }
  for (unsigned int i = 0; i < 100; ++i) {
    ControlCHandler handler;
    (void)ControlCHandler::getGotSignal();
  }
}
}  // namespace

TEST_CASE("ControlCHandler preserves nested signal handlers") {
  previousSigIntCount = 0;
  const auto previousHandler = std::signal(SIGINT, countPreviousSigInt);
  REQUIRE(previousHandler != SIG_ERR);
  {
    ControlCHandler outer;
    {
      ControlCHandler inner;
      std::raise(SIGINT);
      CHECK(ControlCHandler::getGotSignal());
      std::raise(SIGINT);
      CHECK(previousSigIntCount == 1);
      ControlCHandler::reset();
      CHECK_FALSE(ControlCHandler::getGotSignal());
      std::raise(SIGINT);
      CHECK(ControlCHandler::getGotSignal());
    }
    CHECK(ControlCHandler::getGotSignal());
    ControlCHandler::reset();
    std::raise(SIGINT);
    CHECK(ControlCHandler::getGotSignal());
  }
  const auto restoredHandler = std::signal(SIGINT, previousHandler);
  CHECK(restoredHandler == countPreviousSigInt);
  CHECK_FALSE(ControlCHandler::getGotSignal());
}

TEST_CASE("ControlCHandler reset before a scope retains the original handler") {
  previousSigIntCount = 0;
  const auto previousHandler = std::signal(SIGINT, countPreviousSigInt);
  REQUIRE(previousHandler != SIG_ERR);
  ControlCHandler::reset();
  ControlCHandler::reset();
  {
    ControlCHandler handler;
    std::raise(SIGINT);
    CHECK(ControlCHandler::getGotSignal());
    std::raise(SIGINT);
    CHECK(previousSigIntCount == 1);
  }
  const auto restoredHandler = std::signal(SIGINT, previousHandler);
  CHECK(restoredHandler == countPreviousSigInt);
}

TEST_CASE("ControlCHandler supports concurrent lifetimes") {
  std::atomic<bool> start{false};
  std::vector<std::thread> workers;
  for (unsigned int i = 0; i < 8; ++i) {
    workers.emplace_back(exerciseControlCHandler, std::ref(start));
  }
  start.store(true, std::memory_order_release);
  for (auto &worker : workers) {
    worker.join();
  }
  CHECK_FALSE(ControlCHandler::getGotSignal());
}

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
