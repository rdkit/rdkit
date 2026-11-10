#ifdef RDK_BUILD_THREADSAFE_SSS
#include <catch2/catch_all.hpp>
#include <RDGeneral/RDLog.h>

#include <functional>
#include <chrono>
#include <iomanip>
#include <memory>
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

#ifdef RDK_TEST_WRAPPED_SIGNAL
std::atomic_flag interruptInstallation{};
std::atomic_flag pauseRestoration{};
std::atomic_flag restorationEntered{};
std::atomic_flag releaseRestoration{};
std::atomic_flag replacementReady{};
std::atomic_flag finishReplacement{};

void ignorePreviousSigInt(int) {}

bool waitForFlag(const std::atomic_flag &flag) {
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
  while (!flag.test()) {
    if (std::chrono::steady_clock::now() >= deadline) {
      return false;
    }
    std::this_thread::yield();
  }
  return true;
}

void replaceControlCHandler(std::unique_ptr<ControlCHandler> first) {
  first.reset();
  ControlCHandler second;
  replacementReady.test_and_set();
  while (!finishReplacement.test()) {
    std::this_thread::yield();
  }
}
#endif
}  // namespace

#ifdef RDK_TEST_WRAPPED_SIGNAL
using SignalHandler = void (*)(int);
extern "C" SignalHandler __real_signal(int, SignalHandler);
extern "C" SignalHandler __wrap_signal(int number, SignalHandler handler) {
  if (number == SIGINT && handler == ignorePreviousSigInt &&
      pauseRestoration.test() && !restorationEntered.test_and_set()) {
    while (!releaseRestoration.test()) {
    }
  }
  const auto previous = __real_signal(number, handler);
  if (number == SIGINT && handler == ControlCHandler::signalHandler &&
      interruptInstallation.test()) {
    interruptInstallation.clear();
    std::raise(SIGINT);
  }
  return previous;
}

TEST_CASE("ControlCHandler handles an interrupt during installation") {
  const auto original = __real_signal(SIGINT, ignorePreviousSigInt);
  {
    interruptInstallation.test_and_set();
    ControlCHandler handler;
    const auto restored = __real_signal(SIGINT, ignorePreviousSigInt);
    CHECK(ControlCHandler::getGotSignal());
    CHECK(restored == ignorePreviousSigInt);
  }
  __real_signal(SIGINT, original);
}

TEST_CASE("ControlCHandler finishes restoration before replacing a scope") {
  const auto original = __real_signal(SIGINT, ignorePreviousSigInt);
  restorationEntered.clear();
  releaseRestoration.clear();
  replacementReady.clear();
  finishReplacement.clear();
  auto first = std::make_unique<ControlCHandler>();
  pauseRestoration.test_and_set();
  std::thread interrupted([] { std::raise(SIGINT); });
  const bool entered = waitForFlag(restorationEntered);
  std::thread replacement(replaceControlCHandler, std::move(first));
  // Give a replacement scope a chance to overtake the paused restoration.
  const auto deadline = std::chrono::steady_clock::now() +
                        std::chrono::milliseconds(100);
  while (!replacementReady.test() && std::chrono::steady_clock::now() < deadline) {
    std::this_thread::yield();
  }
  releaseRestoration.test_and_set();
  interrupted.join();
  const bool ready = waitForFlag(replacementReady);
  const auto current = __real_signal(SIGINT, ControlCHandler::signalHandler);
  const bool cancelled = ControlCHandler::getGotSignal();
  finishReplacement.test_and_set();
  replacement.join();
  pauseRestoration.clear();
  __real_signal(SIGINT, original);
  CHECK(entered);
  CHECK(ready);
  CHECK(current == ControlCHandler::signalHandler);
  CHECK_FALSE(cancelled);
}
#endif

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
