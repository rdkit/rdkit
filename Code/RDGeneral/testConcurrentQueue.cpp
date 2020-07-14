#ifdef RDK_THREADSAFE_SSS
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <atomic>
#include <functional>
#include <iomanip>
#include <sstream>
#include <vector>

#include "ConcurrentQueue.h"

using namespace RDKit;

//! method for thread-safe printing, only for debugging
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void testPushAndPop() {
  ConcurrentQueue<int>* q = new ConcurrentQueue<int>(4);
  int e1, e2, e3;
  TEST_ASSERT(q->isEmpty());
  TEST_ASSERT(q->limit() == 4);

  q->push(1);
  q->push(2);
  q->push(3);

  TEST_ASSERT(!q->isEmpty());
  TEST_ASSERT(q->size() == 3);

  TEST_ASSERT(q->pop(e1));
  TEST_ASSERT(q->pop(e2));
  TEST_ASSERT(q->pop(e3));

  TEST_ASSERT(e1 == 1);
  TEST_ASSERT(e2 == 2);
  TEST_ASSERT(e3 == 3);

  TEST_ASSERT(q->isEmpty());

  delete (q);
}

void produce(ConcurrentQueue<int>& q, const size_t id, const int numToProduce) {
  for (int i = 0; i < numToProduce; ++i) {
    q.push(i);
  }
}

std::atomic_int at_sum(0);

void consume(ConcurrentQueue<int>& q, const size_t id) {
  int element;
  while (q.pop(element)) {
    at_sum += element;
  }
}

bool testProducerConsumer(const int numProducerThreads,
                          const int numConsumerThreads) {
  ConcurrentQueue<int> q(5);
  TEST_ASSERT(q.isEmpty());
  const int numToProduce = 10;
  const int expectedSum =
      numProducerThreads * (numToProduce - 1) * (numToProduce) / 2;
  std::vector<std::thread> producers(numProducerThreads);
  std::vector<std::thread> consumers(numConsumerThreads);

  //! start producer threads
  for (int i = 0; i < numProducerThreads; i++) {
    producers[i] = std::thread(produce, std::ref(q), i, numToProduce);
  }
  at_sum = 0;
  //! start consumer threads
  for (int i = 0; i < numConsumerThreads; i++) {
    consumers[i] = std::thread(consume, std::ref(q), i);
  }

  std::for_each(producers.begin(), producers.end(),
                std::mem_fn(&std::thread::join));
  //! the producer is done producing
  q.setDone();

  std::for_each(consumers.begin(), consumers.end(),
                std::mem_fn(&std::thread::join));
  TEST_ASSERT(q.isEmpty());

  return (at_sum == expectedSum);
}

void testMultipleTimes() {
  const int trials = 10000;
  //! Single Producer, Single Consumer
  for (int i = 0; i < trials; i++) {
    bool result = testProducerConsumer(1, 1);
    std::cerr
        << "\rTesting Single Producer, Single Consumer. Iterations Remaining: "
        << (trials - i - 1) << ' ' << std::flush;
    TEST_ASSERT(result);
  }
  std::cout << std::endl;

  //! Single Producer, Multiple Consumer
  for (int i = 0; i < trials; i++) {
    bool result = testProducerConsumer(1, 5);
    std::cerr << "\rTesting Single Producer, Multiple Consumer. Iterations "
                 "Remaining: "
              << (trials - i - 1) << ' ' << std::flush;
    TEST_ASSERT(result);
  }
  std::cout << std::endl;
  //! Multiple Producer, Single Consumer
  for (int i = 0; i < trials; i++) {
    bool result = testProducerConsumer(5, 1);
    std::cerr << "\rTesting Multiple Producer, Single Consumer. Iterations "
                 "Remaining: "
              << (trials - i - 1) << ' ' << std::flush;
    TEST_ASSERT(result);
  }
  std::cout << std::endl;

  //! Multiple Producer, Multiple Consumer
  for (int i = 0; i < trials; i++) {
    bool result = testProducerConsumer(2, 4);
    std::cerr << "\rTesting Multiple Producer, Mutiple Consumer. Iterations "
                 "Remaining: "
              << (trials - i - 1) << ' ' << std::flush;
    TEST_ASSERT(result);
  }
  std::cout << std::endl;
}

int main() {
  RDLog::InitLogs();

  //! test basic push and pop operations
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testPushAndPop();
  BOOST_LOG(rdErrorLog) << "Finished: testPushAndPop() \n";
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
#ifdef RDK_TEST_MULTITHREADED
  //! multiple producers and multiple consumers
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testMultipleTimes();
  BOOST_LOG(rdErrorLog) << "Finished: testMultipleTimes() \n";
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
#endif
  return 0;
}

#endif
