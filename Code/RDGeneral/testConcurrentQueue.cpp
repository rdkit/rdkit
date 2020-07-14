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
namespace io = boost::iostreams;
using namespace std::placeholders;

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
    // PrintThread{} << "   Producer " << id << " --> item " << std::setw(3) <<
    // i << '\n';
  }
  // PrintThread{} << "EXIT: Producer " << id << '\n';
}

void consume(ConcurrentQueue<int>& q, const size_t id, int& sum) {
  int element;
  while (q.pop(element)) {
    sum += element;
    // PrintThread{} << "                  item " << std::setw(3) << element  <<
    // " --> Consumer " << id << '\n';
  }
  // PrintThread{} << "EXIT: Consumer " << id << '\n';
}

bool testSingleProducerMultipleConsumers() {
  ConcurrentQueue<int> q(5);
  TEST_ASSERT(q.isEmpty());
  const int numProducerThreads = 2;
  const int numConsumerThreads = 3;
  const int numToProduce = 10;
  const int expectedSum =
      numProducerThreads * (numToProduce - 1) * (numToProduce) / 2;
  int sum = 0;

  std::vector<std::thread> producers(numProducerThreads);
  std::vector<std::thread> consumers(numConsumerThreads);

  //! start producer threads
  for (int i = 0; i < numProducerThreads; i++) {
    producers[i] = std::thread(produce, std::ref(q), i, numToProduce);
  }
  //! start consumer threads
  for (int i = 0; i < numConsumerThreads; i++) {
    consumers[i] = std::thread(consume, std::ref(q), i, std::ref(sum));
  }

  std::for_each(producers.begin(), producers.end(),
                std::mem_fn(&std::thread::join));
  //! the producer is done producing
  q.setDone();

  std::for_each(consumers.begin(), consumers.end(),
                std::mem_fn(&std::thread::join));
  TEST_ASSERT(q.isEmpty());

  return (sum == expectedSum);
}

void testMultipleTimes() {
  const int runs = 1000;
  for (int i = 0; i < runs; i++) {
    bool result = testSingleProducerMultipleConsumers();
    // std::cerr << "\rIterations remaining : " << (runs - 1 - i) << ' '  <<
    // std::flush;
    TEST_ASSERT(result == true);
  }
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
