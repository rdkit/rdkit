#ifdef RDK_THREADSAFE_SSS
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <functional>

#include "ProducerConsumer.h"

using namespace RDKit;

bool oneTest() {
  const int producerThreads = 2;
  const int consumerThreads = 5;
  const int expectedSum = producerThreads * 10;
  std::vector<std::thread> producers(producerThreads);
  std::vector<std::thread> consumers(consumerThreads);
  ProducerConsumer<int> obj;
  size_t sum = 0;
  for (size_t i = 0; i < producerThreads; ++i) {
    producers[i] =
        std::thread(&ProducerConsumer<int>::produce, std::ref(obj), 5, 3);
  }

  for (size_t i = 0; i < consumerThreads; ++i) {
    consumers[i] = std::thread(&ProducerConsumer<int>::consume, std::ref(obj),
                               std::ref(sum));
  }
  std::for_each(producers.begin(), producers.end(),
                std::mem_fn(&std::thread::join));
  obj.setDone();
  std::for_each(consumers.begin(), consumers.end(),
                std::mem_fn(&std::thread::join));
  return (sum == expectedSum);
}

void concurrentTest() {
  const int runs = 100;
  bool result = true;
  for (int i = 0; i < runs; i++) {
    bool curResult = oneTest();
    result = result & curResult;
    if (!result) break;
  }
  TEST_ASSERT(result);
}

int main() {
  RDLog::InitLogs();

#ifdef RDK_TEST_MULTITHREADED
  //! multiple producers and multiple consumers
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  concurrentTest();
  BOOST_LOG(rdErrorLog) << "Finished: concurrentTest() \n";
  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
#endif

  return 0;
}

#endif
