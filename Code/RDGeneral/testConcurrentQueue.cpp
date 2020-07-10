#ifdef RDK_THREADSAFE_SSS
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <atomic>
#include <functional>
#include <vector>

#include "ConcurrentQueue.h"

using namespace RDKit;
namespace io = boost::iostreams;
using namespace std::placeholders;

const int numToProduce = 100;
const int numToConsume = 25;

void testPushAndPop() {
  ConcurrentQueue<int>* q = new ConcurrentQueue<int>(4);
  q->push(1);
  q->push(2);
  q->push(3);
  q->push(4);
  TEST_ASSERT(q->pop() == 1);
  TEST_ASSERT(q->pop() == 2);
  TEST_ASSERT(q->pop() == 3);
  TEST_ASSERT(q->pop() == 4);
  delete (q);
}

void produce(ConcurrentQueue<int>& q, int& producer_count) {
  for (int i = 0; i < numToProduce; ++i) {
    ++producer_count;
    q.push(producer_count);
  }
}

void consume(ConcurrentQueue<int>& q, int& consumer_count, int& sum) {
  for (int i = 0; i < numToConsume; ++i) {
    auto item = q.pop();
    ++consumer_count;
    sum += item;
  }
}
bool testSingleProducerMultipleConsumers() {
  ConcurrentQueue<int> q(10);
  const int numProducerThreads = 1;
  const int numConsumerThreads =
      (numToProduce * numProducerThreads) / numToConsume;

	int sum = 0;
	int producer_count = 0;
	int consumer_count = 0;

  std::vector<std::thread> producers(numProducerThreads);
  std::vector<std::thread> consumers(numConsumerThreads);

  //! start producer threads
  for (int i = 0; i < numProducerThreads; i++) {
    producers[i] = std::thread(std::bind(produce, std::ref(q), std::ref(producer_count)));
  }
  //! start consumer threads
  for (int i = 0; i < numConsumerThreads; i++) {
    consumers[i] = std::thread(std::bind(consume, std::ref(q), std::ref(consumer_count), std::ref(sum)));
  }

  std::for_each(producers.begin(), producers.end(),
                std::mem_fn(&std::thread::join));
  std::for_each(consumers.begin(), consumers.end(),
                std::mem_fn(&std::thread::join));

  int expectedSum = (numToProduce * numProducerThreads) *
                    (numToProduce * numProducerThreads + 1) / 2;
  int expectedProducerCount = numToProduce * numProducerThreads;
  int expectedConsumerCount = numToConsume * numConsumerThreads;

	return (sum == expectedSum) && (producer_count == expectedProducerCount) && (consumer_count == expectedConsumerCount);
}

void testMultipleTimes(){
	const int runs = 100000;
	for(int i = 0; i < runs; i++){
		bool result = testSingleProducerMultipleConsumers();
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
