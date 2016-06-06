#include "tags.h"
#include "Invariant.h"

using namespace RDKit;

template < typename T > std::string to_string( const T& n )
{
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}

const int MOD=10000;
const int COUNT=5000000;
static std::vector<std::string> nums;

void testTags() {
  for( int idx=0; idx<MOD; ++idx) {
    nums.push_back( to_string(idx) );
  }


  RDTags tags;
  {
    std::clock_t clock1 = std::clock();

    for( int idx=0; idx<COUNT; ++idx) {
      int i = idx % MOD;
      CHECK_INVARIANT(tags.get(nums[i]) == common_properties::MAX+i+1,
                      "Bad map");
    }
    std::clock_t clock2 = std::clock();
    
    std::cout << "tag gets: " << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;

  }
  {
    std::clock_t clock1 = std::clock();

    for( int idx=0; idx<COUNT; ++idx) {
      int i = idx % MOD;
      CHECK_INVARIANT(tags.get(nums[i]) == common_properties::MAX+i+1,
                      "Bad map");
    }
    std::clock_t clock2 = std::clock();
    
    std::cout << "tag gets: " << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;

  }
  
  for( int idx=0; idx<COUNT; ++idx) {
    int i = idx % MOD;
    CHECK_INVARIANT(nums[i] == tags.get(common_properties::MAX+i+1),
                    "Bad lookup");
  }
  
  }

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock_test(RDTags *tag, std::map<std::string, int> *m) {
  for( int idx=0; idx<COUNT; ++idx) {
    int i = idx % MOD;
    int value = tag->get(nums[i]);
    (*m)[nums[i]] = value;
  }
}

void runblock_time(RDTags *tag) {
  for( int idx=0; idx<COUNT; ++idx) {
    int i = idx % MOD;
    tag->get(nums[i]);
  }
}
}

#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread.hpp>
#include <RDGeneral/BoostEndInclude.h>
void testMultiThread() {
  RDTags tagmap;
  boost::thread_group tg;
  int count = 8;

  std::vector<std::map<std::string, int> > maps(count);
  for (int i=0; i<count; ++i) {
    tg.add_thread(new boost::thread(runblock_test, &tagmap, &maps[i]));
  }
  tg.join_all();

  for( int idx=0; idx<MOD; ++idx) {
    std::string k = to_string(idx);
    int i = maps[0][k];
    int i2 = maps[1][k];
    int i3 = maps[2][k];
    int i4 = maps[3][k];
    if (!(i == i2 == i3 == i4)) {
      //std::cout << k << ":" << i << " " << i2 << " " << i3 << " " << i4 << std::endl;
      CHECK_INVARIANT(i, "Bad insertion");
      CHECK_INVARIANT(i == maps[1][k], "Bad insertion 1");
      CHECK_INVARIANT(i == maps[2][k], "Bad insertion 1");
      CHECK_INVARIANT(i == maps[3][k], "Bad insertion 1");
    }
  }
}

void testMultiThreadTime() {
  RDTags tagmap;
  {
    boost::thread_group tg;
    int count = 4;
    std::clock_t clock1 = std::clock();
    
    for (int i=0; i<count; ++i) {
      tg.add_thread(new boost::thread(runblock_time, &tagmap));
    }
    tg.join_all();
    std::clock_t clock2 = std::clock();
    
    std::cout << "4 tagmaps colliding a lot: " << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }
  {
    boost::thread_group tg;
    int count = 4;
    std::clock_t clock1 = std::clock();
    
    for (int i=0; i<count; ++i) {
      tg.add_thread(new boost::thread(runblock_time, &tagmap));
    }
    tg.join_all();
    std::clock_t clock2 = std::clock();
    
    std::cout << "4 tagmaps colliding a lot (start with tls loaded): " << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }
  

}
#endif

int main() {
  std::cerr << "-- running tests -- " << std::endl;
  testTags();
#ifdef RDK_TEST_MULTITHREADED
  std::cerr << "-- running multithreaded -- " << std::endl;
  //testMultiThread();
  testMultiThreadTime();
#endif  
  return 0;
}
