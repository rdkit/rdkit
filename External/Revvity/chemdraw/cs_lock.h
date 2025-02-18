// TEMPORARY PLACE HOLDER
#include <memory>
#include <mutex>

namespace cs {
  typedef std::mutex CriticalSection;
}

struct locker {
  std::mutex &mutex;
  locker(std::mutex &mutex): mutex(mutex) {
    mutex.lock();
  }
  ~locker() {
    mutex.unlock();
  }
};

#define PROTECT_GLOBAL_AND_STATIC_DATA(mutex) locker lock(mutex)
