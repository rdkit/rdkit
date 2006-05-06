
#include "logs.h"
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
using namespace boost;

const int THREADS_COUNT = 20;

mutex m;
int remaining = THREADS_COUNT;

void writer_thread(int thread_id) {
    xtime t;
    xtime_get(&t, TIME_UTC);
    t.sec += 3;
    thread::sleep(t); // let all threads get created...

    for ( int idx = 0; idx < 1000; ++idx) {
        BOOST_LOG(app) << "application msg. at idx " << idx << " from thread " << thread_id << ' ' << std::endl;
        
        xtime t;
        xtime_get(&t, TIME_UTC);
        t.nsec += 1000;
        thread::sleep(t);

//        ::Sleep(5);

        BOOST_LOG(dbg) << "this is a dbg message (" << idx << ") from thread " << thread_id << ' ' << std::endl;
        BOOST_LOG(err) << "this is an error message (" << idx << ") from thread " << thread_id << ' ' << std::endl;
    }

    {
    mutex::scoped_lock lock(m);
    --remaining;
    }

    BOOST_LOG(app) << "EXITING thread " << thread_id << std::endl;
}

int main() {
    init_logs();
    for ( int idx = 0; idx < THREADS_COUNT; ++idx)
        thread( bind(writer_thread,idx) );

    xtime t;
    while (true) {
        xtime_get(&t, TIME_UTC);
        t.sec += 5;
        thread::sleep(t);

        mutex::scoped_lock lock(m);
        if (remaining <= 0) break;
    }    
}
