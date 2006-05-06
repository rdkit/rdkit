// use.cpp - write to the logs
#include "declare.h"

int main() {
    init_logs();

    /* 
        check out the out.txt.1, out.txt.2 , etc. files!
    */
    for ( int idx = 0; idx < 1500; ++idx)
        BOOST_LOG(app) << "this is the " << (idx+1) << " message. Each message on its own line\n";
}
