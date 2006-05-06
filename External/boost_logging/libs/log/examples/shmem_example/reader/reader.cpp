// reader.cpp : Defines the entry point for the console application.
//

#include <boost/shmem/shmem_named_shared_object.hpp>
#include <iostream>
#include <string>
#include <fstream>
using namespace boost::shmem;

typedef char char_t;

int main() {
    named_shared_object segment;
    char_t * memory = 0;
    long * occupied_size = 0;

    const int just_in_case = 8192;
    std::size_t mem_size = 2 * 1024 * 1024 ;
    segment.open_or_create("/shared_test", mem_size + just_in_case);

    // the string
    { typedef std::pair<char_t*, std::size_t> pair;
    pair res = segment.find<char_t>("shared_log_object");
    if ( !res.first)
        // we're creating it right now
        segment.construct<char_t>("shared_log_object")[mem_size](0);

    res = segment.find<char_t>("shared_log_object");
    assert( res.first); // should be created by now
    memory = res.first;
    }

    // the occupied size
    { typedef std::pair<long*, std::size_t> pair;
    pair res = segment.find<long>("shared_occupied_size");
    if ( !res.first) 
        // we're creating it right now
        segment.construct<long>("shared_occupied_size")[1](0);

    res = segment.find<long>("shared_occupied_size");
    assert( res.first); // should be created by now
    occupied_size = res.first;
    }

    std::cout << "Shared Memory Logger Object created. \n"
        "Press a key when you want to copy the contents of the log (to log.txt file)." << std::endl;
    std::cin.get();

    std::ofstream out("log.txt");
    out << std::string(memory, memory + *occupied_size);
	return 0;
}

