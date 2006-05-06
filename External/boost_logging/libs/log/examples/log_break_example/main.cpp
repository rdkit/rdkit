// logbreak_example.cpp : Defines the entry point for the console application.
//

#include <log_break/log_break.hpp>
#include "declare.h"

// FIXME case-sensitive !!!!

bool log_break::matches_condition(const std::string & msg, const hit_type & hit) {
    bool matches = hit.bk.match == contains_str;

    typedef std::string string;
    string::size_type at = 0;
    for ( std::vector<string>::const_iterator begin = hit.searches.begin(), end = hit.searches.end(); begin != end ; ++begin)  {
        string::size_type next = msg.find(*begin, at);
        if ( next == string::npos)
            return !matches; // doesn't match
        at = next + begin->length();
    }

    if ( hit.fixed_start)
        if ( msg.find(*hit.searches.begin()) != 0)
            return !matches; // prefix should be like "line*" (start with 'line')
    if ( hit.fixed_end)
        if ( at != msg.length() )
            return !matches; // prefix should be like "*line" (end in 'line')

    return matches;
}




int main()
{
    init_logs();

    std::ifstream in("example.txt");
    std::string line;
    while ( std::getline(in, line) ) {
        BOOST_LOG(app) << "Read line'" <<  line << "' - size " << line.length() << std::endl;
        if ( line.size() > 0 && toupper(line[0]) == line[0] )
            BOOST_LOG(dbg) << "line '" << line << "' starts in upper case" << std::endl;

        if ( line.size() > 70)
            BOOST_LOG(warn) << "line too long!" << std::endl;

        Sleep(100); // don't print too fast, let the user see something ;)
    }
	return 0;
}

