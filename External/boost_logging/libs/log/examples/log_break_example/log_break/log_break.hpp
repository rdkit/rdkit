// Log Break
//
// Copyright (C) 2004 John Torjo (john@torjo.com)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
//
//
// You can find the latest version of this application at http://www.torjo.com/logbreak/

#ifndef JT_LOGBREAK_BREAKPOINT_HPP
#define JT_LOGBREAK_BREAKPOINT_HPP

#pragma once

#include <iostream>
#include <fstream>
#include <windows.h>
#include <assert.h>
#include <vector>
#include <string>
#include <time.h>

namespace {
    const char LOGBREAK_SIGNATURE[] = "LogBreak_v1.0 - see http://www.torjo.com/logbreak/ for details";
}

typedef enum breakpoint_type {
    contains_str,
    does_not_contain_str,
    starts_with_str,
    ends_with_str
};

/*
    Represents a possible breakpoint
*/
struct breakpoint {
    breakpoint() : match(contains_str), break_after(0), break_every(1) {}

    // the condition we're monitoring
    std::string condition;

    // type of matching
    breakpoint_type match;

    // when the breakpoint is hit (condition is met), how many hits
    // before breaking into Debug mode?
    //
    // if zero, as the first breakpoint is hit, it'll break.
    // if one, it will only break into Debug mode the second time this breakpoint is hit
    // etc.
    int break_after;

    // after breaking into Debug mode once, when should it break again?
    // if one, it'll break each time this breakpoint is hit
    // if two, it'll break each other time this breakpoint is hit
    // etc.
    int break_every;
};

inline bool operator<(const breakpoint & first, const breakpoint &second) {
    if ( first.condition != second.condition) return first.condition < second.condition;
    if ( first.match != second.match) return first.match < second.match;
    if ( first.break_after != second.break_after) return first.break_after < second.break_after;
    return first.break_every < second.break_every;
}


inline std::ostream& operator<<( std::ostream& out, const breakpoint & val) {
    // ... escape each quote by doubling it
    out << '"';
    for ( std::string::const_iterator begin = val.condition.begin(), end = val.condition.end(); begin != end; ++begin)
        if ( *begin != '"') out << *begin;
        else out << "\"\"";

    out << '"' << ' ' << val.match << ' ' << val.break_after << ' ' << val.break_every;
    return out;
}

inline std::istream& operator>>( std::istream& in, breakpoint & val) {
    // unescape string
    std::string str;
    char ch;
    in >> ch;
    if (in)
        assert(ch == '"');
    while ( in) {
        ch = in.get(); 
        if ( ch != '"') str += ch;
        else {
            char next = in.get();
            if ( next == '"') str += '"';
            else 
                break; // read end of string
        }
    }

    val.condition = str;
    int match = 0;
    in >> match >> val.break_after >> val.break_every;
    val.match = (breakpoint_type)match;
    return in;
}



/*
    FIXME - make log_break_ext(), which re-reads the "-info" files once per second - in case we're updating it
    with the LogBreak program!
*/


/*
    Monitors when certain messages are written to the log.
    If a breakpoint matches, it will break into debug mode.

    See http://www.torjo.com/logbreak/ for more details
*/
struct log_break {
    log_break(const std::string & file_name, bool output_to_dbg_wnd = true, bool output_to_console = true) 
            : m_file_name(file_name + ".logbreak"), 
              output_to_dbg_wnd(output_to_dbg_wnd), output_to_console(output_to_console) {
        std::ifstream in( m_file_name.c_str() );
        std::string line;
        std::getline(in,line);
        assert( (line == LOGBREAK_SIGNATURE) || line.empty() ); // should be our own file ;)
        breakpoint bk;
        while ( in >> bk) m_hits.push_back(bk);
    }

private: struct hit_type;
public:
    static bool matches_condition(const std::string & msg, const hit_type & hit);

    void operator()(const std::string&, const std::string & msg) {
        bool matched_already = false;
        for ( hits_array::const_iterator begin = m_hits.begin(), end = m_hits.end(); begin != end; ++begin) 
            if ( matches_condition( msg, *begin) ) {
                LONG count = ::InterlockedIncrement( &(begin->hits_count));

                // write that we've hit this - useful for the LogBreak program
                { std::ofstream out((m_file_name + "-info").c_str(), std::ios::app);
                    out << time(0) << ' ' << begin->bk << std::endl; }

                if ( (count - begin->bk.break_after) % begin->bk.break_every == 0) {
                    std::string breakpoint_condition = begin->bk.condition;
                    // now that we've hit a breakpoint, even if multiple conditions match it, don't break
                    // into Debug multiple times - useless
                    if ( !matched_already) {
                        // watch the 'breakpoint_condition' to see what breakpoint got hit
                        std::string msg_no_enter = msg;
                        while ( !msg_no_enter.empty() && msg_no_enter.end()[-1] == '\n') 
                            msg_no_enter.erase( msg_no_enter.end() -1);
                        std::string user_friendly_msg = "***** Breakpoint Hit: message '" + msg_no_enter 
                            + "'\n                      "
                            + (begin->bk.match == contains_str ? "contains" : "does not contain")
                            + " '" + breakpoint_condition + "'\n";
                        if (output_to_dbg_wnd) ::OutputDebugString(user_friendly_msg.c_str());
                        if (output_to_console) std::cout << user_friendly_msg << std::endl;
                        ::DebugBreak();
                        matched_already = true;
                    }
                }
            }
    }

private:
    struct hit_type {
        hit_type(const breakpoint & bk_) : hits_count(0), bk(bk_), fixed_start(false), fixed_end(false) {
            switch ( bk.match) {
            case starts_with_str: bk.match = contains_str; bk.condition = "*" + bk.condition; break;
            case ends_with_str:   bk.match = contains_str; bk.condition += "*"; break;
            }

            std::string condition = bk.condition;
            std::string::size_type pos = 0, next;
            while ( pos < condition.size() ) {
                next = condition.find('*', pos);
                if ( next == std::string::npos) next = condition.size();
                if ( next - pos > 0)
                    searches.push_back( condition.substr(pos, next - pos) );
                pos = next + 1;
            }
            if ( searches.empty()) return; // string was like "****" (only wildcards)

            fixed_start = *condition.begin() != '*';
            fixed_end = *condition.rbegin() != '*';
        }
        // the breakpoint
        breakpoint bk;

        // optimized version of 'condition' - faster to see if a string matches the condition
        std::vector<std::string> searches;
        // if fixed, when search starts/ends, the string must start with searches.begin()/end with searches.rbegin()
        bool fixed_start, fixed_end;

        // how many times has this breakpoint been hit?
        mutable volatile LONG hits_count;
    };

    std::string m_file_name;
    typedef std::vector<hit_type> hits_array;
    hits_array m_hits;

    // In case the breakpoint gets hit, Should we output to the Debug Window and/or to console
    //
    // These are useful in case you want to easily see what breakpoint got hit
    bool output_to_dbg_wnd;
    bool output_to_console;
};




#endif

