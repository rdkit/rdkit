#pragma once
#include <map>
#include <algorithm>
#include "Seed.h"

namespace RDKit
{
 namespace FMCS
 {
    class SeedSet // sorted by amount of bonds
    {
    protected:
        typedef std::list<Seed> ValueSet;
        ValueSet            Seeds;
    public:
        typedef Seed        Value;
        typedef ValueSet::iterator       iterator;
        typedef ValueSet::const_iterator const_iterator;

    public:
        void        clear() { Seeds.clear(); }
        void        erase(iterator where) { Seeds.erase(where); }
        size_t      size () { return Seeds.size (); }  // for statistics only
        bool        empty() { return Seeds.empty(); }
        iterator    begin() { return Seeds.begin(); }
        iterator    end  () { return Seeds.end  (); }
        Value&      front() { return Seeds.front(); }
        const_iterator begin()const { return Seeds.begin(); }
        const_iterator end  ()const { return Seeds.end  (); }
        const Value&   front()const { return Seeds.front(); }

        Value& push_back(const Value& seed) { Seeds.push_back(seed); return Seeds.back();}
        Value& add(const Value& seed)
        {
            iterator where;
            for(where = Seeds.begin(); where != Seeds.end(); where++) // find position in sorted list
                if(where->getNumBonds() < seed.getNumBonds())
                    break;
            iterator it = Seeds.insert(where, seed);
            Value& val = *it;
            return val;
        }
    };

}}

