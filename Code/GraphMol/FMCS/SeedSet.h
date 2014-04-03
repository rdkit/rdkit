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
        unsigned            MaxBonds;
        unsigned            MaxAtoms;

        typedef ValueSet::iterator       iterator;
        typedef ValueSet::const_iterator const_iterator;

    public:
        SeedSet() : MaxBonds(0), MaxAtoms(0) {}
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

//        Value& push_back(const Value& seed) { Seeds.push_back(seed); return Seeds.back();} // initial fast fill with smallest seeds have equal size
        Value& add(const Value& seed)
        {
            iterator where;
            for(where = Seeds.begin(); where != Seeds.end(); where++) // find position in sorted list
                if(where->getNumBonds() < seed.getNumBonds())
                    break;
            iterator it = Seeds.insert(where, seed);  // swap seed.{arrays} would be faster then copy them to new seed
//            iterator it = Seeds.insert(where, Seed());
            Value& val = *it;
//            val.swap((Seed*)&seed); // for insert EMPTY Seed !!!

            // update best sizes (MaxBonds, MaxAtoms)
            if(val.getNumAtoms() > MaxAtoms)
                MaxAtoms = val.getNumAtoms();
            if(val.getNumBonds() > MaxBonds)
            {
                MaxBonds = val.getNumBonds();
/* opt. OPTIMISATION OF MEMORY USAGE (and a bit performance) - do not necessary
                // remove all seed that are too small for future growing
                if(MaxBonds > 12)    // skip some initial steps for performance
                {
                    for(it++; it != Seeds.end(); )
                    {
                        if(!it->canGrowBiggerThan(MaxBonds, MaxAtoms))
                            it = Seeds.erase(it);
                        else
                            it++;
                    }
                }
*/
            }
            return val;
        }
    };

}}

