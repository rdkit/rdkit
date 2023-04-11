#include <RDGeneral/export.h>
#include "Atom.h"
#include "Bond.h"
#ifndef RD_ENUMSTEREO_H
#define RD_ENUMSTEREO_H

namespace RDKit {

    class StereoEnumerationOptions {
        public:
            bool try_embedding, only_unassigned, only_stereo_groups, unique;
            unsigned int max_isomers, rand;
            StereoEnumerationOptions(
                bool try_embedding=false, bool only_unassigned=true,
                bool only_stereo_groups=false, bool unique=true,
                unsigned int max_isomers=1024, unsigned int rand=0);
            ~StereoEnumerationOptions() {};
    };

    class _BondFlipper {
        public:
            _BondFlipper(Bond *bond);
            ~_BondFlipper();
            void flip(bool flag);
        private:
            Bond* bond;
    };

    class _AtomFlipper {
        public:
            _AtomFlipper(Atom *atom);
            ~_AtomFlipper();
            void flip(bool flag);
        private:
            Atom* atom;
    };

    class _StereoGroupFlipper;

    class _RangeBitsGenerator;

    class _UniqueRandomBitsGenerator;
}

#endif // RD_ENUMSTEREO_H