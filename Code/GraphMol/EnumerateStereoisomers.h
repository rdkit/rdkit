// RDKit includes
#include <RDGeneral/export.h>
#include "Atom.h"
#include "Bond.h"
#include "StereoGroup.h"

// std includes
#include <tuple>
#include <vector>
#include <iostream>

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

    class _Flipper {
        public:
            virtual void flip(bool flag);
            virtual ~_Flipper() {};
    };

    class _BondFlipper : public _Flipper {
        public:
            _BondFlipper(Bond *bond);
            void flip(bool flag);
        private:
            Bond* bond;
    };

    class _AtomFlipper : public _Flipper {
        public:
            _AtomFlipper(Atom *atom);
            void flip(bool flag);
        private:
            Atom* atom;
    };

    class _StereoGroupFlipper : public _Flipper {
        public:
            _StereoGroupFlipper(RDKit::StereoGroup* group);
            void flip(bool flag);
        private:
            std::vector<std::tuple<Atom*, RDKit::Atom::ChiralType> > _original_parities;
    };

    std::vector<_Flipper*> _get_flippers(ROMol* mol, const StereoEnumerationOptions& options);

    class _RangeBitsGenerator;

    class _UniqueRandomBitsGenerator;
}

#endif // RD_ENUMSTEREO_H