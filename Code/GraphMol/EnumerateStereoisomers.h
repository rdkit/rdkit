// RDKit includes
#include "Atom.h"
#include "Bond.h"
#include "ROMol.h"
#include "StereoGroup.h"
#include <RDGeneral/export.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

// std includes
#include <tuple>
#include <cmath>
#include <vector>
#include <utility>
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

    std::vector<_Flipper*> _get_flippers(ROMol* mol, const StereoEnumerationOptions options=StereoEnumerationOptions());

    class _RangeBitsGenerator;

    class _UniqueRandomBitsGenerator;

    unsigned int get_stereoisomer_count(ROMol* mol, const StereoEnumerationOptions options=StereoEnumerationOptions());

    std::vector<ROMol*> enumerate_stereoisomers(ROMol* mol, const StereoEnumerationOptions options=StereoEnumerationOptions(), bool verbose=false);

}

#endif // RD_ENUMSTEREO_H