#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>

#ifndef _RD_FINGERPRINTGEN_H_
#define _RD_FINGERPRINTGEN_H_

namespace RDKit {
    class ROMol;
    

    class AtomEnvironmentGenerator{};

    class AtomInvariants{
        // arguments

        public:
        std::vector<boost::uint32_t> getAtomInvariants(const ROMol &mol);
    };

    class BondInvariants{
        //arguments
        
        public:
        std::vector<boost::uint32_t> getBondInvariants(const ROMol &mol);
    };

    class FingerprintGenerator{
        AtomEnvironmentGenerator atomEnvironmentGenerator;
        AtomInvariants atomInvariant;
        BondInvariants bondInvariant;
        //more data to hold arguments and fp type

        //extra return
        std::vector<std::vector<boost::uint32_t> > *atomBits;
        std::map<boost::uint64_t,std::vector<std::vector<int> > > *bitInfo;

        public:
        SparseIntVect<boost::int32_t> *getFingerPrint(const ROMol &mol);
        SparseIntVect<boost::int32_t> *getCountFingerPrint(const ROMol &mol);
        SparseIntVect<boost::int32_t> *getHashedFingerPrint(const ROMol &mol);
        ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol);

        std::vector<std::vector<boost::uint32_t> > *getLastAtomBits();
        std::map<boost::uint64_t,std::vector<std::vector<int> > > *getLastBitInfo();

    }
}

#endif
