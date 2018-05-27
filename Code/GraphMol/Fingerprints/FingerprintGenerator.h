#include <DataStructs/SparseIntVect.h>

#ifndef _RD_FINGERPRINTGEN_H_
#define _RD_FINGERPRINTGEN_H_

namespace RDKit {
    class ROMol;
    

    class AtomEnvironmentGenerator{};

    class AtomInvariants{
       std::vector<boost::uint32_t> getAtomInvariants(const ROMol &mol);
    };

    class BondInvariants{
        std::vector<boost::uint32_t> getBondInvariants(const ROMol &mol);
    };

    class FingerprintGenerator{
        AtomEnvironmentGenerator atomEnvironmentGenerator;
        AtomInvariants atomInvariant;
        BondInvariants BondInvariant;

        //extra return
        std::vector<std::vector<boost::uint32_t> > *atomBits;
        std::map<boost::uint64_t,std::vector<std::vector<int> > > *bitInfo;

        public:
        SparseIntVect<boost::int32_t> *getFingerPrint(const ROMol &mol);

        std::vector<std::vector<boost::uint32_t> > *getLastAtomBits();
        std::map<boost::uint64_t,std::vector<std::vector<int> > > * getLastBitInfo();

    }


}

#endif
