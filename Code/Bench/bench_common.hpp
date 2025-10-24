#pragma once

#include <vector>

#include <GraphMol/ROMol.h>

namespace bench_common {

constexpr const char *SAMPLES[] = {
    // tricky stereochem
    "O=C1N2[C@H]3N4CN5C(=O)N6C7C5N(C2)C(=O)N7CN2C(=O)N5CN1[C@@H]3N(CN1C5C2N(C1=O)C6)C4=O",

    // randomly selected
    "C[C@H](NS(/C=C/c1ccccc1)(=O)=O)C(OCC(N1CCC(C)CC1)=O)=O",
    "COc1ccc(-n2cc(CNC(C)c3ncncc3)c(-c3cc(F)ccc3)n2)cc1",
    "O=C1N=C([n+]2ccccc2)/C(=C/[O-])N1c1ccccc1",
    "c1csc(CNC(=O)C2CCC(CNS(c3ccccc3)(=O)=O)CC2)c1",
    "CC1=C2C(c3ccc(C)cc3)C3=C(N=C2NN1)CC(C)(C)CC3=O",
    "NC(N)=NN/C=C1\\C=CC=C([N+]([O-])=O)C1=O",
    "CCc1c2c(oc(=O)c1)c(CN1CCN(C)CC1)c(O)c(Cl)c2",
    "COc1ccc(C2CC(C(F)(F)F)N3NC=C(C(Nc4c(C)n(C)n(-c5ccccc5)c4=O)=O)C3=N2)cc1",
    "Cc1ccc(CC(=O)O/N=C(\\N)Cc2ccc([N+]([O-])=O)cc2)cc1",
    "Br.COc1ccc(/N=C/C=C2/OC(C)(C)OC(c3ccccc3)=C2)cc1",
};

std::vector<RDKit::ROMol> load_samples();

uint64_t nth_random(uint64_t n) noexcept;

}  // namespace bench_common
