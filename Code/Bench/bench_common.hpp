#pragma once

#include <initializer_list>
#include <string_view>

namespace bench_common {
constexpr std::initializer_list<std::string_view> CASES = {
    "COC1/C=C/OC2(C)Oc3c(C)c(O)c4c(O)c(c(/C=N/OC(c5ccccc5)c5ccccc5)cc4c3C2=O)NC(=O)/C(C)=C\\C=C\\C(C)C(O)C(C)C(O)C(C)C(OC(C)=O)C1C",
    "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
    "c1ccc2c(c1)c3ccccc3c4ccccc24",
    "F[C@@H](Cl)[C@H](Br)[C@](I)(O)[C@@](N)(C#N)C(=O)O",
    "C[S+](C)(C)[O-]",
};
}  // namespace bench_common
