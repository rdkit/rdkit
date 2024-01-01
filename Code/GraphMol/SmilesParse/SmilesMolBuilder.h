#pragma once

#include <memory>
#include <string_view>

namespace RDKit {
class RWMol;
class Atom;
class Bond;
}  // namespace RDKit

namespace smiles_parser {
namespace mol_builder {

// private apis to construct mol components
[[nodiscard]] std::unique_ptr<::RDKit::RWMol> SmilesToMol(
    std::string_view smiles);

[[nodiscard]] std::unique_ptr<::RDKit::Atom> SmilesToAtom(
    std::string_view atom_token);

[[nodiscard]] std::unique_ptr<::RDKit::Bond> SmilesToBond(
    std::string_view bond_token);
}  // namespace mol_builder
}  // namespace smiles_parser
