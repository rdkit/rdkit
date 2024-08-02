#include "SmilesMolBuilder.h"

#include <array>
#include <charconv>
#include <unordered_map>

#include <GraphMol/RWMol.h>
#include <GraphMol/Bond.h>
#include <GraphMol/QueryBond.h>
#include "SmilesASTParser.h"
#include "SmilesParse.h"
#include "SmilesParseOps.h"

// helper class for std::visit
template <class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

namespace smiles_parser {

static std::unique_ptr<::RDKit::Atom> get_atom(
    const ast_parser::atom_t& atom_info);
static std::unique_ptr<::RDKit::Bond> get_bond(std::string_view bond_token);

namespace {

void add_atom_chirality(::RDKit::Atom& atom,
                        const ast_parser::chirality_t& chirality) {
  using ChiralType = ::RDKit::Atom::ChiralType;

  static constexpr std::string_view supported_tags{"@  @@ TH AL SP TB OH"};
  static constexpr std::array<ChiralType, 7> supported_chiralities{
      ChiralType::CHI_TETRAHEDRAL_CCW, ChiralType::CHI_TETRAHEDRAL_CW,
      ChiralType::CHI_TETRAHEDRAL,     ChiralType::CHI_ALLENE,
      ChiralType::CHI_SQUAREPLANAR,    ChiralType::CHI_TRIGONALBIPYRAMIDAL,
      ChiralType::CHI_OCTAHEDRAL};

  auto& tag = chirality.tag.empty() ? chirality.cclass : chirality.tag;
  atom.setChiralTag(supported_chiralities[supported_tags.find(tag) / 3]);
  if (chirality.cclass.empty()) {
    return;
  }

  atom.setProp(::RDKit::common_properties::_chiralPermutation,
               chirality.permutation);
}

void add_bond(::RDKit::RWMol& mol, int begin_atom, int end_atom,
              std::string_view bond_token, unsigned int num_bonds) {
  auto is_aromatic_bond = [&](auto atom1, auto atom2) {
    return mol.getAtomWithIdx(atom1)->getIsAromatic() &&
           mol.getAtomWithIdx(atom2)->getIsAromatic();
  };

  if (bond_token.empty()) {
    bond_token = is_aromatic_bond(begin_atom, end_atom) ? ":" : "-";
  }

  auto bond = get_bond(bond_token);

  using BondType = ::RDKit::Bond::BondType;
  auto bond_type = bond->getBondType();
  if (bond_type == BondType::DATIVER || bond_type == BondType::DATIVEL) {
    bond->setBondType(BondType::DATIVE);
  }

  if (bond_token == "<-") {
    bond->setBeginAtomIdx(end_atom);
    bond->setEndAtomIdx(begin_atom);
  } else {
    bond->setBeginAtomIdx(begin_atom);
    bond->setEndAtomIdx(end_atom);
  }

  bond->setProp("_cxsmilesBondIdx", num_bonds);
  mol.addBond(bond.release(), true);
}

size_t add_atom(::RDKit::RWMol& mol, const ast_parser::mol_events_t& mol_info,
                const ast_parser::event_t& atom_event, size_t prev_atom,
                unsigned int& num_bonds) {
  auto atom = get_atom(std::get<ast_parser::atom_t>(atom_event));
  if (mol.getNumAtoms() == 0) {
    atom->setProp(::RDKit::common_properties::_SmilesStart, 1);
    return mol.addAtom(atom.release(), true, true);
  }

  auto idx = mol.addAtom(atom.release(), true, true);
  auto& prev_token = mol_info[std::distance(mol_info.data(), &atom_event) - 1];
  std::visit(
      overloaded{
          [&](const ast_parser::dot_t&) {
            mol.getAtomWithIdx(idx)->setProp(
                RDKit::common_properties::_SmilesStart, 1, true);
          },
          [&](const ast_parser::bond_t& bond) {
            add_bond(mol, prev_atom, idx, bond.token, num_bonds++);
          },
          [&](const auto&) { add_bond(mol, prev_atom, idx, "", num_bonds++); },
      },
      prev_token);
  return idx;
}

void add_ring_bond(::RDKit::RWMol& mol,
                   const ast_parser::mol_events_t& mol_info,
                   const ast_parser::event_t& ring_event, size_t prev_atom,
                   unsigned int& num_bonds) {
  auto& prev_token = mol_info[std::distance(mol_info.data(), &ring_event) - 1];
  std::string_view bond_token = "UNSPECIFIED";
  if (std::holds_alternative<ast_parser::bond_t>(prev_token)) {
    bond_token = std::get<ast_parser::bond_t>(prev_token).token;
  }

  auto bond = get_bond(bond_token);
  auto atom = mol.getAtomWithIdx(prev_atom);

  auto ring_bond = mol.createPartialBond(prev_atom, bond->getBondType());
  if (bond->hasProp(RDKit::common_properties::_unspecifiedOrder)) {
    ring_bond->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
  }

  ring_bond->setBondDir(bond->getBondDir());
  int ring_number = [](auto s) {
    int result = 0;
    std::from_chars(s.data(), s.data() + s.size(), result);
    return result;
  }(std::get<ast_parser::ring_t>(ring_event).token);

  mol.setAtomBookmark(atom, ring_number);
  mol.setBondBookmark(ring_bond, ring_number);
  if (!(mol.getAllBondsWithBookmark(ring_number).size() % 2)) {
    ring_bond->setProp("_cxsmilesBondIdx", num_bonds++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom, &mol);

  std::vector<int> ring_closures;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,
                         ring_closures);
  ring_closures.push_back(-(ring_number + 1));
  atom->setProp(RDKit::common_properties::_RingClosures, ring_closures);
}

void construct_mol(::RDKit::RWMol& mol,
                   const ast_parser::mol_events_t& mol_info) {
  std::vector<int> branch_points;
  unsigned int num_bonds = 0;
  size_t prev_atom = 0;
  for (auto& info : mol_info) {
    if (std::holds_alternative<ast_parser::atom_t>(info)) {
      prev_atom = add_atom(mol, mol_info, info, prev_atom, num_bonds);
    } else if (std::holds_alternative<ast_parser::ring_t>(info)) {
      add_ring_bond(mol, mol_info, info, prev_atom, num_bonds);
    } else if (std::holds_alternative<ast_parser::branch_t>(info)) {
      auto& v = std::get<ast_parser::branch_t>(info);
      if (v.token == "(") {
        branch_points.push_back(prev_atom);
      } else {
        prev_atom = branch_points.back();
        branch_points.pop_back();
      }
    }
  }
}
}  // namespace

namespace mol_builder {
std::unique_ptr<::RDKit::RWMol> SmilesToMol(std::string_view smiles) {
  auto smiles_ast = ast_parser::parse(smiles);
  if (smiles_ast == std::nullopt) {
    return nullptr;
  }

  auto mol = std::make_unique<::RDKit::RWMol>();
  construct_mol(*mol, *smiles_ast);

  try {
    SmilesParseOps::CloseMolRings(mol.get(), false);
    SmilesParseOps::CheckChiralitySpecifications(mol.get(), true);
    SmilesParseOps::SetUnspecifiedBondTypes(mol.get());
    SmilesParseOps::AdjustAtomChiralityFlags(mol.get());
  } catch (const ::RDKit::SmilesParseException& e) {
    SmilesParseOps::CleanupAfterParseError(mol.get());
    mol.reset();
  }
  return mol;
}

std::unique_ptr<::RDKit::Atom> SmilesToAtom(std::string_view atom_token) {
  auto smiles_ast = ast_parser::parse(atom_token);
  if (smiles_ast == std::nullopt || smiles_ast->size() != 1) {
    return nullptr;
  }
  return get_atom(std::get<ast_parser::atom_t>(smiles_ast->front()));
}

std::unique_ptr<::RDKit::Bond> SmilesToBond(std::string_view bond_token) {
  return get_bond(bond_token);
}
}  // namespace mol_builder

static std::unique_ptr<::RDKit::Atom> get_atom(
    const ast_parser::atom_t& atom_info) {
  auto atom = std::make_unique<::RDKit::Atom>(atom_info.atomic_number);

  if (atom_info.isotope != std::nullopt) {
    atom->setIsotope(*atom_info.isotope);
  }

  if (atom_info.chirality != std::nullopt) {
    add_atom_chirality(*atom, *(atom_info.chirality));
  }

  if (atom_info.explicit_h_count != std::nullopt) {
    atom->setNumExplicitHs(*atom_info.explicit_h_count);
  }

  if (atom_info.formal_charge != std::nullopt) {
    atom->setFormalCharge(*atom_info.formal_charge);
  }

  atom->setIsAromatic(std::all_of(atom_info.name.begin(), atom_info.name.end(),
                                  [](auto c) { return std::islower(c); }));

  if (atom_info.name == "*") {
    atom->setProp(::RDKit::common_properties::dummyLabel, "*");
  }

  atom->setNoImplicit(atom_info.no_implicit_hs);
  if (atom_info.map_number != std::nullopt) {
    atom->setProp(RDKit::common_properties::molAtomMapNumber,
                  *atom_info.map_number);
  }

  return atom;
}

static std::unique_ptr<::RDKit::Bond> get_bond(std::string_view bond_token) {
  using BondType = ::RDKit::Bond::BondType;

  static constexpr std::string_view supported_tokens{
      "-  =  #  :  $  -> <- \\  \\\\ /"};
  static constexpr std::array<BondType, 10> supported_bonds{
      BondType::SINGLE,      BondType::DOUBLE,      BondType::TRIPLE,
      BondType::AROMATIC,    BondType::QUADRUPLE,   BondType::DATIVER,
      BondType::DATIVEL,     BondType::UNSPECIFIED, BondType::UNSPECIFIED,
      BondType::UNSPECIFIED,
  };

  std::unique_ptr<::RDKit::Bond> bond;
  if (auto idx = supported_tokens.find(bond_token);
      idx != std::string_view::npos && !bond_token.empty()) {
    bond.reset(new ::RDKit::Bond(supported_bonds[idx / 3]));
  } else if (bond_token == "~") {
    bond = std::make_unique<::RDKit::QueryBond>();
    bond->setQuery(::RDKit::makeBondNullQuery());
    // this is a special case for unspecified ring bond closure types.
    // we don't want to assign a default bond type at this point
  } else if (bond_token == "UNSPECIFIED") {
    bond.reset(new ::RDKit::Bond(BondType::UNSPECIFIED));
    bond->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
  }

  if (bond_token[0] == '\\') {
    bond->setBondDir(::RDKit::Bond::BondDir::ENDDOWNRIGHT);
    bond->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
  } else if (bond_token == "/") {
    bond->setBondDir(::RDKit::Bond::BondDir::ENDUPRIGHT);
    bond->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
  } else if (bond_token == ":") {
    bond->setIsAromatic(true);
  }

  return bond;
}

}  // namespace smiles_parser
