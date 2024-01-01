#pragma once

#include <cstdint>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

namespace smiles_parser {
namespace ast_parser {

struct [[nodiscard]] chirality_t {
  std::string_view tag;
  std::string_view cclass;
  int permutation;
};

struct [[nodiscard]] atom_t {
  std::string_view name;
  std::optional<chirality_t> chirality = std::nullopt;
  std::optional<int> map_number = std::nullopt;
  std::optional<int> isotope = std::nullopt;
  std::optional<int> formal_charge = std::nullopt;
  bool no_implicit_hs = false;
  int atomic_number = 0;
  std::optional<int> explicit_h_count = std::nullopt;
};

struct [[nodiscard]] bond_t {
  std::string_view token;
};

struct [[nodiscard]] ring_t {
  std::string_view token;
};

struct [[nodiscard]] branch_t {
  std::string_view token;
};

struct [[nodiscard]] dot_t {};

using event_t = std::variant<atom_t, bond_t, branch_t, ring_t, dot_t>;
using mol_events_t = std::vector<event_t>;

class [[nodiscard]] Builder {
 public:
  [[nodiscard]] Builder(std::string_view input) : d_input(input){};

  void add_atom(std::string_view atom_name);
  void open_branch();
  void close_branch();
  void add_dot();
  void add_bond(std::string_view bond_token);
  void add_ring(std::string_view ring_number);

  // atom event helpers
  void add_explicit_h(size_t count);
  void add_explicit_h(std::string_view token);
  void add_atom_charge(int atom_charge);
  void add_atom_charge(std::string_view token);
  void add_isotope_num(std::string_view token);
  void add_chirality_tag(std::string_view chirality_tag);
  void add_chirality_class(std::string_view chirality_class,
                           std::string_view chiral_permutation = {});
  void add_atom_map_number(std::string_view token);
  void set_no_implicit_hs();

  void save_error_message(size_t bad_token_idx, std::string_view error_message);
  void save_error_message(std::string_view bad_token, std::string_view err_msg);

  [[nodiscard]] std::optional<mol_events_t> result();

 private:
  mol_events_t d_events;
  std::string_view d_input;
  std::vector<std::string> d_error_messages;

  [[nodiscard]] atom_t& get_last_atom();
  [[nodiscard]] std::int32_t svtoi(std::string_view token);

  void validate_atoms();
  void validate_rings();
  void print_error_messages();
};

[[nodiscard]] std::optional<ast_parser::mol_events_t> parse(
    std::string_view val);

}  // namespace ast_parser
}  // namespace smiles_parser
