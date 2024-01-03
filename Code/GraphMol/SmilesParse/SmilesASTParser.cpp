#include "SmilesASTParser.h"

#include <array>
#include <charconv>
#include <optional>
#include <RDGeneral/RDLog.h>
#include <sstream>
#include <unordered_map>
#include <utility>

#include "smiles.tab.hpp"
#include "SmilesScanner.h"

namespace smiles_parser {

namespace ast_parser {
int Builder::svtoi(std::string_view s) {
  std::int32_t result = 0;

  auto [unused, ec] = std::from_chars(s.data(), s.data() + s.size(), result);
  if (ec == std::errc::result_out_of_range) {
    save_error_message(s, "number too large");
  }

  return result;
}

void Builder::add_atom(std::string_view atom_name) {
  d_events.push_back(atom_t{atom_name});
}

void Builder::add_explicit_h(size_t count) {
  get_last_atom().explicit_h_count = count;
}

void Builder::add_explicit_h(std::string_view token) {
  get_last_atom().explicit_h_count = svtoi(token);
}
void Builder::add_isotope_num(std::string_view token) {
  get_last_atom().isotope = svtoi(token);
}

void Builder::add_chirality_tag(std::string_view chirality_tag) {
  get_last_atom().chirality = chirality_t{chirality_tag, "", 0};
}

void Builder::add_chirality_class(std::string_view chirality_class,
                                  std::string_view chiral_permutation) {
  get_last_atom().chirality =
      chirality_t{"",
                  chirality_class.substr(1),  // remove @ prefix
                  svtoi(chiral_permutation)};
}

void Builder::add_atom_charge(int atom_charge) {
  get_last_atom().formal_charge = atom_charge;
}

void Builder::add_atom_charge(std::string_view token) {
  get_last_atom().formal_charge = svtoi(token);
  if (d_input[std::distance(d_input.data(), token.data()) - 1] == '-') {
    *(get_last_atom().formal_charge) *= -1;
  }
}

void Builder::add_atom_map_number(std::string_view token) {
  get_last_atom().map_number = svtoi(token);
}

void Builder::add_bond(std::string_view bond_token) {
  d_events.push_back(bond_t{bond_token});
}

void Builder::add_ring(std::string_view ring_number) {
  auto prev_char =
      d_input[std::distance(d_input.data(), ring_number.data()) - 1];
  // if ring number is of the form %(N), use the number as is
  if (prev_char == '(') {
    d_events.push_back(ring_t{ring_number});
    return;
  }

  // if ring number is of the form %N, first 2 digits encode a single ring
  if (prev_char == '%') {
    d_events.push_back(ring_t{ring_number.substr(0, 2)});
  }

  // if there are additional ring numbers, skip the right number of digits
  size_t start_idx = prev_char == '%' ? 2 : 0;
  for (size_t i = start_idx; i < ring_number.size(); ++i) {
    d_events.push_back(ring_t{ring_number.substr(i, 1)});
  }
}

void Builder::open_branch() {
  std::string_view token = "(";
  d_events.push_back(branch_t{token});
}

void Builder::close_branch() {
  std::string_view token = ")";
  d_events.push_back(branch_t{token});
}

void Builder::set_no_implicit_hs() { get_last_atom().no_implicit_hs = true; }

void Builder::add_dot() { d_events.push_back(dot_t{}); }

atom_t& Builder::get_last_atom() {
  // assume last event is an atom
  return std::get<atom_t>(d_events.back());
}

static int get_atomic_number(std::string_view& atom_symbol) {
  static constexpr std::array<std::pair<std::string_view, int>, 138> info{{
      {"*", 0},     {"H", 1},     {"He", 2},    {"Li", 3},    {"Be", 4},
      {"B", 5},     {"b", 5},     {"C", 6},     {"c", 6},     {"N", 7},
      {"n", 7},     {"O", 8},     {"o", 8},     {"F", 9},     {"Ne", 10},
      {"Na", 11},   {"Mg", 12},   {"Al", 13},   {"Si", 14},   {"si", 14},
      {"P", 15},    {"p", 15},    {"S", 16},    {"s", 16},    {"Cl", 17},
      {"Ar", 18},   {"K", 19},    {"Ca", 20},   {"Sc", 21},   {"Ti", 22},
      {"V", 23},    {"Cr", 24},   {"Mn", 25},   {"Fe", 26},   {"Co", 27},
      {"Ni", 28},   {"Cu", 29},   {"Zn", 30},   {"Ga", 31},   {"Ge", 32},
      {"As", 33},   {"as", 33},   {"Se", 34},   {"se", 34},   {"Br", 35},
      {"Kr", 36},   {"Rb", 37},   {"Sr", 38},   {"Y", 39},    {"Zr", 40},
      {"Nb", 41},   {"Mo", 42},   {"Tc", 43},   {"Ru", 44},   {"Rh", 45},
      {"Pd", 46},   {"Ag", 47},   {"Cd", 48},   {"In", 49},   {"Sn", 50},
      {"Sb", 51},   {"Te", 52},   {"te", 52},   {"I", 53},    {"Xe", 54},
      {"Cs", 55},   {"Ba", 56},   {"La", 57},   {"Ce", 58},   {"Pr", 59},
      {"Nd", 60},   {"Pm", 61},   {"Sm", 62},   {"Eu", 63},   {"Gd", 64},
      {"Tb", 65},   {"Dy", 66},   {"Ho", 67},   {"Er", 68},   {"Tm", 69},
      {"Yb", 70},   {"Lu", 71},   {"Hf", 72},   {"Ta", 73},   {"W", 74},
      {"Re", 75},   {"Os", 76},   {"Ir", 77},   {"Pt", 78},   {"Au", 79},
      {"Hg", 80},   {"Tl", 81},   {"Pb", 82},   {"Bi", 83},   {"Po", 84},
      {"At", 85},   {"Rn", 86},   {"Fr", 87},   {"Ra", 88},   {"Ac", 89},
      {"Th", 90},   {"Pa", 91},   {"U", 92},    {"Np", 93},   {"Pu", 94},
      {"Am", 95},   {"Cm", 96},   {"Bk", 97},   {"Cf", 98},   {"Es", 99},
      {"Fm", 100},  {"Md", 101},  {"No", 102},  {"Lr", 103},  {"Rf", 104},
      {"Db", 105},  {"Sg", 106},  {"Bh", 107},  {"Hs", 108},  {"Mt", 109},
      {"Ds", 110},  {"Rg", 111},  {"Cn", 112},  {"Nh", 113},  {"Fl", 114},
      {"Mc", 115},  {"Lv", 116},  {"Ts", 117},  {"Og", 118},  {"Uun", 110},
      {"Uuu", 111}, {"Uub", 112}, {"Uut", 113}, {"Uuq", 114}, {"Uup", 115},
      {"Uuh", 116}, {"Uus", 117}, {"Uuo", 118},
  }};

  if (atom_symbol[0] == '\'') {
    atom_symbol = atom_symbol.substr(1, atom_symbol.size() - 2);
  }

  auto entry = std::find_if(info.begin(), info.end(), [&](auto target) {
    return target.first == atom_symbol;
  });
  return entry == info.end() ? std::numeric_limits<int>::max() : entry->second;
}

void Builder::validate_atoms() {
  for (auto& info : d_events) {
    if (!std::holds_alternative<atom_t>(info)) {
      continue;
    }

    auto& atom = std::get<atom_t>(info);
    if (std::isdigit(atom.name[0])) {
      atom.atomic_number = svtoi(atom.name);
    } else {
      atom.atomic_number = get_atomic_number(atom.name);
      if (atom.atomic_number == std::numeric_limits<int>::max()) {
        save_error_message(atom.name, "unsupported atom symbol");
      }
    }

    // hydrogen shouldn't have chirality
    if (atom.atomic_number == 1 && atom.chirality != std::nullopt) {
      save_error_message(atom.name, "chirality specified for hydrogen atom");
    }

    // source: http://opensmiles.org/opensmiles.html 3.1.3
    if (atom.formal_charge != std::nullopt &&
        (atom.formal_charge < -15 || atom.formal_charge > 15)) {
      save_error_message(atom.name, "charges should be -15 <= charge <= 15");
    }

    if (atom.chirality == std::nullopt || atom.chirality->cclass.empty()) {
      continue;
    }

    static constexpr std::string_view supported_classes{"TH AL SP TB OH"};
    auto cc = atom.chirality->cclass;
    if (supported_classes.find(cc) == std::string_view::npos) {
      save_error_message(cc, "unsupported chirality class");
    }
  }
}

void Builder::validate_rings() {
  for (auto& info : d_events) {
    if (!std::holds_alternative<ring_t>(info)) {
      continue;
    }

    auto& ring = std::get<ring_t>(info);
    auto prev_char =
        d_input[std::distance(d_input.data(), ring.token.data()) - 1];
    if (prev_char == '%' && ring.token.size() != 2) {
      save_error_message(ring.token,
                         "%N tokens must define two-digit ring numbers ");
    }

    // check whether number is out of bonds
    static constexpr int MAX_RING_NUMBER{99999};
    if (svtoi(ring.token) > MAX_RING_NUMBER) {
      save_error_message(ring.token, "ring number should be < 100,000");
    }
  }
}
void Builder::print_error_messages() {
  for (auto& error_message : d_error_messages) {
    BOOST_LOG(rdErrorLog) << error_message << "\n";
  }
}

std::optional<mol_events_t> Builder::result() {
  validate_atoms();
  validate_rings();
  print_error_messages();

  return d_error_messages.empty()
             ? std::make_optional<mol_events_t>(std::move(d_events))
             : std::nullopt;
}

void Builder::save_error_message(std::string_view bad_token,
                                 std::string_view err_msg) {
  save_error_message(std::distance(d_input.data(), bad_token.data()), err_msg);
}

void Builder::save_error_message(size_t bad_token_idx,
                                 std::string_view err_msg) {
  std::stringstream ss;
  ss << "SMILES Parse Error: check for mistakes around position "
     << bad_token_idx << "\n";

  // NOTE: If the input is very long, the pointer to the failed location
  // becomes less useful. We should truncate the length of the error message
  // to 101 chars.
  static constexpr unsigned int error_size{101};
  static constexpr unsigned int prefix_size{error_size / 2};
  static auto truncate_input = [](const auto& input, const unsigned int pos) {
    if ((pos >= prefix_size) && (pos + prefix_size) < input.size()) {
      return input.substr(pos - prefix_size, error_size);
    } else if (pos >= prefix_size) {
      return input.substr(pos - prefix_size);
    } else {
      return input.substr(
          0, std::min(input.size(), static_cast<size_t>(error_size)));
    }
  };

  size_t num_dashes =
      (bad_token_idx >= prefix_size ? prefix_size : bad_token_idx);

  ss << truncate_input(d_input, bad_token_idx) << "\n"
     << std::string(num_dashes, '-') << "^\n"
     << err_msg;
  d_error_messages.push_back(ss.str());
}

std::optional<mol_events_t> parse(std::string_view val) {
  std::istringstream iss(val.data());
  Scanner sc(iss, val);
  Builder ast_builder(val);

  static_cast<void>(Parser(sc, ast_builder).parse());
  return ast_builder.result();
}

}  // namespace ast_parser
}  // namespace smiles_parser
