#pragma once

#include <istream>
#include <string_view>

#undef yyFlexLexer
#define yyFlexLexer smilesFlexLexer
#include <flex/FlexLexer.h>

#include "smiles.tab.hpp"

namespace smiles_parser {
namespace ast_parser {
class [[nodiscard]] Scanner : public smilesFlexLexer {
 public:
  [[nodiscard]] Scanner(std::istream& input_stream, std::string_view ss)
      : smilesFlexLexer(&input_stream), d_input(ss) {}

  [[nodiscard]] int lex(Parser::semantic_type* const lval,
                        Parser::location_type* location);

  [[nodiscard]] std::string_view input() { return d_input; }

 private:
  std::string_view d_input;
};
}  // namespace ast_parser
}  // namespace smiles_parser

// don't know why these aren't generated
inline int smilesFlexLexer::yylex() { return 1; }
inline int smilesFlexLexer::yywrap() { return 1; }
