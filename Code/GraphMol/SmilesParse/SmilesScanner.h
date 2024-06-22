#pragma once

#include <istream>
#include <string_view>

#ifndef __FLEX_LEXER_H
#undef yyFlexLexer
#define yyFlexLexer yysmiles_FlexLexer
#include <FlexLexer.h>
#endif

#include "smiles.tab.hpp"

namespace smiles_parser {
namespace ast_parser {
class [[nodiscard]] Scanner : public yysmiles_FlexLexer {
 public:
  [[nodiscard]] Scanner(std::istream& input_stream, std::string_view ss)
      : yysmiles_FlexLexer(&input_stream), d_input(std::move(ss)) {}

  [[nodiscard]] int lex(Parser::semantic_type* const lval,
                        Parser::location_type* location);

  [[nodiscard]] const std::string_view& input() { return d_input; }

 private:
  std::string_view d_input;
};
}  // namespace ast_parser
}  // namespace smiles_parser
