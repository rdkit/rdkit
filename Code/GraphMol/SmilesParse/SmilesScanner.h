#pragma once

#include <istream>
#include <string_view>

#ifndef __FLEX_LEXER_H
#undef yyFlexLexer
#define yyFlexLexer yysmiles_FlexLexer
#include <FlexLexer.h>
#endif

namespace RDKit {
class SmilesParser;

class [[nodiscard]] SmilesScanner : public yysmiles_FlexLexer {
 public:
  [[nodiscard]] SmilesScanner(std::istream& input_stream, std::string_view ss)
      : yysmiles_FlexLexer(&input_stream), d_input(std::move(ss)) {}

  [[nodiscard]] int lex(SmilesParser::semantic_type* const lval,
                        SmilesParser::location_type* location,
                        int& start_token);

  [[nodiscard]] const std::string_view& input() { return d_input; }

 private:
  std::string_view d_input;
};
}  // namespace RDKit
