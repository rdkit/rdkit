#pragma once

#include <cstdio>
#include <memory>
#include <string>

#include "MaeParserConfig.hpp"

namespace schrodinger
{
namespace mae
{

class Block;

class EXPORT_MAEPARSER Writer
{
  private:
    std::shared_ptr<std::ostream> m_out = nullptr;

    void write_opening_block();

  public:
    Writer() = delete;
    explicit Writer(const std::string& fname);
    Writer(std::shared_ptr<std::ostream> stream);

    void write(const std::shared_ptr<Block>& block);
};

} // namespace mae
} // namespace schrodinger
