#pragma once

#include <cstdio>
#include <memory>
#include <string>

#include "Buffer.hpp"
#include "MaeBlock.hpp"
#include "MaeParser.hpp"
#include "MaeParserConfig.hpp"

namespace schrodinger
{
namespace mae
{

class EXPORT_MAEPARSER Reader
{
  private:
    std::shared_ptr<MaeParser> m_mae_parser;

  public:
    Reader() = delete;
    Reader(FILE* file, size_t buffer_size = BufferLoader::DEFAULT_SIZE);

    Reader(const std::shared_ptr<std::istream>& stream,
           size_t buffer_size = BufferLoader::DEFAULT_SIZE);

    Reader(const std::string& fname,
           size_t buffer_size = BufferLoader::DEFAULT_SIZE);

    // Should be made private if we conclude there's no need for the
    // DirectParser. The only current purpose of allowing construction from a
    // MaeParser is to allow direct/buffered behavior difference.
    Reader(std::shared_ptr<MaeParser> mae_parser);

    std::shared_ptr<Block> next(const std::string& outer_block_name);
};

} // namespace mae
} // namespace schrodinger
