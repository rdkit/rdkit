#ifndef MAE_READER_HPP_
#define MAE_READER_HPP_

// Visual Studio versions prior to 2015 don't support noexcept
#if defined(_MSC_FULL_VER) && _MSC_FULL_VER < 190023026
#define NOEXCEPT
#else
#define NOEXCEPT noexcept
#endif

#include <cerrno>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>

#include <boost/dynamic_bitset.hpp>
#include <utility>

#include "Buffer.hpp"
#include "MaeBlock.hpp"
#include "MaeParserConfig.hpp"

#define MAEPARSER_EXCEPTION_BUFFER_SIZE 256

namespace schrodinger
{

namespace mae
{

/**
 * Parse (and throw away) a comment of the form '# comment #'.
 */
void comment(Buffer& buffer);

/**
 * Parse (and throw away) zero or more characters of whitespace including \t,
 * \r, \n and ' ', along with any embedded comments.
 */
EXPORT_MAEPARSER void whitespace(Buffer& buffer);

/**
 * Parse a triple colon or raise an exception.
 */
void triple_colon(Buffer& buffer);

/**
 * Parse the specific character requested. Return true if successful, false if
 * not.
 */
bool character(char c, Buffer& buffer);

/**
 * Parse the specific character requested. Return true if successful, false if
 * not. Update any save points if buffer reload is required.
 */
bool character(char c, Buffer& buffer, char*& save);

/**
 * Parse a full (b|i|r|s)_<author>_<name> property key.
 *
 * Return NULL if a starting character of ':' is found (the beginning of the
 * ':::' property name terminator).
 *
 * Raise a read_exception in any other situation.
 */
std::shared_ptr<std::string> property_key(Buffer& buffer);

/**
 * Read through the opening '{' of a named or unnamed outer block.
 *
 * Set the name of the block in the provided argument; an empty string
 * if unnamed. If EOF is reached, return false, otherwise return true.
 */
EXPORT_MAEPARSER std::string outer_block_beginning(Buffer& buffer);

template <typename T> T parse_value(Buffer& buffer);

class EXPORT_MAEPARSER read_exception : public std::exception
{
  private:
    char m_msg[MAEPARSER_EXCEPTION_BUFFER_SIZE];

    void format(size_t line_number, size_t column, const char* msg);

  public:
    read_exception(const Buffer& buffer, const char* msg)
    {
        format(buffer.line_number, buffer.getColumn(), msg);
    }

    read_exception(size_t line_number, size_t column, const char* msg)
    {
        format(line_number, column, msg);
    }

    const char* what() const NOEXCEPT override { return m_msg; }
};

/**
 * A pure virtual base class for parsers. Allows us to store these in
 * collections, so we can do things like easily invoke different parsers on
 * each element while parsing a row of values.
 */
class EXPORT_MAEPARSER Parser
{
  public:
    virtual void parse(Buffer& buffer) = 0;
    virtual ~Parser() = default;
};

class EXPORT_MAEPARSER IndexedBlockParser
{
    std::vector<std::string> m_property_names;

  public:
    virtual ~IndexedBlockParser() = default;

    virtual void parse(const std::string& name, size_t size,
                       Buffer& buffer) = 0;

    virtual std::shared_ptr<IndexedBlockMapI> getIndexedBlockMap() = 0;
};

class EXPORT_MAEPARSER IndexedBlockBuffer
{
  private:
    std::vector<std::string> m_property_names;
    std::string m_name;
    TokenBufferList m_tokens_list;
    size_t m_rows;

  public:
    IndexedBlockBuffer(std::string name, size_t rows)
        : m_property_names(), m_name(std::move(name)), m_rows(rows)
    {
    }

    virtual ~IndexedBlockBuffer() = default;

    void addPropertyName(std::string&& name)
    {
        m_property_names.push_back(name);
    }

    /**
     * Parse the indexed block values, store them in a linked list of buffers.
     */
    virtual void parse(Buffer& buffer);

    /**
     * Parse a value.
     */
    void value(Buffer& buffer);

    void getData(size_t ix, const char** data, size_t* const len) const
    {
        m_tokens_list.getData(ix, data, len);
    }

    std::string getName() const { return m_name; }

    size_t size() const { return m_rows; }

    IndexedBlock* getIndexedBlock();
};

class EXPORT_MAEPARSER BufferedIndexedBlockParser : public IndexedBlockParser
{
  private:
    std::shared_ptr<BufferedIndexedBlockMap> m_indexed_block_map;

  public:
    BufferedIndexedBlockParser();

    std::shared_ptr<IndexedBlockMapI> getIndexedBlockMap() override;

    void parse(const std::string& name, size_t size, Buffer& buffer) override;
};

class EXPORT_MAEPARSER DirectIndexedBlockParser : public IndexedBlockParser
{
    std::shared_ptr<IndexedBlockMap> m_indexed_block_map;

  public:
    void parse(const std::string& name, size_t size, Buffer& buffer) override;

    std::shared_ptr<IndexedBlockMapI> getIndexedBlockMap() override;
};

class EXPORT_MAEPARSER IndexedValueParser : public Parser
{
  public:
    virtual void addToIndexedBlock(IndexedBlock* block) = 0;
};

template <typename T> class IndexedValueCollector : public IndexedValueParser
{
  public:
    std::string m_name;
    std::vector<T> m_values;
    boost::dynamic_bitset<>* m_is_null;

  public:
    explicit IndexedValueCollector(std::string name, size_t size)
        : m_name(std::move(name)), m_values(), m_is_null(nullptr)
    {
        m_values.reserve(size);
        m_is_null = nullptr;
    }

    ~IndexedValueCollector() override
    {
        if (m_is_null) {
            delete m_is_null;
        }
    }

    void parse(Buffer& buffer) override
    {
        if (buffer.current >= buffer.end) {
            if (!buffer.load()) {
                throw read_exception(buffer, "Unexpected EOF.");
            }
        }
        if (*buffer.current == '<') {
            char* save = buffer.current;
            ++buffer.current;
            if (buffer.current >= buffer.end) {
                if (!buffer.load(save)) {
                    throw read_exception(buffer, "Unexpected EOF.");
                }
            }
            // TODO: not sure, but I assume that unquoted strings like
            // <foo> are allowed as values. Ugh. This requires saving the
            // starting point of '<' in case we need to back up.
            if (*buffer.current != '>') {
                // Back up and parse as a normal value.
                --buffer.current;
            } else {
                ++buffer.current;
                if (m_is_null == nullptr) {
                    m_is_null =
                        new boost::dynamic_bitset<>(m_values.capacity());
                }
                m_is_null->set(m_values.size());
                m_values.push_back(T());
                return;
            }
        }
        m_values.push_back(parse_value<T>(buffer));
    }

    void addToIndexedBlock(IndexedBlock* block) override
    {
        auto ptr = std::shared_ptr<IndexedProperty<T>>(
            new IndexedProperty<T>(m_values, m_is_null));
        block->setProperty<T>(m_name, ptr);
        m_is_null = nullptr;
    }
};

class EXPORT_MAEPARSER MaeParser
{
  protected:
    Buffer m_buffer;
    std::shared_ptr<std::istream> m_stream;

    virtual IndexedBlockParser* getIndexedBlockParser()
    {
        return new BufferedIndexedBlockParser();
    }

  public:
    explicit MaeParser(const std::shared_ptr<std::istream>& stream,
                       size_t buffer_size = BufferLoader::DEFAULT_SIZE)
        : m_buffer(*stream, buffer_size), m_stream(stream)
    {
        m_buffer.load();
    }

    explicit MaeParser(FILE* file,
                       size_t buffer_size = BufferLoader::DEFAULT_SIZE)
        : m_buffer(file, buffer_size)
    {
        if (file == nullptr) {
            std::string msg("Bad file argument");
            if (errno != 0) {
                msg += ": ";
                msg += strerror(errno);
            } else {
                msg += ".";
            }
            throw std::runtime_error(msg);
        }
        m_buffer.load();
    }

    // TODO: finish big three (four)
    virtual ~MaeParser() = default;

    std::shared_ptr<Block> blockBody(const std::string& name);

    IndexedBlock* indexedBlock(const std::string& name, size_t size);

    // TODO: private
    IndexedBlockBuffer* indexedBlockBuffer(const std::string& name,
                                           size_t size);

    std::shared_ptr<Block> outerBlock();

    /**
     * Read a block name or a closing '}'. The argument 'indexed' is set to
     * a positive integer value indicating the number of rows, or zero if
     * the block is not indexed.
     *
     * Return the block name or NULL if the closing '}' was found.
     */
    std::string blockBeginning(int* indexed);

    std::shared_ptr<std::string> property();

    /**
     * Read a list of properties, ending at the name/value separator.
     * Populate the provided std::vector of std::string shared pointers.
     */
    void properties(std::vector<std::shared_ptr<std::string>>* property_names);

    /**
     * Read (and throw away) any whitespace.
     */
    void whitespace() { schrodinger::mae::whitespace(m_buffer); }
};

class EXPORT_MAEPARSER DirectMaeParser : public MaeParser
{
  public:
    explicit DirectMaeParser(const std::shared_ptr<std::istream>& stream,
                             size_t buffer_size = BufferLoader::DEFAULT_SIZE)
        : MaeParser(stream, buffer_size)
    {
    }

    explicit DirectMaeParser(FILE* file,
                             size_t buffer_size = BufferLoader::DEFAULT_SIZE)
        : MaeParser(file, buffer_size)
    {
    }

  private:
    IndexedBlockParser* getIndexedBlockParser() override
    {
        return new DirectIndexedBlockParser();
    }
};

} // end namespace mae

} // end namespace schrodinger

#endif
