#ifndef _BUFFERED_READER_HPP
#define _BUFFERED_READER_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <iostream>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "MaeParserConfig.hpp"

namespace schrodinger
{

/**
 * A simple data class to hold unchanging character buffer data. Copies are
 * reference counted.
 */
class EXPORT_MAEPARSER BufferData
{
  private:
    std::vector<char> m_data;
    size_t m_size;

  public:
    explicit BufferData(size_t size = 0);

    /**
     * Return access to the beginning of the data buffer for loading.
     */
    char* begin() { return m_data.data(); }

    /**
     * Return a pointer to the beginning of the character buffer.
     */
    const char* begin() const { return m_data.data(); }

    /**
     * Return the logical size of the buffer.
     */
    size_t size() const { return m_size; }

    /**
     * Reduce (but not increase) logical size of the buffer (to keep track of
     * size).
     *
     * Specifying a size larger than the current one throws a runtime_error.
     *
     * This doesn't actually free up any memory or modify the underlying
     * character buffer.
     */
    void resize(size_t size);
};

/**
 * Base class for loading BufferData objects from some source.
 */
class BufferLoader
{
  private:
    size_t m_default_size;

  public:
    static const size_t DEFAULT_SIZE = 131072;

    explicit BufferLoader(size_t default_size = DEFAULT_SIZE)
        : m_default_size(default_size)
    {
    }

    virtual ~BufferLoader() = default;

    /**
     * Return the default buffer size for this BufferLoader. This should be
     * used to construct the BufferData object if no other size is preferred.
     */
    virtual size_t getDefaultSize() const { return m_default_size; }

    /**
     * Load the next chunk of data into the BufferData object.
     *
     * The amount of data loaded will be a maximum of the current size of
     * BufferData.
     *
     * Returns true if new characters were loaded.
     */
    bool load(BufferData& data) const
    {
        char* begin = nullptr;
        char* end = nullptr;
        return load(data, begin, end);
    }

    /**
     * Copy everything from 'begin' through 'end' to the beginning of the
     * BufferData object (to deal with tokens that span buffer boundaries) and
     * then load the next chunk of data into the remainder of the BufferData.
     *
     * Returns true if new characters were loaded.
     */
    virtual bool load(BufferData& data, const char* begin,
                      const char* end) const;

  protected:
    /**
     * Read 'size' bytes and dump them into 'ptr'. Return the number of bytes
     * loaded.
     */
    virtual size_t readData(char* ptr, size_t size) const = 0;
};

/**
 * A BufferLoader that reads from an input stream.
 *
 * Note that the input stream is not owned by the StreamLoader (and so
 * ifstreams must be closed by the caller).
 */
class StreamLoader : public BufferLoader
{
  private:
    std::istream& m_stream;

  public:
    StreamLoader(std::istream& stream) : m_stream(stream) {}

    StreamLoader() = delete;
    StreamLoader(const StreamLoader&) = delete;
    StreamLoader& operator=(const StreamLoader&) = delete;

    size_t readData(char* ptr, size_t size) const override;
};

/**
 * A BufferLoader that reads from a FILE pointer.
 *
 * Experiments on OS X show this to be measurably faster than input stream
 * reading. While it seems this isn't expected in theory, no way to make them
 * equivalent has yet been found in practice.
 *
 * Note that the FILE pointer is not owned by the FileLoader and must be closed
 * by the caller.
 */
class FileLoader : public BufferLoader
{
  private:
    FILE* m_file;

  public:
    FileLoader(FILE* file) : m_file(file) {}

    size_t readData(char* ptr, size_t size) const override;
};

/**
 * Character buffer.
 *
 * This is a shared resource manager for an allocated character buffer along
 * with iterator location information.
 *
 * The character buffer is stored in a BufferData object, which is reference
 * counted and can be retrieved by the data() method.
 */
class EXPORT_MAEPARSER Buffer
{
  private:
    BufferData m_data;
    BufferLoader* m_loader{nullptr};
    size_t m_starting_column{1};

  public:
    char* begin{nullptr};
    char* end{nullptr};
    char* current{nullptr};
    size_t line_number{1};

    explicit Buffer(size_t buffer_size = 0);

    /**
     * Construct an empty buffer that can be loaded from the provided input
     * stream.
     */
    explicit Buffer(std::istream& stream, size_t buffer_size = 0);

    /**
     * Construct an empty buffer that can be loaded from the provided FILE
     * pointer.
     */
    explicit Buffer(FILE* file, size_t buffer_size = 0);

    /**
     * Create a buffer from a string.
     *
     * This makes a copy of the string data. Calls to load() always return
     * false.
     */
    explicit Buffer(const std::string& str);

    ~Buffer();

    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

  public:
    void setBufferLoader(BufferLoader* loader) { m_loader = loader; }

    BufferLoader* getBufferLoader() { return m_loader; }

    /**
     * Load new BufferData from the BufferLoader. Update Buffer pointers.
     */
    bool load()
    {
        char* save = nullptr;
        return load(save);
    }

    /**
     * Save data from the 'save' pointer to the end of the current BufferData
     * into a new BufferData instance, then load the remainder of the
     * BufferData instance with data from the BufferLoader.
     *
     * This allows us to deal with tokens that cross buffer boundaries.
     *
     * Update Buffer pointers.
     */
    bool load(char*& save);

    inline size_t size() const { return m_data.size(); }

    BufferData data() const { return m_data; }

    inline bool operator==(const Buffer& other) const
    {
        return this->m_data.begin() == other.m_data.begin();
    }

    inline bool operator!=(const Buffer& other) const
    {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const Buffer& b);

    /**
     * Return the column number of the last character read.
     */
    size_t getColumn() const { return getColumn(current); }

    /**
     * Return the column number of the provided character.
     */
    size_t getColumn(const char* ptr) const;
};

/**
 * Allow Buffer objects to be written to a stream.
 *
 * This purpose of this is only to allow Boost testing assertions to work. It
 * doesn't do anything terribly useful.
 */
std::ostream& operator<<(std::ostream& os, const Buffer& b);

/**
 * A class to collect tokens with minimal copying by saving the buffer and
 * token start and end indices.
 */
class EXPORT_MAEPARSER TokenBufferList
{
  public:
    /// A simple data class to keep the info about the buffers and tokens
    // straight.
    class EXPORT_MAEPARSER TokenBuffer
    {
      public:
        BufferData buffer_data;
        /// The index of the first token stored in this buffer.
        size_t first_value;
        /// One greater than the index of the last token stored in this buffer.
        size_t last_value;

        TokenBuffer(BufferData data, size_t next_index)
            : buffer_data(std::move(data)), first_value(next_index),
              last_value(next_index)
        {
        }
    };

  private:
    /// List of TokenBuffer objects.
    std::list<TokenBuffer> m_token_buffer_list;

    /// Buffer indices for the beginnings of collected tokens.
    std::vector<size_t> m_begin;

    /// Buffer indices for one past the end of collected tokens.
    std::vector<size_t> m_end;

  public:
    TokenBufferList() : m_token_buffer_list(), m_begin(), m_end() {}

    void reserve(size_t size)
    {
        m_begin.reserve(size);
        m_end.reserve(size);
    }

    inline void setTokenIndices(size_t begin, size_t end)
    {
        m_begin.push_back(begin);
        m_end.push_back(end);
        m_token_buffer_list.back().last_value = m_end.size();
    }

    void appendBufferData(const BufferData& buffer_data);

    /**
     * Return token data as char pointer and length.
     *
     * Data is owned by the TokenBufferList and does not need to be freed.
     *
     * No trailing '\0' is present; data length must be observed.
     */
    void getData(size_t index, const char** const data,
                 size_t* const length) const;
};

/**
 * A class to modify a Buffer's loading behavior through RAII.
 *
 * It collects all loaded BufferData instances and stores them in a
 * TokenBufferList for parsing later.
 */
class BufferDataCollector : public BufferLoader
{
  private:
    Buffer* m_buffer;
    BufferLoader* m_loader;
    TokenBufferList* m_tokens_list;

  public:
    BufferDataCollector(Buffer* buffer, TokenBufferList* tokens_list)
        : m_loader(nullptr), m_tokens_list(tokens_list)
    {
        m_buffer = buffer;
        m_loader = m_buffer->getBufferLoader();
        m_buffer->setBufferLoader(this);
    }

    ~BufferDataCollector() override { m_buffer->setBufferLoader(m_loader); }

    bool load(BufferData& data, const char* begin,
              const char* end) const override;

    size_t readData(char* ptr, size_t size) const override;
};

} // end namespace schrodinger

#endif
