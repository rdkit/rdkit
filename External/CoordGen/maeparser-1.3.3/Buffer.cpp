#include <algorithm>
#include <cerrno>
#include <cstring>
#include <ios>
#include <stdexcept>

#include "Buffer.hpp"

namespace schrodinger
{

Buffer::Buffer(size_t buffer_size) : m_data(buffer_size)
{
    begin = current = m_data.begin();
    // Set end to begin since no data has been loaded.
    end = begin;
}

Buffer::Buffer(std::istream& stream, size_t buffer_size) : Buffer(buffer_size)
{
    m_loader = new StreamLoader(stream);
}

Buffer::Buffer(FILE* file, size_t buffer_size) : Buffer(buffer_size)
{
    m_loader = new FileLoader(file);
}

Buffer::Buffer(const std::string& str) : Buffer(str.size())
{
    const char* c = str.c_str();
    std::copy(c, c + str.size(), m_data.begin());
}

Buffer::~Buffer()
{
    if (m_loader != nullptr) {
        delete m_loader;
    }
}

size_t Buffer::getColumn(const char* ptr) const
{
    assert(ptr >= begin && ptr <= end);

    const char* save = ptr;
    while (ptr > begin) {
        --ptr;
        if (*ptr == '\n') {
            return save - ptr;
        }
    };
    return (save - ptr) + m_starting_column;
}

bool Buffer::load(char*& save)
{
    // begin, current, end are public member variables
    if (current < end) {
        return true;
    }

    if (m_loader == nullptr) {
        return false;
    }

    size_t new_size = m_data.size();
    if (new_size == 0) {
        new_size = m_loader->getDefaultSize();
    }

    size_t saved_chars = 0;
    if (save != nullptr) {
        saved_chars = end - save;
        if (saved_chars > new_size / 2) {
            new_size = saved_chars * 2;
        }
    }

    BufferData data(new_size);

    if (m_loader->load(data, save, end)) {
        m_starting_column = this->getColumn();
        // line_number stays the same
        m_data = data;
        begin = save = m_data.begin();
        current = begin + saved_chars;
        end = begin + m_data.size();
        return true;
    } else {
        // If nothing new got loaded, keep the Buffer unchanged and return
        // false
        return false;
    }
}

std::ostream& operator<<(std::ostream& os, const Buffer& b)
{
    size_t max_length = 10;
    size_t head_length = std::min(max_length, b.size());
    std::string head(b.begin, b.begin + head_length);
    os << "Buffer(" << head << "...)";
    return os;
}

BufferData::BufferData(size_t size) : m_size(size)
{
    // Allocate space and add a trailing null character.
    m_data.resize(m_size + 1);
    m_data[m_size] = '\0';
}

void BufferData::resize(size_t size)
{
    if (size >= m_data.size()) {
        throw std::runtime_error("BufferData size can't be increased.");
    }
    m_size = size;
    m_data[m_size + 1] = '\0';
}

bool BufferDataCollector::load(BufferData& data, const char* begin,
                               const char* end) const
{
    bool succeeded = m_loader->load(data, begin, end);
    if (succeeded) {
        m_tokens_list->appendBufferData(data);
    }
    return succeeded;
}

size_t BufferDataCollector::readData(char*, size_t) const
{
    // char* ptr, size_t size are unnamed to avoid compilation warning messages.

    // The BufferDataCollector doesn't actually read any data directly; it
    // delegates that to its member BufferLoader instance.
    return 0;
}

// C++ wart; class static const definitions. See Item 2 of Effective C++
// (3rd ed).
const size_t BufferLoader::DEFAULT_SIZE;

bool BufferLoader::load(BufferData& data, const char* begin,
                        const char* end) const
{
    ptrdiff_t saved_chars = 0;
    if (begin != nullptr && end != nullptr) {
        saved_chars = end - begin;
    }

    if (saved_chars && begin != data.begin()) {
        std::copy(begin, end, data.begin());
    }

    size_t bytes =
        readData(data.begin() + saved_chars, data.size() - saved_chars);
    if (bytes < (data.size() - saved_chars)) {
        data.resize(bytes + saved_chars);
    }

    return bytes > 0;
}

size_t FileLoader::readData(char* ptr, size_t size) const
{
    size_t bytes = fread(ptr, sizeof(char), size, m_file);
    if (bytes < size && ferror(m_file) != 0) {
        std::string err(strerror(errno));
        throw std::runtime_error("An error occurred: " + err);
    }
    return bytes;
}

size_t StreamLoader::readData(char* ptr, size_t size) const
{
    m_stream.read(ptr, size);
    if (m_stream) {
        return size;
    } else {
        if (m_stream.bad()) {
            throw std::runtime_error("Error in reading stream.");
        } else {
            return m_stream.gcount();
        }
    }
}

void TokenBufferList::appendBufferData(const BufferData& buffer_data)
{
    if (m_token_buffer_list.empty()) {
        m_token_buffer_list.emplace_back(buffer_data, 0);
        return;
    }

    TokenBuffer& previous_buffer = m_token_buffer_list.back();
    size_t next_index = m_begin.size();

    // If the previous buffer stored no values, throw it away.
    if (previous_buffer.first_value == previous_buffer.last_value) {
        m_token_buffer_list.pop_back();
    }

    m_token_buffer_list.emplace_back(buffer_data, next_index);
}

void TokenBufferList::getData(size_t index, const char** const data,
                              size_t* const length) const
{
    // If this isn't true it means that getData was called before data
    // collection was complete, which is an incorrect use of the API.
    assert(m_begin.size() == m_end.size());
    auto token_buffer_iter = m_token_buffer_list.begin();

    while (index >= token_buffer_iter->last_value) {
        ++token_buffer_iter;
        // This indicates a problem in the token buffer logic - that the
        // index being searched for was not part of the BufferStringList.
        assert(token_buffer_iter != m_token_buffer_list.end());
    }

    *length = m_end[index] - m_begin[index];
    *data = token_buffer_iter->buffer_data.begin() + m_begin[index];
}

} // namespace schrodinger
