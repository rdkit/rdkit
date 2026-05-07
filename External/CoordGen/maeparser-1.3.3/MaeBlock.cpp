#include "MaeBlock.hpp"
#include <cmath>
#include <utility>

#include "MaeParser.hpp"

using namespace std;

namespace schrodinger
{
namespace mae
{

namespace
{
const double tolerance = 0.00001; // Tolerance to match string cutoff

// Wrap to-string to allow it to take strings and be a no-op
template <typename T> inline string local_to_string(T val)
{
    return to_string(val);
}

inline bool char_requires_escaping(char c)
{
    return c == '"' || c == '\\';
}

string local_to_string(const string& val)
{
    if (val.empty()) {
        return R"("")";
    }

    // Quotes are required if any character needs escaping, or there are
    // spaces in the string (spaces do not require escaping)
    bool quotes_required = false;
    for (const char& c : val) {
        if (char_requires_escaping(c) || c == ' ') {
            quotes_required = true;
            break;
        }
    }

    if (!quotes_required) {
        return val;
    }

    std::stringstream new_string;
    new_string << '\"';
    for (const char& c : val) {
        if (char_requires_escaping(c)) {
            new_string << '\\';
        }
        new_string << c;
    }
    new_string << '\"';

    return new_string.str();
}

template <typename T>
inline void output_property_names(ostream& out, const string& indentation,
                                  const map<string, T>& properties)
{
    for (const auto& p : properties) {
        out << indentation << p.first << "\n";
    }
}

template <typename T>
inline void output_property_values(ostream& out, const string& indentation,
                                   const map<string, T>& properties)
{
    for (const auto& p : properties) {
        out << indentation << local_to_string(p.second) << "\n";
    }
}

template <typename T>
void output_indexed_property_values(ostream& out,
                                    const map<string, T>& properties,
                                    unsigned int index)
{
    for (const auto& p : properties) {
        const auto& property = p.second;
        if (property->isDefined(index)) {
            out << ' ' << local_to_string(property->at(index));
        } else {
            out << " <>";
        }
    }
}

template <typename T>
bool maps_indexed_props_equal(const T& lmap, const T& rmap)
{
    if (rmap.size() != lmap.size())
        return false;
    auto diff = std::mismatch(
        lmap.begin(), lmap.end(), rmap.begin(),
        [](decltype(*begin(lmap)) l, decltype(*begin(lmap)) r) {
            return l.first == r.first && *(l.second) == *(r.second);
        });
    if (diff.first != lmap.end())
        return false;
    return true;
}
} // namespace

void Block::write(ostream& out, unsigned int current_indentation) const
{

    string root_indentation = string(current_indentation, ' ');
    current_indentation += 2;
    string indentation = string(current_indentation, ' ');

    const bool has_data = !m_bmap.empty() || !m_rmap.empty() ||
                          !m_imap.empty() || !m_smap.empty();

    out << root_indentation << getName() << " {\n";

    if (has_data) {
        output_property_names(out, indentation, m_bmap);
        output_property_names(out, indentation, m_rmap);
        output_property_names(out, indentation, m_imap);
        output_property_names(out, indentation, m_smap);

        out << indentation + ":::\n";

        output_property_values(out, indentation, m_bmap);
        output_property_values(out, indentation, m_rmap);
        output_property_values(out, indentation, m_imap);
        output_property_values(out, indentation, m_smap);
    }

    if (hasIndexedBlockData()) {
        const auto block_names = m_indexed_block_map->getBlockNames();
        for (const auto& name : block_names) {
            const auto& indexed_block =
                m_indexed_block_map->getIndexedBlock(name);
            indexed_block->write(out, current_indentation);
        }
    }

    for (const auto& p : m_sub_block) {
        const auto& sub_block = p.second;
        sub_block->write(out, current_indentation);
    }

    out << root_indentation << "}\n\n";

    return;
}

string Block::toString() const
{
    ostringstream stream;
    write(stream);

    return stream.str();
}

shared_ptr<const IndexedBlock> Block::getIndexedBlock(const string& name) const
{
    if (!hasIndexedBlockData()) {
        throw out_of_range("Indexed block not found: " + name);
    }
    return const_pointer_cast<const IndexedBlock>(
        m_indexed_block_map->getIndexedBlock(name));
}

bool real_map_equal(const map<string, double>& rmap1,
                    const map<string, double>& rmap2)
{
    if (rmap1.size() != rmap2.size())
        return false;
    for (const auto& p : rmap1) {
        if (rmap2.count(p.first) != 1)
            return false;
        if ((float) abs(p.second - rmap2.at(p.first)) > tolerance)
            return false;
    }

    return true;
}

bool Block::operator==(const Block& rhs) const
{
    if (m_bmap != rhs.m_bmap)
        return false;
    if (!real_map_equal(m_rmap, rhs.m_rmap))
        return false;
    if (m_imap != rhs.m_imap)
        return false;
    if (m_smap != rhs.m_smap)
        return false;
    if (m_sub_block != rhs.m_sub_block)
        return false;
    if (!(*m_indexed_block_map == *(rhs.m_indexed_block_map)))
        return false;
    return true;
}

bool IndexedBlockMapI::operator==(const IndexedBlockMapI& rhs) const
{
    const auto& block_names = getBlockNames();
    for (const auto& name : block_names) {
        if (!rhs.hasIndexedBlock(name))
            return false;
        const auto& block1 = rhs.getIndexedBlock(name);
        const auto& block2 = getIndexedBlock(name);
        if (*block1 != *block2)
            return false;
    }
    return true;
}

bool IndexedBlockMap::hasIndexedBlock(const string& name) const
{
    return m_indexed_block.find(name) != m_indexed_block.end();
}

shared_ptr<const IndexedBlock>
IndexedBlockMap::getIndexedBlock(const string& name) const
{
    auto block_iter = m_indexed_block.find(name);
    if (block_iter != m_indexed_block.end()) {
        return const_pointer_cast<const IndexedBlock>(block_iter->second);
    } else {
        throw out_of_range("Indexed block not found: " + name);
    }
}

bool BufferedIndexedBlockMap::hasIndexedBlock(const string& name) const
{
    if (m_indexed_buffer.find(name) != m_indexed_buffer.end()) {
        return true;
    } else if (m_indexed_block.find(name) != m_indexed_block.end()) {
        return true;
    } else {
        return false;
    }
}

shared_ptr<const IndexedBlock>
BufferedIndexedBlockMap::getIndexedBlock(const string& name) const
{
    auto itb = m_indexed_block.find(name);
    if (itb != m_indexed_block.end()) {
        return itb->second;
    }

    auto itbb = m_indexed_buffer.find(name);
    if (itbb == m_indexed_buffer.end()) {
        throw out_of_range("Indexed block not found: " + name);
    } else {
        shared_ptr<const IndexedBlock> ib(itbb->second->getIndexedBlock());
        return ib;
    }
}

template <>
EXPORT_MAEPARSER void
IndexedBlock::setProperty<BoolProperty>(const string& name,
                                        shared_ptr<IndexedBoolProperty> value)
{
    set_indexed_property<IndexedBoolProperty>(m_bmap, name, std::move(value));
}

template <>
EXPORT_MAEPARSER void
IndexedBlock::setProperty<double>(const string& name,
                                  shared_ptr<IndexedProperty<double>> value)
{
    set_indexed_property<IndexedProperty<double>>(m_rmap, name,
                                                  std::move(value));
}

template <>
EXPORT_MAEPARSER void
IndexedBlock::setProperty<int>(const string& name,
                               shared_ptr<IndexedProperty<int>> value)
{
    set_indexed_property<IndexedProperty<int>>(m_imap, name, std::move(value));
}

template <>
EXPORT_MAEPARSER void
IndexedBlock::setProperty<string>(const string& name,
                                  shared_ptr<IndexedProperty<string>> value)
{
    set_indexed_property<IndexedProperty<string>>(m_smap, name,
                                                  std::move(value));
}

size_t IndexedBlock::size() const
{
    size_t count = 0;
    // To save memory, not all maps will have max index count for the block,
    // so we must find the max size of all maps in the block.
    for (const auto& p : m_bmap)
        count = max(p.second->size(), count);
    for (const auto& p : m_imap)
        count = max(p.second->size(), count);
    for (const auto& p : m_rmap)
        count = max(p.second->size(), count);
    for (const auto& p : m_smap)
        count = max(p.second->size(), count);

    return count;
}

void IndexedBlock::write(ostream& out, unsigned int current_indentation) const
{
    string root_indentation = string(current_indentation, ' ');
    string indentation = string(current_indentation + 2, ' ');

    const bool has_data = !m_bmap.empty() || !m_rmap.empty() ||
                          !m_imap.empty() || !m_smap.empty();

    out << root_indentation << getName() << "[" << to_string((int) size())
        << "] {\n";

    if (has_data) {
        out << indentation + "# First column is Index #\n";

        output_property_names(out, indentation, m_bmap);
        output_property_names(out, indentation, m_rmap);
        output_property_names(out, indentation, m_imap);
        output_property_names(out, indentation, m_smap);

        out << indentation + ":::\n";

        for (unsigned int i = 0; i < size(); ++i) {
            out << indentation << i + 1;
            output_indexed_property_values(out, m_bmap, i);
            output_indexed_property_values(out, m_rmap, i);
            output_indexed_property_values(out, m_imap, i);
            output_indexed_property_values(out, m_smap, i);
            out << endl;
        }

        out << indentation + ":::\n";
    }

    out << root_indentation << "}\n";

    return;
}

string IndexedBlock::toString() const
{
    ostringstream stream;
    write(stream);

    return stream.str();
}

template <typename T>
bool IndexedProperty<T>::operator==(const IndexedProperty<T>& rhs) const
{
    if (m_is_null == nullptr || rhs.m_is_null == nullptr) {
        if ((m_is_null == nullptr) != (rhs.m_is_null == nullptr))
            return false;
    } else if (*m_is_null != *(rhs.m_is_null)) {
        return false;
    }
    if (m_data != rhs.m_data)
        return false;
    return true;
}

// For doubles we need to implement our own comparator for the vectors to
// take precision into account
template <>
bool IndexedProperty<double>::operator==(
    const IndexedProperty<double>& rhs) const
{
    if (m_is_null == nullptr || rhs.m_is_null == nullptr) {
        if ((m_is_null == nullptr) != (rhs.m_is_null == nullptr))
            return false;
    } else if (*m_is_null != *(rhs.m_is_null))
        return false;

    for (size_t i = 0; i < m_data.size(); ++i)
        if ((float) abs(m_data[i] - rhs.m_data[i]) > tolerance)
            return false;

    return true;
}

bool IndexedBlock::operator==(const IndexedBlock& rhs) const
{
    if (!maps_indexed_props_equal(m_bmap, rhs.m_bmap))
        return false;
    if (!maps_indexed_props_equal(m_imap, rhs.m_imap))
        return false;
    if (!maps_indexed_props_equal(m_rmap, rhs.m_rmap))
        return false;
    if (!maps_indexed_props_equal(m_smap, rhs.m_smap))
        return false;

    return true;
}

} // namespace mae
} // namespace schrodinger
