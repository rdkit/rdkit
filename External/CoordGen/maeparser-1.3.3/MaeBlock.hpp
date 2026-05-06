#pragma once

#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include "MaeParserConfig.hpp"

namespace schrodinger
{
namespace mae
{
using BoolProperty = uint8_t;

template <typename T>
inline const T& get_property(const std::map<std::string, T>& map,
                             const std::string& name)
{
    auto iter = map.find(name);
    if (iter == map.end()) {
        throw std::out_of_range("Key not found: " + name);
    } else {
        return iter->second;
    }
}

// Forward declaration.
class IndexedBlockBuffer;
class IndexedBlock;

class EXPORT_MAEPARSER IndexedBlockMapI
{
  public:
    virtual ~IndexedBlockMapI() = default;

    virtual bool hasIndexedBlock(const std::string& name) const = 0;

    virtual std::shared_ptr<const IndexedBlock>
    getIndexedBlock(const std::string& name) const = 0;

    virtual std::vector<std::string> getBlockNames() const = 0;
    bool operator==(const IndexedBlockMapI& rhs) const;
};

class EXPORT_MAEPARSER IndexedBlockMap : public IndexedBlockMapI
{
    std::map<std::string, std::shared_ptr<IndexedBlock>> m_indexed_block;

  public:
    bool hasIndexedBlock(const std::string& name) const override;

    std::shared_ptr<const IndexedBlock>
    getIndexedBlock(const std::string& name) const override;

    std::vector<std::string> getBlockNames() const override
    {
        std::vector<std::string> rval;
        for (const auto& p : m_indexed_block) {
            rval.push_back(p.first);
        }

        return rval;
    }

    /**
     * Add an IndexedBlock to the map.
     */
    void addIndexedBlock(const std::string& name,
                         std::shared_ptr<IndexedBlock> indexed_block)
    {
        m_indexed_block[name] = std::move(indexed_block);
    }
};

class EXPORT_MAEPARSER BufferedIndexedBlockMap : public IndexedBlockMapI
{
  private:
    std::map<std::string, std::shared_ptr<IndexedBlock>> m_indexed_block;
    std::map<std::string, std::shared_ptr<IndexedBlockBuffer>> m_indexed_buffer;

  public:
    bool hasIndexedBlock(const std::string& name) const override;

    std::shared_ptr<const IndexedBlock>
    getIndexedBlock(const std::string& name) const override;

    std::vector<std::string> getBlockNames() const override
    {
        std::vector<std::string> rval;
        for (const auto& p : m_indexed_buffer) {
            rval.push_back(p.first);
        }

        return rval;
    }

    /**
     * Add an IndexedBlockBuffer to the map, which can be used to retrieve an
     * IndexedBlock.
     */
    void addIndexedBlockBuffer(const std::string& name,
                               std::shared_ptr<IndexedBlockBuffer> block_buffer)
    {
        m_indexed_buffer[name] = std::move(block_buffer);
    }
};

class EXPORT_MAEPARSER Block
{
  private:
    const std::string m_name;

    std::map<std::string, BoolProperty> m_bmap;
    std::map<std::string, double> m_rmap;
    std::map<std::string, int> m_imap;
    std::map<std::string, std::string> m_smap;
    std::map<std::string, std::shared_ptr<Block>> m_sub_block;
    std::shared_ptr<IndexedBlockMapI> m_indexed_block_map;

  public:
    // Prevent copying.
    Block(const Block&) = delete;
    Block& operator=(const Block&) = delete;

    Block(std::string name)
        : m_name(std::move(name)), m_bmap(), m_rmap(), m_imap(), m_smap(),
          m_indexed_block_map(nullptr)
    {
    }

    const std::string& getName() const { return m_name; }

    std::string toString() const;

    void write(std::ostream& out, unsigned int current_indentation = 0) const;

    void setIndexedBlockMap(std::shared_ptr<IndexedBlockMapI> indexed_block_map)
    {
        m_indexed_block_map = std::move(indexed_block_map);
    }

    bool hasIndexedBlockData() const { return m_indexed_block_map != nullptr; }
    bool hasIndexedBlock(const std::string& name)
    {
        return hasIndexedBlockData() &&
               m_indexed_block_map->hasIndexedBlock(name);
    }

    std::shared_ptr<const IndexedBlock>
    getIndexedBlock(const std::string& name) const;

    void addBlock(std::shared_ptr<Block> b)
    {
        m_sub_block[b->getName()] = std::move(b);
    }

    /**
     * Check whether this block has a sub-block of the provided name.
     */
    bool hasBlock(const std::string& name)
    {
        std::map<std::string, std::shared_ptr<Block>>::const_iterator iter =
            m_sub_block.find(name);
        return (iter != m_sub_block.end());
    }

    /**
     * Retrieve a shared pointer to the named sub-block.
     */
    std::shared_ptr<Block> getBlock(const std::string& name) const
    {
        std::map<std::string, std::shared_ptr<Block>>::const_iterator iter =
            m_sub_block.find(name);
        if (iter == m_sub_block.end()) {
            throw std::out_of_range("Sub-block not found: " + name);
        } else {
            return iter->second;
        }
    }

    /**
     * Get the names of all non-indexed sub-blocks
     */
    std::vector<std::string> getBlockNames() const
    {
        std::vector<std::string> names;
        for (auto& n : m_sub_block) {
            names.push_back(n.first);
        }
        return names;
    }

    /**
     * Get the names of all indexed sub-blocks
     */
    std::vector<std::string> getIndexedBlockNames() const
    {
        return m_indexed_block_map->getBlockNames();
    }

    bool operator==(const Block& rhs) const;

    bool hasRealProperty(const std::string& name) const
    {
        return (m_rmap.find(name) != m_rmap.end());
    }

    double getRealProperty(const std::string& name) const
    {
        return get_property<double>(m_rmap, name);
    }

    void setRealProperty(const std::string& name, double value)
    {
        m_rmap[name] = value;
    }

    bool hasIntProperty(const std::string& name) const
    {
        return (m_imap.find(name) != m_imap.end());
    }

    int getIntProperty(const std::string& name) const
    {
        return get_property<int>(m_imap, name);
    }

    void setIntProperty(const std::string& name, int value)
    {
        m_imap[name] = value;
    }

    bool hasBoolProperty(const std::string& name) const
    {
        return (m_bmap.find(name) != m_bmap.end());
    }

    bool getBoolProperty(const std::string& name) const
    {
        return 1u == get_property<BoolProperty>(m_bmap, name);
    }

    void setBoolProperty(const std::string& name, bool value)
    {
        m_bmap[name] = static_cast<BoolProperty>(value);
    }

    bool hasStringProperty(const std::string& name) const
    {
        return (m_smap.find(name) != m_smap.end());
    }

    const std::string& getStringProperty(const std::string& name) const
    {
        return get_property<std::string>(m_smap, name);
    }

    void setStringProperty(const std::string& name, std::string value)
    {
        m_smap[name] = std::move(value);
    }

    template <typename T> const std::map<std::string, T>& getProperties() const;
};

template <typename T> class IndexedProperty
{
  private:
    std::vector<T> m_data;
    boost::dynamic_bitset<>* m_is_null;

  public:
    // Prevent copying.
    IndexedProperty(const IndexedProperty<T>&) = delete;
    IndexedProperty& operator=(const IndexedProperty<T>&) = delete;

    using size_type = typename std::vector<T>::size_type;

    /**
     * Construct an IndexedProperty from a reference to a vector of data.
     * This swaps out the data of the input vector.
     *
     * The optional boost::dynamic_bitset is owned by the created object.
     */
    explicit IndexedProperty(std::vector<T>& data,
                                boost::dynamic_bitset<>* is_null = nullptr)
        : m_data(), m_is_null(is_null)
    {
        m_data.swap(data);
    }

    ~IndexedProperty()
    {
        if (m_is_null != nullptr) {
            delete m_is_null;
        }
    }

    bool operator==(const IndexedProperty<T>& rhs) const;

    size_type size() const { return m_data.size(); }

    bool hasUndefinedValues() const
    {
        return (m_is_null != NULL && m_is_null->any());
    }

    bool isDefined(size_type index) const
    {
        if (m_is_null == nullptr) {
            // Use of assert matches out-of-bounds behavior for dynamic_bitset.
            assert(index < m_data.size());
            return true;
        } else {
            return !m_is_null->test(index);
        }
    }

    void undefine(size_type index)
    {
        if (m_is_null == NULL) {
            m_is_null = new boost::dynamic_bitset<>(m_data.size());
        }
        m_is_null->set(index);
    }

    inline T& operator[](size_type index)
    {
        if (m_is_null && m_is_null->test(index)) {
            throw std::runtime_error("Indexed property value undefined.");
        }
        return m_data[index];
    }

    inline const T& operator[](size_type index) const
    {
        if (m_is_null && m_is_null->test(index)) {
            throw std::runtime_error("Indexed property value undefined.");
        }
        return m_data[index];
    }

    inline T& at(size_type index) { return operator[](index); }

    inline const T& at(size_type index) const { return operator[](index); }

    inline const T& at(size_type index, const T& default_) const
    {
        if (m_is_null && m_is_null->test(index)) {
            return default_;
        }
        return m_data[index];
    }

    void set(size_type index, const T& value)
    {
        m_data[index] = value;
        if (m_is_null != NULL && m_is_null->test(index)) {
            m_is_null->reset(index);
        }
    }

    const std::vector<T>& data() const { return m_data; }
    const boost::dynamic_bitset<>* nullIndices() const { return m_is_null; }
};

using IndexedRealProperty = IndexedProperty<double>;
using IndexedIntProperty = IndexedProperty<int>;
using IndexedBoolProperty = IndexedProperty<BoolProperty>;
using IndexedStringProperty = IndexedProperty<std::string>;

template <typename T>
inline std::shared_ptr<T>
get_indexed_property(const std::map<std::string, std::shared_ptr<T>>& map,
                     const std::string& name)
{
    auto iter = map.find(name);
    if (iter == map.end()) {
        return std::shared_ptr<T>(nullptr);
    } else {
        return iter->second;
    }
}

template <typename T>
inline void set_indexed_property(std::map<std::string, std::shared_ptr<T>>& map,
                                 const std::string& name,
                                 std::shared_ptr<T> value)

{
    map[name] = std::move(value);
}

class EXPORT_MAEPARSER IndexedBlock
{
  private:
    const std::string m_name;

    std::map<std::string, std::shared_ptr<IndexedBoolProperty>> m_bmap;
    std::map<std::string, std::shared_ptr<IndexedIntProperty>> m_imap;
    std::map<std::string, std::shared_ptr<IndexedRealProperty>> m_rmap;
    std::map<std::string, std::shared_ptr<IndexedStringProperty>> m_smap;

  public:
    // Prevent copying.
    IndexedBlock(const IndexedBlock&) = delete;
    IndexedBlock& operator=(const IndexedBlock&) = delete;

    /**
     * Create an indexed block.
     */
    IndexedBlock(std::string name)
        : m_name(std::move(name)), m_bmap(), m_imap(), m_rmap(), m_smap()
    {
    }

    size_t size() const;

    const std::string& getName() const { return m_name; }

    std::string toString() const;

    void write(std::ostream& out, unsigned int current_indentation = 0) const;

    bool operator==(const IndexedBlock& rhs) const;

    bool operator!=(const IndexedBlock& rhs) const
    {
        return !(operator==(rhs));
    }

    template <typename T>
    void setProperty(const std::string& name,
                     std::shared_ptr<IndexedProperty<T>> value);

    bool hasBoolProperty(const std::string& name) const
    {
        return (m_bmap.find(name) != m_bmap.end());
    }

    std::shared_ptr<IndexedBoolProperty>
    getBoolProperty(const std::string& name) const
    {
        return get_indexed_property<IndexedBoolProperty>(m_bmap, name);
    }

    void setBoolProperty(const std::string& name,
                         std::shared_ptr<IndexedBoolProperty> value)
    {
        set_indexed_property<IndexedBoolProperty>(m_bmap, name,
                                                  std::move(value));
    }

    bool hasIntProperty(const std::string& name) const
    {
        return (m_imap.find(name) != m_imap.end());
    }

    std::shared_ptr<IndexedIntProperty>
    getIntProperty(const std::string& name) const
    {
        return get_indexed_property<IndexedIntProperty>(m_imap, name);
    }

    void setIntProperty(const std::string& name,
                        std::shared_ptr<IndexedIntProperty> value)
    {
        set_indexed_property<IndexedIntProperty>(m_imap, name,
                                                 std::move(value));
    }

    bool hasRealProperty(const std::string& name) const
    {
        return (m_rmap.find(name) != m_rmap.end());
    }

    std::shared_ptr<IndexedRealProperty>
    getRealProperty(const std::string& name) const
    {
        return get_indexed_property<IndexedRealProperty>(m_rmap, name);
    }

    void setRealProperty(const std::string& name,
                         std::shared_ptr<IndexedRealProperty> value)
    {
        set_indexed_property<IndexedRealProperty>(m_rmap, name,
                                                  std::move(value));
    }

    bool hasStringProperty(const std::string& name) const
    {
        return (m_smap.find(name) != m_smap.end());
    }

    std::shared_ptr<IndexedStringProperty>
    getStringProperty(const std::string& name) const
    {
        return get_indexed_property<IndexedStringProperty>(m_smap, name);
    }

    void setStringProperty(const std::string& name,
                           std::shared_ptr<IndexedStringProperty> value)
    {
        set_indexed_property<IndexedStringProperty>(m_smap, name,
                                                    std::move(value));
    }

    template <typename T>
    const std::map<std::string, std::shared_ptr<IndexedProperty<T>>>&
    getProperties() const;
};

// Template specializations

template <>
inline const std::map<std::string, BoolProperty>&
Block::getProperties<BoolProperty>() const
{
    return m_bmap;
}

template <>
inline const std::map<std::string, int>& Block::getProperties<int>() const
{
    return m_imap;
}

template <>
inline const std::map<std::string, double>& Block::getProperties<double>() const
{
    return m_rmap;
}

template <>
inline const std::map<std::string, std::string>&
Block::getProperties<std::string>() const
{
    return m_smap;
}

template <>
inline const std::map<std::string,
                      std::shared_ptr<IndexedProperty<BoolProperty>>>&
IndexedBlock::getProperties() const
{
    return m_bmap;
}

template <>
inline const std::map<std::string, std::shared_ptr<IndexedProperty<int>>>&
IndexedBlock::getProperties() const
{
    return m_imap;
}

template <>
inline const std::map<std::string, std::shared_ptr<IndexedProperty<double>>>&
IndexedBlock::getProperties() const
{
    return m_rmap;
}

template <>
inline const std::map<std::string,
                      std::shared_ptr<IndexedProperty<std::string>>>&
IndexedBlock::getProperties() const
{
    return m_smap;
}

} // namespace mae
} // namespace schrodinger
