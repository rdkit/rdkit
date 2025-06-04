#pragma once

#include <optional>
#include <string>
#include <string_view>

#include "boost/noncopyable.hpp"

namespace RDKit
{

enum class ChainType;
class [[nodiscard]] MonomerDatabase : public boost::noncopyable
{
  public:
    using monomer_smiles_t = std::optional<std::string>;
    using helm_info_t = std::optional<std::pair<std::string, ChainType>>;

    MonomerDatabase(std::string_view database_path);

    ~MonomerDatabase();

    [[nodiscard]] monomer_smiles_t getMonomerSmiles(std::string monomer_id,
                                                      ChainType monomer_type);

    [[nodiscard]] helm_info_t
    getHelmInfo(const std::string& three_letter_code);

  // private:
    // sqlite3* m_db;
};
} // namespace RDKit
