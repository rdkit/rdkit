//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once

#include <optional>
#include <string>
#include <string_view>

#include "boost/noncopyable.hpp"

#include "MonomerMol.h"

namespace RDKit
{

class [[nodiscard]] MonomerDatabase : public boost::noncopyable
{
  public:
    using helm_info_t = std::optional<std::tuple<std::string, std::string, ChainType>>;

    MonomerDatabase(std::string_view database_path);
    MonomerDatabase();

    ~MonomerDatabase();

    [[nodiscard]] std::optional<std::string> getMonomerSmiles(std::string monomer_id,
                                                      ChainType monomer_type);

    [[nodiscard]] helm_info_t
    getHelmInfo(const std::string& three_letter_code);

    [[nodiscard]] std::optional<std::string>
    getPdbCode(const std::string& helm_symbol, ChainType monomer_type);

  // private:
    // sqlite3* m_db;
};
} // namespace RDKit
