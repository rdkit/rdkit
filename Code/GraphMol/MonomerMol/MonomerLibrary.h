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
#include <tuple>

#include "boost/noncopyable.hpp"

#include <RDGeneral/export.h>

namespace RDKit
{

class RDKIT_MONOMERMOL_EXPORT MonomerLibrary : public boost::noncopyable
{
  public:
    // Returns: {helm_symbol, smiles, monomer_class}
    using helm_info_t = std::optional<std::tuple<std::string, std::string, std::string>>;

    MonomerLibrary(std::string_view database_path);
    MonomerLibrary();

    ~MonomerLibrary();

    [[nodiscard]] std::optional<std::string> getMonomerSmiles(
        const std::string& monomer_id,
        const std::string& monomer_class);

    [[nodiscard]] helm_info_t
    getHelmInfo(const std::string& three_letter_code);

    [[nodiscard]] std::optional<std::string>
    getPdbCode(const std::string& helm_symbol, const std::string& monomer_class);

  // private:
    // sqlite3* m_db;
};
} // namespace RDKit
