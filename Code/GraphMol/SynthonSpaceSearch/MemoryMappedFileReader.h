//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef MEMORYMAPPEDFILEREADER_H
#define MEMORYMAPPEDFILEREADER_H

#include <string>

#include <RDGeneral/export.h>

namespace RDKit::SynthonSpaceSearch::details
{
  struct RDKIT_SYNTHONSPACESEARCH_EXPORT MemoryMappedFileReader
  {
    MemoryMappedFileReader() = delete;
    MemoryMappedFileReader(const std::string& filePath);
    MemoryMappedFileReader(const MemoryMappedFileReader&) = delete;
    MemoryMappedFileReader(MemoryMappedFileReader&& other);

    ~MemoryMappedFileReader();

    MemoryMappedFileReader& operator=(const MemoryMappedFileReader&) = delete;
    MemoryMappedFileReader& operator=(MemoryMappedFileReader&& other);

    char* d_mappedMemory{nullptr};
    size_t d_size{0};
  };
} // namespace RDKit::SynthonSpaceSearch::details

#endif  // MEMORYMAPPEDFILEREADER_H
