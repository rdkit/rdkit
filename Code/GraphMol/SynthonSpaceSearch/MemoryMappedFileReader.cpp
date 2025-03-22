//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/SynthonSpaceSearch/MemoryMappedFileReader.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>

#include <iostream>
#include <string>

#include <RDGeneral/RDLog.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif

namespace RDKit::SynthonSpaceSearch::details {
// This code is a lightly modified version of something provided by
// ChatGPT in response to the prompt:
// "in c++ can I use the same code for mmap on windows and linux?"
// ChatGPT assures me there are no license implications in using it.
// Accessed 26/2/2025.
MemoryMappedFileReader::MemoryMappedFileReader(const std::string &filePath) {
#ifdef _WIN32
  HANDLE hFile = CreateFile(filePath.c_str(), GENERIC_READ, 0, NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (hFile == INVALID_HANDLE_VALUE) {
    throw(std::runtime_error("Error opening file " + filePath + "."));
  }

  // Get file size
  LARGE_INTEGER fileSize;
  if (!GetFileSizeEx(hFile, &fileSize)) {
    CloseHandle(hFile);
    throw(std::runtime_error("Error reading file " + filePath + "."));
  }
  d_size = static_cast<size_t>(fileSize.QuadPart);

  // Create a file mapping
  HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
  if (hMapping == NULL) {
    CloseHandle(hFile);
    throw(std::runtime_error("Error reading file " + filePath + "."));
  }

  // Map the file into memory
  d_mappedMemory =
      static_cast<char *>(MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0));
  if (d_mappedMemory == NULL) {
    CloseHandle(hMapping);
    CloseHandle(hFile);
    throw(std::runtime_error("Error reading file " + filePath + "."));
  }

  CloseHandle(hMapping);  // Handle is no longer needed once mapped
  CloseHandle(hFile);     // File handle is no longer needed
#else
  int fd = open(filePath.c_str(), O_RDWR);
  if (fd == -1) {
    BOOST_LOG(rdErrorLog) << "Error opening file.\n";
    throw(std::runtime_error("Error opening file " + filePath + "."));
  }

  // Get file size
  struct stat fileStat;
  if (fstat(fd, &fileStat) == -1) {
    close(fd);
    throw(std::runtime_error("Error reading file " + filePath + "."));
  }
  d_size = static_cast<size_t>(fileStat.st_size);

  // Memory map the file
  d_mappedMemory =
      static_cast<char *>(mmap(NULL, d_size, PROT_READ, MAP_SHARED, fd, 0));
  if (d_mappedMemory == MAP_FAILED) {
    BOOST_LOG(rdErrorLog) << "Error mapping file.\n";
    close(fd);
    throw(std::runtime_error("Error reading file " + filePath + "."));
  }

  close(fd); // File descriptor is no longer needed
#endif
}

MemoryMappedFileReader::MemoryMappedFileReader(MemoryMappedFileReader &&other) {
  d_mappedMemory = other.d_mappedMemory;
  other.d_mappedMemory = nullptr;
  d_size = other.d_size;
  other.d_size = 0;
}

MemoryMappedFileReader::~MemoryMappedFileReader() {
#ifdef _WIN32
  // Windows-specific unmapping
  UnmapViewOfFile(d_mappedMemory);
#else
  // Linux-specific unmapping
  munmap(d_mappedMemory, d_size);
#endif
}

MemoryMappedFileReader &MemoryMappedFileReader::operator=(
  MemoryMappedFileReader &&other) {
  if (this != &other) {
#ifdef _WIN32
    // Windows-specific unmapping
    UnmapViewOfFile(d_mappedMemory);
#else
    // Linux-specific unmapping
    munmap(d_mappedMemory, d_size);
#endif
    d_mappedMemory = other.d_mappedMemory;
    d_size = other.d_size;
    other.d_mappedMemory = nullptr;
    other.d_size = 0;
  }
  return *this;
}
}; // namespace RDKit::SynthonSpaceSearch::details
