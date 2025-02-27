//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
// ChatGPT assures me there are no license implications in using it, although
// that brings Mandy Rice-Davies to mind.
// Accessed 26/2/2025.
void *createReadOnlyMemoryMapping(const std::string &filePath, size_t &size) {
#ifdef _WIN32
  HANDLE hFile = CreateFile(filePath.c_str(), GENERIC_READ, 0, NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (hFile == INVALID_HANDLE_VALUE) {
    BOOST_LOG(rdErrorLog) << "Error opening file.\n";
    return nullptr;
  }

  // Get file size
  LARGE_INTEGER fileSize;
  if (!GetFileSizeEx(hFile, &fileSize)) {
    BOOST_LOG(rdErrorLog) << "Error getting file size.\n";
    CloseHandle(hFile);
    return nullptr;
  }
  size = static_cast<size_t>(fileSize.QuadPart);

  // Create a file mapping
  HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
  if (hMapping == NULL) {
    BOOST_LOG(rdErrorLog) << "Error creating file mapping.\n";
    CloseHandle(hFile);
    return nullptr;
  }

  // Map the file into memory
  void *mappedMemory = MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0);
  if (mappedMemory == NULL) {
    BOOST_LOG(rdErrorLog) << "Error mapping view of file.\n";
    CloseHandle(hMapping);
    CloseHandle(hFile);
    return nullptr;
  }

  CloseHandle(hMapping);  // Handle is no longer needed once mapped
  CloseHandle(hFile);     // File handle is no longer needed

  return mappedMemory;

#else
  int fd = open(filePath.c_str(), O_RDWR);
  if (fd == -1) {
    BOOST_LOG(rdErrorLog) << "Error opening file.\n";
    return nullptr;
  }

  // Get file size
  struct stat fileStat;
  if (fstat(fd, &fileStat) == -1) {
    BOOST_LOG(rdErrorLog) << "Error getting file size.\n";
    close(fd);
    return nullptr;
  }
  size = static_cast<size_t>(fileStat.st_size);

  // Memory map the file
  void *mappedMemory = mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
  if (mappedMemory == MAP_FAILED) {
    BOOST_LOG(rdErrorLog) << "Error mapping file.\n";
    close(fd);
    return nullptr;
  }

  close(fd);  // File descriptor is no longer needed
  return mappedMemory;
#endif
}

void unmapMemory(void *mappedMemory, size_t size) {
#ifdef _WIN32
  // Windows-specific unmapping
  UnmapViewOfFile(mappedMemory);
#else
  // Linux-specific unmapping
  munmap(mappedMemory, size);
#endif
}

}