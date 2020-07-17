//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <fstream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <chrono>

using namespace RDKit;
std::string pathName = std::string(getenv("RDBASE")) +
                       std::string("/Code/GraphMol/ossfuzz_testdata/");

void molFileParsingTests() {
  {  // Issue 23694
    std::string fname = pathName +
                        "clusterfuzz-testcase-minimized-mol_data_stream_to_mol_"
                        "fuzzer-5986616752734208";
    std::ifstream ins(fname, std::ifstream::in | std::ifstream::binary);
    unsigned int line = 0;
    bool ok;
    try {
      MolDataStreamToMol(&ins, line);
      ok = false;
    } catch (const FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {  // Issue 24074
    std::string fname = pathName +
                        "clusterfuzz-testcase-minimized-mol_data_stream_to_mol_"
                        "fuzzer-5318835625000960";
    std::ifstream ins(fname, std::ifstream::in | std::ifstream::binary);
    unsigned int line = 0;
    bool ok;
    try {
      MolDataStreamToMol(&ins, line);
      ok = false;
    } catch (const FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
}

void depicklingTests() {
  {  // Issue 23896
    std::string fname = pathName +
                        "clusterfuzz-testcase-minimized-mol_deserialization_"
                        "fuzzer-6283983225880576";
    std::ifstream ins(fname, std::ifstream::in | std::ifstream::binary);
    unsigned int line;
    ROMol res;
    bool ok;
    try {
      MolPickler::molFromPickle(ins, res);
      ok = false;
    } catch (const std::runtime_error &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
}

int main() {
  RDLog::InitLogs();

#if 1
  molFileParsingTests();
  depicklingTests();
#else

  unsigned int nReps = 10000000;
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      int val = strtol("1234", nullptr, 10);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "strtol: " << elapsed_seconds.count() << "s\n";
  }
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      const char *txt = "1234";
      const char *ptr = txt;
      while (*ptr != '\x00') {
        if ((*ptr >= '0' && *ptr <= '9') || *ptr == '+' || *ptr == '-') {
          ++ptr;
          continue;
        } else {
          break;
        }
      }
      int val = strtol(txt, nullptr, 10);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "strtol2: " << elapsed_seconds.count() << "s\n";
  }
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      int val = boost::lexical_cast<int>("1234");
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "lexical: " << elapsed_seconds.count() << "s\n";
  }
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      const char *txt = "1234.4";
      double val = atof(txt);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "atof: " << elapsed_seconds.count() << "s\n";
  }
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      const char *txt = "1234.4";
      const char *ptr = txt;
      for (unsigned int j = 0; j < strlen(txt); ++j, ++ptr) {
        if ((*ptr >= '0' && *ptr <= '9') || *ptr == '+' || *ptr == '-' ||
            *ptr == '.') {
          continue;
        } else {
          break;
        }
      }
      double val = atof(txt);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "atof check: " << elapsed_seconds.count() << "s\n";
  }
  {
    auto start = std::chrono::steady_clock::now();
    for (auto i = 0u; i < nReps; ++i) {
      double val = boost::lexical_cast<double>("1234.4");
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "lexical: " << elapsed_seconds.count() << "s\n";
  }
#endif
}
