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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;
std::string pathName = std::string(getenv("RDBASE")) +
                       std::string("/Code/GraphMol/ossfuzz_testdata/");

void molFileParsingTests() {
  {  // Issue 24074
    std::string fname = pathName +
                        "clusterfuzz-testcase-minimized-mol_data_stream_to_mol_"
                        "fuzzer-5318835625000960";
    std::ifstream ins(fname, std::ifstream::in | std::ifstream::binary);
    unsigned int line;
    bool ok;
    try {
      MolDataStreamToMol(&ins, line);
      ok = false;
    } catch (const Invar::Invariant &) {
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
  molFileParsingTests();
  depicklingTests();
}
