// $Id: MMPA_UnitTest.cpp $
//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef _MSC_VER
#include <RDGeneral/test.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/resource.h>
#endif

#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>
#include <iostream>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/test_fixtures.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h"  //MOL single molecule !
#include "../FileParsers/MolSupplier.h"  //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "../Substruct/SubstructMatch.h"

#include "MMPA.h"

using namespace RDKit;

static unsigned n_failed = 0;
static unsigned long long T0;
static unsigned long long t0;

#ifdef _MSC_VER
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL

struct timezone {
  int tz_minuteswest;  // minutes W of Greenwich
  int tz_dsttime;      // type of dst correction
};

static inline int gettimeofday(struct timeval *tv, struct timezone *tz) {
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;

  if (nullptr != tv) {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    // converting file time to unix epoch
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
    tmpres /= 10;  // convert into microseconds
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (nullptr != tz) {
    if (!tzflag) {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
  return 0;
}
#endif

static inline unsigned long long nanoClock(
    void) {  // actually returns microseconds
  struct timeval t;
  gettimeofday(&t, (struct timezone *)nullptr);
  return t.tv_usec + t.tv_sec * 1000000ULL;
}

void printTime() {
  unsigned long long t1 = nanoClock();
  double sec = double(t1 - t0) / 1000000.;
  printf("Time elapsed %.6lf seconds\n", sec);
  t0 = nanoClock();
}

std::string getSmilesOnly(
    const char *smiles,
    std::string *id = nullptr) {  // remove label, because RDKit parse FAILED
  const char *sp = strchr(smiles, ' ');
  unsigned n = (sp ? sp - smiles + 1 : strlen(smiles));
  if (id) {
    *id = std::string(smiles + (n--));
  }
  return std::string(smiles, n);
}

void debugTest1(const char *init_mol) {
  std::unique_ptr<RWMol> m(SmilesToMol(init_mol));
  std::cout << "INIT MOL: " << init_mol << "\n";
  std::cout << "CONV MOL: " << MolToSmiles(*m, true) << "\n";
}
/*
 * Work-around functions for RDKit canonical
 */
std::string createCanonicalFromSmiles(const char *smiles) {
  std::unique_ptr<RWMol> m(SmilesToMol(smiles));
  std::string res = MolToSmiles(*m, true);
  //    replaceAllMap(res);
  return res;
}

std::string createCanonicalFromSmiles(std::string smiles) {
  return createCanonicalFromSmiles(smiles.c_str());
}

// UNIT Test Set:
//=========================================================================
void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "MMPA test1()\n" << std::endl;

  // DEBUG PRINT.      MolToSmiles() fails with new RDKit version

  //-----

  const char *smi[] = {
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-] ZINC21984717",
  };
  /*
  C[*:1].O=C(NCCO)c1c(n([O-])c2ccccc2[n+]1=O)[*:1]
  Cc1c([n+](=O)c2ccccc2n1[O-])[*:1].O=C(NCCO)[*:1]
   *
  */

  const char *fs_sm[] = {
      // 15 reference result's SMILES.
      "",
      "C[*:1].O=C(NCCO)c1c(n([O-])c2ccccc2[n+]1=O)[*:1]",
      // ORIGINAL:        "","[*:1]C.[*:1]c1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]",

      "[*:1]c1c([*:2])[n+](=O)c2ccccc2n1[O-]",
      "[*:1]C.[*:2]C(=O)NCCO",
      "[*:2]CNC([*:1])=O",
      "[*:2]CO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]C[*:1]",
      "[*:2]CO.[*:1]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]NC([*:1])=O",
      "[*:2]CCO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]NC(=O)c1c([*:1])n([O-])c2ccccc2[n+]1=O",
      "[*:1]C.[*:2]CCO",
      "[*:2]CNC(=O)c1c([*:1])n([O-])c2ccccc2[n+]1=O",
      "[*:1]C.[*:2]CO",
      "[*:2]CCNC(=O)c1c([*:1])n([O-])c2ccccc2[n+]1=O",
      "[*:1]C.[*:2]O",

      "",
      "Cc1c([n+](=O)c2ccccc2n1[O-])[*:1].O=C(NCCO)[*:1]",
      // ORIGINAL:        "","[*:1]C(=O)NCCO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",

      "",
      "[*:1]CCO.[*:1]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "",
      "[*:1]CO.[*:1]CNC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "",
      "[*:1]O.[*:1]CCNC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]CC[*:1]",
      "[*:2]O.[*:1]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]CCNC([*:1])=O",
      "[*:2]O.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "[*:2]C[*:1]",
      "[*:2]O.[*:1]CNC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
  };
  char fs[15][256];  // 15 reference results with updated RDKit's SMARTS Writer
  const char *fs1[] = {
      // 15+1dup reordered reference results
      // Fix for updated RDKit. new SMILES for the same molecule:
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,C[*:1].O=C(NCCO)c1c("
      "n([O-])c2ccccc2[n+]1=O)[*:1]",
      // was
      // "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,[*:1]C.[*:1]c1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]",

      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]c1c([*:2])[n+](="
      "O)c2ccccc2n1[O-],[*:1]C.[*:2]C(=O)NCCO",
      // exchanged atom mapping labels (1<->2):
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]CNC([*:1])=O,[*:"
      "2]CO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      // was
      // "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]CNC([*:2])=O,[*:1]CO.[*:2]c1c(C)n([O-])c2ccccc2[n+]1=O",
      // exchanged atom mapping labels (1<->2):
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]C[*:1],[*:2]CO.["
      "*:1]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      // was
      // "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]C[*:2],[*:1]CO.[*:2]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      // exchanged atom mapping labels (1<->2):
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]NC([*:1])=O,[*:"
      "2]CCO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      // was
      // "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]NC([*:2])=O,[*:1]CCO.[*:2]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]NC(=O)c1c([*:1])"
      "n([O-])c2ccccc2[n+]1=O,[*:1]C.[*:2]CCO",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]CNC(=O)c1c([*:1]"
      ")n([O-])c2ccccc2[n+]1=O,[*:1]C.[*:2]CO",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]CCNC(=O)c1c([*:"
      "1])n([O-])c2ccccc2[n+]1=O,[*:1]C.[*:2]O",

      // Fix for updated RDKit. new SMILES for the same molecule:
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,Cc1c([n+](=O)"
      "c2ccccc2n1[O-])[*:1].O=C(NCCO)[*:1]",
      // was
      // "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,[*:1]C(=O)NCCO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",

      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,[*:1]CCO.[*:1]NC(=O)"
      "c1c(C)n([O-])c2ccccc2[n+]1=O",
      // #10
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,[*:1]CO.[*:1]CNC(=O)"
      "c1c(C)n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,[*:1]O.[*:1]CCNC(=O)"
      "c1c(C)n([O-])c2ccccc2[n+]1=O",
      // exchanged atom mapping labels (1<->2):
      ///-- dup of #5:
      ///"Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]NC([*:1])=O,[*:2]CCO.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]CC[*:1],[*:2]O.["
      "*:1]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]CCNC([*:1])=O,[*"
      ":2]O.[*:1]c1c(C)n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:2]C[*:1],[*:2]O.[*"
      ":1]CNC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      /* ORIGINALS of exchanged atom mapping labels (1<->2):
      ///dup of #5:
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]NC([*:2])=O,[*:1]CCO.[*:2]c1c(C)n([O-])c2ccccc2[n+]1=O",
              "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]CC[*:2],[*:1]O.[*:2]NC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
              "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]CCNC([*:2])=O,[*:1]O.[*:2]c1c(C)n([O-])c2ccccc2[n+]1=O",
              "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,[*:1]C[*:2],[*:1]O.[*:2]CNC(=O)c1c(C)n([O-])c2ccccc2[n+]1=O",
      */
  };

  // FIX RDKit update. But MolToSmiles() fails with new RDKit version
  for (size_t r = 0; r < sizeof(fs) / sizeof(fs[0]); r++) {
    //        strcpy(fs[r], fs1[r]);
    // MolToSmiles() fails with new RDKit version in DEBUG Build
    char core[256] = "", side[256];
    if (fs_sm[2 * r][0]) {
      RWMol *m = SmilesToMol(fs_sm[2 * r]);
      strcpy(core, MolToSmiles(*m, true).c_str());
      delete m;
    }
    if (fs_sm[2 * r + 1][0]) {
      RWMol *m = SmilesToMol(fs_sm[2 * r + 1]);
      strcpy(side, MolToSmiles(*m, true).c_str());
      delete m;
    }
    sprintf(fs[r], "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,%s,%s",
            core, side);
  }
  strcpy(fs[0], fs1[0]);
  strcpy(fs[8], fs1[8]);
  //-----------------------------------

  for (unsigned int i = 0; i < sizeof(smi) / sizeof(smi[0]); i++) {
    static const std::string es("NULL");
    std::string id;
    std::string smiles = getSmilesOnly(smi[i], &id);
    ROMol *mol = SmilesToMol(smiles);
    std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> res;
    t0 = nanoClock();
    RDKit::MMPA::fragmentMol(*mol, res, 3);
    printTime();
    delete mol;
    std::cout << "TEST FINISHED "
                 "!!!\n*******************************************************"
                 "\n"
                 "  --- VERIFY RESULTS ---\n";

    // === VERIFY RESULTS ===
    std::map<size_t, size_t> fs2res;
    std::cout << "\nTEST " << i + 1 << " mol: " << smi[i] << "\n";
    bool test_failed = false;
    for (size_t j = 0; j < res.size(); j++) {
      //            std::cout <<"   "<< j+1 << ": ";
      //            std::cout << (res[j].first.get()  ?
      //            MolToSmiles(*res[j].first ) : es) <<", ";
      //            std::cout << (res[j].second.get() ?
      //            MolToSmiles(*res[j].second) : es) <<"\n";
      std::stringstream ss;
      ss << smiles << "," << id << ",";
      ss << (res[j].first.get() ? MolToSmiles(*res[j].first, true) : "") << ",";
      ss << (res[j].second.get() ? MolToSmiles(*res[j].second, true) : "");
      size_t matchedRefRes = -1;
      bool failed = true;
      for (size_t r = 0; r < sizeof(fs) / sizeof(fs[0]); r++) {
        if (0 == strcmp(std::string(ss.str()).c_str(), fs[r])) {  // PASSED
          failed = false;
          matchedRefRes = r;
          fs2res[r] = j;
          break;
        }
      }
      if (j < 9) {
        std::cout << " ";
      }
      if (failed) {
        test_failed = true;
        std::cout << j + 1 << ": NOREF. Reference data NOT LISTED in test case."
                  << ss.str() << "\n";  // << "FS: " << fs[j] <<"\n";
      } else {
        std::cout << j + 1 << ": PASSED. matchedRefRes = " << matchedRefRes + 1
                  << "\n";  // ok: << "ss: " << ss.str() <<"\n";
      }
      std::cout.flush();
    }
    std::cout << "\n --- UNMATCHED Reference RESULTS: --- \n";
    for (size_t r = 0; r < sizeof(fs) / sizeof(fs[0]); r++) {
      if (fs2res.end() == fs2res.find(r)) {
        test_failed = true;
        std::cout << (r < 9 ? " " : "") << r + 1 << ": " << fs[r] << "\n";
      }
    }
    std::cout << " -----------------------------------\n"
              << "DO TEST_ASSERT():\n";
    if (test_failed) {
      n_failed++;
    }
    TEST_ASSERT(!test_failed);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "MMPA test2()\n" << std::endl;

  const char *smi[] = {
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-] ZINC21984717",
  };

  const char *fs[] = {
      // 15 reordered reference results
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=[n+]1c([*:2])c([*:"
      "1])n([O-])c2ccccc21,C[*:1].O=C(NCCO)[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,Cc1c(C(=O)NCC[*:1])["
      "n+](=O)c2ccccc2n1[O-].O[*:1]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,OCC[*:1].Cc1c(C(=O)"
      "N[*:1])[n+](=O)c2ccccc2n1[O-]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,C[*:1].O=C(NCCO)c1c("
      "[*:1])n([O-])c2ccccc2[n+]1=O",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C([*:1])[*:2],Cc1c(C("
      "=O)N[*:1])[n+](=O)c2ccccc2n1[O-].OC[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(N[*:2])c1c([*:1])"
      "n([O-])c2ccccc2[n+]1=O,C[*:1].OCC[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(N[*:2])[*:1],"
      "Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].OCC[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C(C[*:2])[*:1],Cc1c("
      "C(=O)N[*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,Cc1c(C(=O)NC[*:1])["
      "n+](=O)c2ccccc2n1[O-].OC[*:1]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,O=C(NCCO)[*:1].Cc1c("
      "[*:1])[n+](=O)c2ccccc2n1[O-]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C([*:1])[*:2],Cc1c(C("
      "=O)NC[*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NCC[*:2])c1c([*:"
      "1])n([O-])c2ccccc2[n+]1=O,C[*:1].O[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NC[*:2])c1c([*:1]"
      ")n([O-])c2ccccc2[n+]1=O,C[*:1].OC[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NCC[*:2])[*:1],"
      "Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NC[*:2])[*:1],"
      "Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].OC[*:2]"};

  for (unsigned int i = 0; i < sizeof(smi) / sizeof(smi[0]); i++) {
    static const std::string es("NULL");
    std::string id;
    std::string smiles = getSmilesOnly(smi[i], &id);
    std::unique_ptr<ROMol> mol(SmilesToMol(smiles));

    std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> res;

    t0 = nanoClock();
    /*
     * General test fragment
     */
    RDKit::MMPA::fragmentMol(*mol, res);

    printTime();
    std::map<size_t, size_t> fs2res;
    std::map<std::string, size_t> ref_map;
    std::stringstream ref_str;
    std::string s_token;
    /*
     * Create reference map
     */
    for (size_t r = 0; r < sizeof(fs) / sizeof(fs[0]); r++) {
      std::stringstream ss_token(fs[r]);

      int token_num = 0;
      ref_str.str("");
      while (getline(ss_token, s_token, ',')) {
        if (token_num == 2) {
          ref_str << createCanonicalFromSmiles(s_token.c_str()) << ",";
        }
        if (token_num == 3) {
          ref_str << createCanonicalFromSmiles(s_token.c_str());
        }
        token_num++;
      }
      ref_map[ref_str.str()] = r;
    }
    std::cout << "\nTEST " << i + 1 << " mol: " << smi[i] << "\n";
    bool has_failed = false;
    for (size_t res_idx = 0; res_idx < res.size(); res_idx++) {
      std::cout << "   " << res_idx + 1 << ": ";
      /*
       * Somehow canonical smiles does not return the same result after just
       * saving.
       * Workaround is: save -> load -> save
       */
      std::string first_res =
          (res[res_idx].first.get() ? createCanonicalFromSmiles(MolToSmiles(
                                          *res[res_idx].first, true))
                                    : "");
      std::string second_res =
          (res[res_idx].second.get() ? createCanonicalFromSmiles(MolToSmiles(
                                           *res[res_idx].second, true))
                                     : "");

      std::cout << (res[res_idx].first.get() ? first_res : es) << ", ";
      std::cout << (res[res_idx].second.get() ? second_res : es) << "\n";

      std::stringstream res_str;
      res_str << first_res << "," << second_res;

      if (res_idx < 9) {
        std::cout << " ";
      }

      if (ref_map.find(res_str.str()) != ref_map.end()) {
        size_t matchedRefRes = ref_map[res_str.str()];
        fs2res[matchedRefRes] = res_idx;
        std::cout << res_idx + 1
                  << ": PASSED. matchedRefRes = " << matchedRefRes + 1
                  << "\n";  // ok: << "ss: " << ss.str() <<"\n";
      } else {
        std::cout << res_idx + 1
                  << ": NOREF. Reference data NOT LISTED in test case."
                  << res_str.str() << "\n";  //<< "FS: " << fs[j] <<"\n";
        has_failed = true;
      }
      std::cout.flush();
    }
    if (has_failed && fs2res.size() < sizeof(fs) / sizeof(fs[0])) {
      n_failed++;
      std::cout << "\n --- UNMATCHED Reference RESULTS: --- \n";
      for (size_t r = 0; r < sizeof(fs) / sizeof(fs[0]); r++) {
        if (fs2res.end() == fs2res.find(r)) {
          std::cout << (r < 9 ? " " : "") << r + 1 << ": " << fs[r] << "\n";
        }
      }
    } else {
      std::cout << "\n --- ALL PASSED --- \n";
    }
  }
  std::cout << " -----------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\tDone" << std::endl;
}

//====================================================================================================

void doTest(const char *smi, const char *fs[], unsigned fs_size) {
  static const std::string es("NULL");
  std::string id;
  std::string smiles = getSmilesOnly(smi, &id);
  std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> res;

  std::cout << "\nTEST mol: " << id << " " << smi << "\n";
  t0 = nanoClock();
  RDKit::MMPA::fragmentMol(*mol, res);
  printTime();

  // Create reference map
  std::map<size_t, size_t> fs2res;
  std::map<std::string, size_t> ref_map;
  std::stringstream ref_str;
  std::string s_token;
  for (size_t r = 0; r < fs_size; r++) {
    std::stringstream ss_token(fs[r]);

    int token_num = 0;
    ref_str.str("");
    while (getline(ss_token, s_token, ',')) {
      if (token_num == 2) {
        ref_str << createCanonicalFromSmiles(s_token.c_str()) << ",";
      }
      if (token_num == 3) {
        ref_str << createCanonicalFromSmiles(s_token.c_str());
      }
      token_num++;
    }
    ref_map[ref_str.str()] = r;
  }

  bool has_failed = false;
  for (size_t res_idx = 0; res_idx < res.size(); res_idx++) {
    if (res_idx < 9) {
      std::cout << " ";
    }
    std::cout << res_idx + 1 << ": res= ";
    /*
     * Somehow canonical smiles does not return the same result after just
     * saving.
     * Workaround is: save -> load -> save
     */
    std::string first_res =
        (res[res_idx].first.get()
             ? createCanonicalFromSmiles(MolToSmiles(*res[res_idx].first, true))
             : "");
    std::string second_res =
        (res[res_idx].second.get() ? createCanonicalFromSmiles(MolToSmiles(
                                         *res[res_idx].second, true))
                                   : "");

    std::cout << (res[res_idx].first.get() ? first_res : es) << ",";
    std::cout << (res[res_idx].second.get() ? second_res : es) << "\n";

    std::stringstream res_str;
    res_str << first_res << "," << second_res;

    if (res_idx < 9) {
      std::cout << " ";
    }

    if (ref_map.find(res_str.str()) != ref_map.end()) {
      size_t matchedRefRes = ref_map[res_str.str()];
      fs2res[matchedRefRes] = res_idx;
      std::cout << res_idx + 1
                << ": PASSED. matchedRefRes = " << matchedRefRes + 1
                << "\n";  // ok: << "ss: " << ss.str() <<"\n";
    } else {
      std::cout << res_idx + 1
                << ": NOREF. Reference data NOT LISTED in test case."
                << "\n";  // res_str.str() << "\n"; //<< "FS: " << fs[j] <<"\n";
      has_failed = true;
    }
    std::cout.flush();
  }
  if (has_failed && fs2res.size() < fs_size) {
    n_failed++;
    std::cout << "\n --- UNMATCHED Reference RESULTS: --- \n";
    for (size_t r = 0; r < fs_size; r++) {
      if (fs2res.end() == fs2res.find(r)) {
        std::cout << (r < 9 ? " " : "") << r + 1 << ": " << fs[r] << "\n";
      }
    }
  } else {
    std::cout << "\n --- ALL PASSED --- \n";
  }
  std::cout << " -----------------------------------\n";
}

void test3() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "MMPA test3()\n" << std::endl;

  {  // test2:
    const char *smi = "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-] ZINC21984717";
    // clang-format off
    const char* fs[] = {
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=[n+]1c([*:2])c([*:1])n([O-])c2ccccc21,C[*:1].O=C(NCCO)[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,Cc1c(C(=O)NCC[*:1])[n+](=O)c2ccccc2n1[O-].O[*:1]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,OCC[*:1].Cc1c(C(=O)N[*:1])[n+](=O)c2ccccc2n1[O-]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,C[*:1].O=C(NCCO)c1c([*:1])n([O-])c2ccccc2[n+]1=O",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C([*:1])[*:2],Cc1c(C(=O)N[*:1])[n+](=O)c2ccccc2n1[O-].OC[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(N[*:2])c1c([*:1])n([O-])c2ccccc2[n+]1=O,C[*:1].OCC[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(N[*:2])[*:1],Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].OCC[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C(C[*:2])[*:1],Cc1c(C(=O)N[*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,Cc1c(C(=O)NC[*:1])[n+](=O)c2ccccc2n1[O-].OC[*:1]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,,O=C(NCCO)[*:1].Cc1c([*:1])[n+](=O)c2ccccc2n1[O-]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,C([*:1])[*:2],Cc1c(C(=O)NC[*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NCC[*:2])c1c([*:1])n([O-])c2ccccc2[n+]1=O,C[*:1].O[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NC[*:2])c1c([*:1])n([O-])c2ccccc2[n+]1=O,C[*:1].OC[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NCC[*:2])[*:1],Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].O[*:2]",
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NC[*:2])[*:1],Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].OC[*:2]"};
    // clang-format on
    doTest(smi, fs, sizeof(fs) / sizeof(fs[0]));
  }

  {  // Case1 SIMPLE: (PASSED)
    const char *smi = "CC(N1CC1)C(=O) Case1-SIMPLE";
    // clang-format off
    const char* fs[] = {
        // from Greg's message
        "CC(N1CC1)C(=O),Case1-SIMPLE,C([*:1])([*:2])[*:3],C1CN1[*:2].C[*:1].O=C[*:3]",
        // from the results:
        "CC(N1CC1)C(=O),Case1-SIMPLE,,C[*:1].O=CC(N1CC1)[*:1]",
        "CC(N1CC1)C(=O),Case1-SIMPLE,O=CC([*:1])[*:2],C1CN1[*:2].C[*:1]",
        "CC(N1CC1)C(=O),Case1-SIMPLE,C1CN1C([*:1])[*:2],C[*:1].O=C[*:2]",
        "CC(N1CC1)C(=O),Case1-SIMPLE,,C1CN1[*:1].CC(C=O)[*:1]",
        "CC(N1CC1)C(=O),Case1-SIMPLE,CC([*:1])[*:2],C1CN1[*:1].O=C[*:2]",
        "CC(N1CC1)C(=O),Case1-SIMPLE,,CC(N1CC1)[*:1].O=C[*:1]",
    };
    // clang-format on
    doTest(smi, fs, sizeof(fs) / sizeof(fs[0]));
  }

  {  // Case1 (with additionally labeled central carbon [C:7]):
    const char *smi = "Cc1ccccc1NC(=O)[C:7](C)[NH+]1CCCC1 Case1";
    const char *fs[] = {
        "Cc1ccccc1NC(=O)[C:7](C)[NH+]1CCCC1,Case1,[*:1][C:7]([*:2])[*:3],C1CC["
        "NH+]([*:3])C1.C[*:2].Cc1ccccc1NC(=O)[*:1]"
        //            "Cc1ccccc1NC(=O)[C:7](C)[NH+]1CCCC1,Case1,[C:7]([*:1])([*:2])[*:3],C1CC[NH+]([*:1])C1.C[*:2].Cc1ccccc1NC(=O)[*:3]",
    };
    doTest(smi, fs, sizeof(fs) / sizeof(fs[0]));
  }

  {  // Case2:
    const char *smi = "O=C(OCc1ccccc1)C(O)c1ccccc1 Case2";
    const char *fs[] = {
        "O=C(OCc1ccccc1)C(O)c1ccccc1,Case2,C([*:1])([*:2])[*:3],O=C(OCc1ccccc1)"
        "[*:1].O[*:2].c1ccc([*:3])cc1",
    };
    doTest(smi, fs, sizeof(fs) / sizeof(fs[0]));
  }
  BOOST_LOG(rdInfoLog) << "\tDone" << std::endl;
}

void testCase_1() {
  // Case1 (with additionally labeled central carbon [C:7]):
  const char *smi = "Cc1ccccc1NC(=O)[C:7](C)[NH+]1CCCC1 Case1";
  const char *fs[] = {
      "Cc1ccccc1NC(=O)[C:7](C)[NH+]1CCCC1,Case1,[*:1][C:7]([*:2])[*:3],C1CC[NH+"
      "]([*:3])C1.C[*:2].Cc1ccccc1NC(=O)[*:1]"};
  doTest(smi, fs, sizeof(fs) / sizeof(fs[0]));
}
/*
void test4() {
    BOOST_LOG(rdInfoLog) << "-------------------------------------" <<
std::endl;
    BOOST_LOG(rdInfoLog) << "MMPA test4()\n" << std::endl;

    {
        const char* smi = "Cc1ccccc1NC(=O)[C:9](C)[NH+]1CCCC1 CASE-4-1";
        const char* fs[] = { ""
        };
        doTest(smi, fs, sizeof(fs)/sizeof(fs[0]));
    }

    {
        const char* smi = "c1ccccc1NC(=O)[C:9](C)[NH+]1CCCC1 CASE-4-2";
        const char* fs[] = { ""
        };
        doTest(smi, fs, sizeof(fs)/sizeof(fs[0]));
    }
    BOOST_LOG(rdInfoLog) << "\tDone" << std::endl;
}
*/
//====================================================================================================
//====================================================================================================

void testGithub6900() {
  UseLegacyStereoPerceptionFixture useLegacy(false);
  // auto mol = "CN1CCCN=C1/C=C/c1cccs1"_smiles;
  auto mol = "N/C=C/C"_smiles;
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> res;
  RDKit::MMPA::fragmentMol(*mol, res, 3);
}

void testMorganCodeHashChargeShift() {
  auto m = "Cc1ccccc1"_smiles;
  TEST_ASSERT(m);

  auto at = m->getAtomWithIdx(0);
  std::vector<unsigned long long> hashes;
  for (auto charge : {-2, -1, 0, 1, 2}) {
    at->setFormalCharge(charge);
    hashes.push_back(MMPA::detail::computeMorganCodeHash(*m));
  }

  for (unsigned i = 0; i < hashes.size() - 1; ++i) {
    for (unsigned j = i + 1; j < hashes.size(); ++j) {
      TEST_ASSERT(hashes[i] != hashes[j]);
    }
  }
}

int main() {
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  BOOST_LOG(rdInfoLog) << "MMPA Unit Test \n";

// use maximum CPU resoures to increase time measuring accuracy and stability in
// multi process environment
#ifdef _WIN32
  //    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
#else
  setpriority(PRIO_PROCESS, getpid(), -20);
#endif
  testGithub6900();
  T0 = nanoClock();
  t0 = nanoClock();

  testCase_1();
  // /*
  test2();
  test3();

  testMorganCodeHashChargeShift();

  //    test4();
  // */
  //    debugTest1("C[*:1].O=C(NCCO)c1c([*:1])n([O-])c2ccccc2[n+]1=O");
  //    debugTest1("C[*:1].O=C(NCCO)c1c(n([O-])c2ccccc2[n+]1=O)[*:1]");
  /*
      unsigned long long t1 = nanoClock();
      double sec = double(t1-T0) / 1000000.;
      printf("TOTAL Time elapsed %.4lf seconds\n", sec);
  */
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  if (0 != n_failed) {
    std::cout << n_failed << " TEST CASES FAILED \n";
    TEST_ASSERT(0 != n_failed);
  } else {
    std::cout << " --- ALL TEST CASES PASSED --- \n";
  }
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return 0;
}
