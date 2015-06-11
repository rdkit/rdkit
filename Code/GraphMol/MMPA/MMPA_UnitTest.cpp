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
//       products derived from this software without specific prior written permission.
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
#ifdef WIN32
#include <Windows.h>
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <string>
#include <iostream>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h" //MOL single molecule !
#include "../FileParsers/MolSupplier.h" //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "../Substruct/SubstructMatch.h"

#include "MMPA.h"

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

#ifdef WIN32
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL

struct timezone {
    int  tz_minuteswest; // minutes W of Greenwich
    int  tz_dsttime;     // type of dst correction
};

static inline int gettimeofday(struct timeval *tv, struct timezone *tz) {
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv) {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        //converting file time to unix epoch
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tmpres /= 10;  //convert into microseconds
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz) {
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

static inline unsigned long long nanoClock (void) { // actually returns microseconds
    struct timeval t;
    gettimeofday(&t, (struct timezone*)0);
    return t.tv_usec + t.tv_sec * 1000000ULL;
}


void printTime() {
    unsigned long long t1 = nanoClock();
    double sec = double(t1-t0) / 1000000.;
    printf("Time elapsed %.4lf seconds\n", sec);
    t0 = nanoClock();
}

std::string getSmilesOnly(const char* smiles, std::string* id=0) { // remove label, because RDKit parse FAILED
    const char* sp = strchr(smiles,' ');
    unsigned n = (sp ? sp-smiles+1 : strlen(smiles));
    if(id)
        *id = std::string(smiles+(n--));
    return std::string(smiles, n);
}

// UNIT Test Set:
//=========================================================================

void test1() {
    BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
    BOOST_LOG(rdInfoLog) << "MMPA test1()\n" << std::endl;
/*
// DEBUG PRINT
    RWMol *m = SmartsToMol("[*:1]C.[*:1]c1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]");
    Atom  *a = m->getAtomWithIdx(0);
    std::cout<<"DEBUG: "<< MolToSmiles(*m, true) <<"\n";
    delete m;
//-----
*/
    const char* smi[] = {
        "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-] ZINC21984717",
    };

    const char* fs[] = { // 16 reordered reference results
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
      "Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-],ZINC21984717,O=C(NC[*:2])[*:1],Cc1c([*:1])[n+](=O)c2ccccc2n1[O-].OC[*:2]"
    };

    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++) {
        static const std::string es("NULL");
        std::string id;
        std::string smiles = getSmilesOnly(smi[i], &id);
        ROMol* mol = SmilesToMol( smiles );
        std::vector< std::pair<ROMOL_SPTR,ROMOL_SPTR> > res;
        t0 = nanoClock();
        RDKit::MMPA::fragmentMol(*mol, res);
        printTime();
        delete mol;
        std::map<size_t, size_t> fs2res;
        std::cout << "\nTEST "<< i+1 << " mol: " << smi[i] <<"\n";
        for(size_t j=0; j<res.size(); j++) {
//            std::cout <<"   "<< j+1 << ": ";
//            std::cout << (res[j].first.get()  ? MolToSmiles(*res[j].first ) : es) <<", ";
//            std::cout << (res[j].second.get() ? MolToSmiles(*res[j].second) : es) <<"\n";
            std::stringstream ss;
            ss << smiles << "," << id << ",";
            ss << (res[j].first.get()  ? MolToSmiles(*res[j].first, true) : "") <<",";
            ss << (res[j].second.get() ? MolToSmiles(*res[j].second,true) : "");
            bool failed = true;
            size_t matchedRefRes = -1;
            for(size_t r=0; r < sizeof(fs)/sizeof(fs[0]); r++) {
                if(0==strcmp(std::string(ss.str()).c_str(), fs[r])) { // PASSED
                    failed = false;
                    matchedRefRes = r;
                    fs2res[r] = j;
                    break;
                }
            }
            if(j<9)
                std::cout << " ";
            if(failed) { //0!=strcmp(std::string(ss.str()).c_str(), fs[j])) { // FAILED
                std::cout << j+1 << ":*FAILED* " << ss.str() <<"\n";//<< "FS: " << fs[j] <<"\n";
//tmp            TEST_ASSERT(0==strcmp(std::string(ss.str()).c_str(), fs[j]));
            }
            else
                std::cout << j+1 << ": PASSED. matchedRefRes = "<< matchedRefRes+1 <<"\n";//ok: << "ss: " << ss.str() <<"\n";
            std::cout.flush();
        }
        std::cout << "\n --- UNMATCHED Reference RESULTS: --- \n";
        for(size_t r=0; r < sizeof(fs)/sizeof(fs[0]); r++) {
            if(fs2res.end() == fs2res.find(r))
                std::cout <<(r<9?" ":"")<< r+1 << ": " << fs[r] <<"\n";
        }
    }
    std::cout << " -----------------------------------\n";
    BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//====================================================================================================
//====================================================================================================

int main(int argc, const char* argv[]) {
    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    BOOST_LOG(rdInfoLog) << "MMPA Unit Test \n";

// use maximum CPU resoures to increase time measuring accuracy and stability in multi process environment
#ifdef WIN32
//    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
    SetThreadPriority(GetCurrentThread (), THREAD_PRIORITY_HIGHEST );
#else
    setpriority(PRIO_PROCESS, getpid(), -20);
#endif

    T0 = nanoClock();
    t0 = nanoClock();

    test1();
/*
    unsigned long long t1 = nanoClock();
    double sec = double(t1-T0) / 1000000.;
    printf("TOTAL Time elapsed %.4lf seconds\n", sec);
*/
    BOOST_LOG(rdInfoLog) << "*******************************************************\n";
    return 0;
}
