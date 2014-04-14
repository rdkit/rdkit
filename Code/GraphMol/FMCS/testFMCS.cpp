// $Id: testFMCS.cpp $
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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
//#include <RDGeneral/RDLog.h>
//#include <RDGeneral/utils.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h" //MOL single molecule !
#include "../FileParsers/MolSupplier.h" //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "FMCS.h"

#include "DebugTrace.h" //#ifdef VERBOSE_STATISTICS_ON

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

#ifdef WIN32
//#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
//  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
//#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
//#endif
 
struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
 
  return 0;
}
#endif

unsigned long long nanoClock (void)   // actually returns microseconds
{
   struct timeval t;
   //struct timezone tz;
   gettimeofday(&t, (struct timezone*)0);
   return t.tv_usec + t.tv_sec * 1000000ULL;
}

void printTime()
{
    unsigned long long t1 = nanoClock();
    double sec = double(t1-t0) / 1000000.;
    printf("\nTime elapsed %.2f seconds\n", sec);
    t0 = nanoClock();
}

std::string getSmilesOnly(const char* smiles, std::string* id=0) // remove label, because RDKit parse FAILED
{
    const char* sp = strchr(smiles,' ');
    unsigned n = (sp ? sp-smiles+1 : strlen(smiles));
    if(id)
        *id = std::string(smiles+n);
    return std::string(smiles, n);
}
//====================================================================================================

void testFileMCSB(const char* test, unsigned timeout=30, std::vector<unsigned> test_N=std::vector<unsigned>())  // optional list of some tests for investigation
{
    std::vector<ROMOL_SPTR> mols; // IT CAN OCCUPY A LOT OF MEMORY. store SMILES only to reduce memory usage.
    char str [4096];
    std::string molFile, id;
    std::map<std::string, size_t>         molIdMap;
    std::vector<std::string>            smilesList;
    std::list< std::vector<std::string> > testCase;
    std::string referenceOutFile(test); referenceOutFile += ".REF.out";
    std::string outFile(test);
    if(!test_N.empty())
    {
        if(1==test_N.size())
            sprintf(str,".%u.out", test_N[0]);
        else
            sprintf(str,".%u-%u.out", test_N[0], test_N.back());
        outFile += str;
    }
    else
    {
        RDKit::FMCS::ConsoleOutputEnabled = false;
        outFile += ".Cpp.out";
    }

    std::string outSmilesFile(test); outSmilesFile += ".smiles";
    unsigned n=0, passed=0, failed=0, failed_less=0, timedout=0;
    double secTotal = 0.;

    std::vector<MCSResult>  referenceResults;
    std::vector<float>      referenceResultsTime;

    FILE* f = fopen(referenceOutFile.c_str(), "rt");
    if(!f)
        perror("Could not open reference test result file");
    else
    {
        std::cout<<"Loading reference test results ... \n";
        while(fgets(str, sizeof(str), f))
         if('#' != str[0])
         {
             char c;
             int  frag;
             float t;
             char mcs [1024];
             MCSResult res;
             sscanf(str, "%u %c %d %d %d %f %s", &n, &c, &frag, &res.NumAtoms, &res.NumBonds, &t, mcs);
             res.Canceled = ('.' != c);
             res.SmartsString = mcs;
             referenceResults.push_back(res);
             referenceResultsTime.push_back(t);
         }
    }
    fclose(f);

    f = fopen(test, "rt");
    if(!f)
    {
        perror("Could not open test case list MCSB file");
        exit(1);
    }
    {
        std::cout<<"Loading MCSB test list ... \n";
        if(fgets(str, sizeof(str), f))
        if(fgets(str, sizeof(str), f))
        {
            char* c = strrchr(str, '\n');   // remove LineFeed
            if(c)
                *c = '\0';
            c = strrchr(str, '\r');
            if(c)
                *c = '\0';
            molFile = str + 6;  // #File filename
        }
        std::cout<<"Molecules file:"<<molFile<<"\n";
        n = 0;
        while(fgets(str, sizeof(str), f))
         if('#' != str[0])
         { //str= "1 CHEMBL526291 CHEMBL498211 ..."
            char name [256];
            unsigned nn, len;
            n++;
            testCase.push_back(std::vector<std::string>());
            sscanf(str, "%u%n", &nn, &len);
            while('\0'!=*(str+len) && 1==sscanf(str+len, "%s%n", name, &nn))
            {
                len += nn;
                testCase.back().push_back(std::string(name));
            }
         }
        std::cout<<n<<" Test cases loaded\n";
    }
    fclose(f);

    f = fopen(molFile.c_str(), "rt");
    if(!f)
    {
        perror("Could not open molecules file");
        exit(1);
    }
    std::cout<<"Loading SMILES ... \n";
    n = 0;
    while(fgets(str, sizeof(str), f))
    {
        std::cout<<"\rLine: "<< ++n <<" ";
        if('#' != str[0] && ' ' != str[0] && '/' != str[0]) // commented to skip
        {
            char* c = strrchr(str, '\n');   // remove LineFeed
            if(c)
                *c = '\0';
            c = strrchr(str, '\r');
            if(c)
                *c = '\0';
           std::string sm = getSmilesOnly(str, &id);
           smilesList.push_back(sm);        // without Id and LineFeed
           mols.push_back(ROMOL_SPTR(SmilesToMol(sm))); // SmartsToMol ???
           molIdMap[id] = mols.size()-1;// index in mols
        }
    }
    fclose(f);

    printTime();

    f  = fopen(outFile.c_str(), "wt");
    if(!f)
    {
        perror("Could not create output file");
        exit(1);
    }
    FILE* fs  = fopen(outSmilesFile.c_str(), "wt");
    FILE* ft  = fopen((outFile+".stat.csv").c_str(), "wt");
    setvbuf(f , 0, _IOFBF, 4*1024);
    setvbuf(fs, 0, _IOFBF, 4*1024);
    setvbuf(ft, 0, _IOFBF,  4*1024); // small file
    if(ft)
        fprintf(ft, "N; Status; dAtoms; dBonds; t(sec); ref.t; Seed; MatchCall; AtomCmp; BondCmp\n"); //CSV Header

    n = 0;
    MCSParameters p;
    p.Timeout = timeout;
    p.Threshold = 1.0;

    fprintf(f, "#software RDKit C++ FMCS \n#options  timeout=%u threshold=%g\n", p.Timeout, p.Threshold);
    std::cout<<"Perform test cases ... \n";
    for(std::list< std::vector<std::string> >::const_iterator tc = testCase.begin(); tc != testCase.end(); tc++, n++)
    {
//        if(test_N != 0 && test_N != n+1)
        if(!test_N.empty() && test_N.end() == std::find(test_N.begin(), test_N.end(), n+1))
            continue;

        std::cout<<"\rTest: "<< n+1 <<" ";
        if(!test_N.empty()) // test case is listed
            std::cout<<"\n";

        std::vector<ROMOL_SPTR> tcmols;
        fprintf(f, "# %u Using ", n+1);
        if(fs)
            fprintf(fs, "\n//TEST %u\n", n+1);
        for(std::vector<std::string>::const_iterator mid = tc->begin(); mid != tc->end(); mid++)
        {
            std::map<std::string, size_t>::const_iterator id = molIdMap.find(*mid);
            if(molIdMap.end() == id)
                continue;
            size_t i = id->second;
            tcmols.push_back(mols[i]);
            fprintf(f, "%s ", mid->c_str());
            if(fs)
                fprintf(fs, "\"%s%s\",\n", smilesList[i].c_str(), mid->c_str());
        }
        fprintf(f, "\n");
//        ExecStatistics curStat = stat;          //to compute the difference for this test only
        unsigned long long tc0 = nanoClock();
        MCSResult res = findMCS(tcmols, &p);    // *** T E S T ***
        unsigned long long tc1 = nanoClock();
        double sec = double(tc1-tc0) / 1000000.; // without time of SMILES to ROMol conversion
        secTotal += sec;
        if(!test_N.empty())
            std::cout<<"\n" << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
        else if(!referenceResults[n].Canceled && !res.Canceled && (/*referenceResults[n].NumAtoms > res.NumAtoms ||*/ referenceResults[n].NumBonds > res.NumBonds))
            std::cout<<" - failed. LESS: "<<res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds ("<<referenceResults[n].NumAtoms - res.NumAtoms<<", "<<referenceResults[n].NumBonds-res.NumBonds<<")\n";

        if(!referenceResults.empty())
        {
            fprintf(f, "# %u REFCMP: time %s %.2f sec.\n", n+1, fabs(referenceResultsTime[n] - sec)<sec/25.?"EQUAL"
                :(referenceResultsTime[n] < sec ? (sec<.7?"slow":"SLOW")
                :(sec<.7 && referenceResultsTime[n]<0.7 ? "fast":"FAST")), referenceResultsTime[n]);
            if(!referenceResults[n].Canceled)// && !res.Canceled)
            {
                if(referenceResults[n].NumAtoms == res.NumAtoms && referenceResults[n].NumBonds == res.NumBonds)
                    fprintf(f, "# %u REFCMP: res  %s %s %u %u %s.\n", n+1, "PASSED"
                    , "-------"
                    , referenceResults[n].NumAtoms, referenceResults[n].NumBonds, referenceResults[n].SmartsString.c_str());
                else
                    fprintf(f, "# %u REFCMP: res  %s %s %u %u %s.\n", n+1, "FAILED"
                    , /*referenceResults[n].NumAtoms > res.NumAtoms ||*/ referenceResults[n].NumBonds > res.NumBonds ? "MISSING":"GREATER"
                    , referenceResults[n].NumAtoms, referenceResults[n].NumBonds, referenceResults[n].SmartsString.c_str());

                if(referenceResults[n].Canceled 
                 ||(referenceResults[n].NumAtoms == res.NumAtoms && referenceResults[n].NumBonds == res.NumBonds))
                    passed++;
                else if(res.Canceled)
                    timedout++;
                else
                {
                    if(referenceResults[n].NumBonds > res.NumBonds)
                        failed_less++;
                    failed++;
                }
            }
            else
                fprintf(f, "# %u REFCMP: res  ABSENT - timeout\n", n+1);
        }
        // 1 . 1 25 28 1.69 F-c1:c:c:
        fprintf(f, "%u %c %d %u %u %.2f %s\n", n+1, (res.Canceled ? 'F':'.'), 1 //number of fragments in the MCS
                                          , res.NumAtoms, res.NumBonds, sec, res.SmartsString.c_str());
        if(fs)
            fprintf(fs, "//# %u %c  %u %u %.2f sec MCS: %s\n", n+1, (res.Canceled ? 'F':'.'), res.NumAtoms, res.NumBonds, sec, res.SmartsString.c_str());
#ifdef xxVERBOSE_STATISTICS_ON
        if(ft) // statistic details
            fprintf(ft, "%u; %s; %d; %d; %.2f; %.2f; %u; %u; %u; %u\n",n+1
                      , !res.Canceled ? "ok" : referenceResults[n].Canceled ? "bad" : "TIMEOUT"
                      , referenceResults[n].Canceled ? 0 : res.NumAtoms - referenceResults[n].NumAtoms
                      , referenceResults[n].Canceled ? 0 : res.NumBonds - referenceResults[n].NumBonds
                      , sec, referenceResultsTime[n]
                      , stat.Seed - curStat.Seed
                      , stat.MatchCall - curStat.MatchCall
                      , stat.AtomCompareCalls - curStat.AtomCompareCalls
                      , stat.BondCompareCalls - curStat.BondCompareCalls
                   );
        stat.AtomCompareCalls = 0;  // 32 bit counter with very big value -> possible overflow
        stat.BondCompareCalls = 0;
#endif
    }
    fprintf(f, "#\n#\n# %u passed, %u failed, %u failed_less, %u timed out.\n# Total %.2f seconds, Average %.2f seconds, Average exclude timeouts about %.2f seconds.\n"
             , passed, failed, failed_less, timedout, secTotal, secTotal/n, (secTotal-30.6*timedout)/n);
#ifdef xxVERBOSE_STATISTICS_ON
    fprintf(f, "#\n# --- STATISTICS:---\n#             Total value   |   Average\n"
        "# Seeds Num %15u | %8u (amount of generated seeds)\n"
        "# BestSizeR %15u | %8u = %d%% (rejected by RemainingSize against BestSize seed)\n"
        "# MatchCall %15u | %8u (SubstructMatch function calls)\n"
        "# MatchTRUE %15u | %8u = %d%%\n"
#ifdef FAST_SUBSTRUCT_CACHE
//        "#HashCache  %15u | %8u keys\n"
//        "#HashCache  %15u | %8u entries\n"
        "# HashCacheFind  %15u | %8u \n"
        "# HashKeysFound  %15u | %8u = %d%% hash keys found \n"
        "# ExactMatchCall %15u | %8u (SubstructMatch function calls)\n"
        "# ExactMatchTRUE %15u | %8u \n"
#endif
        , stat.Seed, stat.Seed/n
        , stat.RemainingSizeRejected, stat.RemainingSizeRejected/n
        , 0==stat.Seed ? 0 : int((double)stat.RemainingSizeRejected / (double)stat.Seed *100.)
        , stat.MatchCall        , stat.MatchCall/n
        , stat.MatchCallTrue        , stat.MatchCallTrue/n
        , int((double)stat.MatchCallTrue / (double)stat.MatchCall *100.)
#ifdef FAST_SUBSTRUCT_CACHE
//        , stat.HashCacheKeysSize, stat.HashCacheKeysSize/n
//        , stat.HashCacheEntries , stat.HashCacheEntries/n
        , stat.FindHashInCache,    stat.FindHashInCache/n
        , stat.HashKeyFoundInCache,stat.HashKeyFoundInCache/n
        , 0==stat.FindHashInCache ? 0 : int((double)stat.HashKeyFoundInCache / (double)stat.FindHashInCache *100.)

        , stat.ExactMatchCall    , stat.ExactMatchCall/n
        , stat.ExactMatchCallTrue, stat.ExactMatchCallTrue/n
#endif
        );
#endif
    fclose(f);
    if(fs)
        fclose(fs);
    if(ft)
        fclose(ft);
    printTime();
}

//=========================================================================


//=========================================================================

void test1Basics()
{
//    BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
//    BOOST_LOG(rdInfoLog) << "FMCS test1Basics()" << std::endl;

    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
//         "CC1CCC(N)CC1", "CC1CC(C)CC(C)C1"  // OK test.sdf
//         "OC1CCC1", "OC1CCC1" // OK ==
//         "CC1CC(C=O=S)CC(C=N)C1", "CC1CC(C=O=S)CC(C=N)C1"  // OK test.sdf ++
//         "CC1CC(C=O=S)CC(N=C)C1", "CC1CC(C)CC(N)C1"  // OK test.sdf ++
//        "CC1CC(C=O=S3)C2C(N=C-O23)C1", "C7C1CC(C8)C7C(NO8)C1"  // OK test.sdf ++
//        "O1CCN(OSC)CO1-CC1CC(N=C-O23)C2C(C=O=S3)C1", "C7C1CC(C8)C7C(NO8)C1"  // OK MCS:{0 9 8 5 7 4 2 3 11}: Smarts=CCCC(N)CC(C)CO

//      "[O:8][C:1]1[C:2][C:3]([N:9])[C:4][C:5][C:6]1", "NC1(S-O-C)C[C:9](O-S-C)[C:7]CC1" // OK MCS:{0 1 7 2 6 3 4 5}: Smarts=NC1CCCC(O)C1 == N[CH:3]1[CH2:4][CH2:5][CH2:6][CH:1](O)[CH2:2]1
        //MCS: N[CH:3]1[CH2:4][CH2:5][CH2:6][CH:1][(OH:8)][CH2:2]1


// PASSED     ///FAILED (missing one bond to close ring):
// Python MCS = 26 bonds : COCc1cncc(c1):n:c1cccc(Oc2ccc(Cl)cc2)c1 WITH AROMATIZATION
// MCS 26: COCc1c-ncc(c1)nc1cccc(c1)Oc1ccc(Cl)cc1 24 atoms, 26 bonds Time elapsed 0.12 seconds
// MCS 25: cc(cc(cn)nc1cccc(Oc2ccc(Cl)cc2)c1)COC 24 atoms, 25 bonds. WITH AROMATIZATION !!!!
// MCS 16: COCC=CN=CC1=CC(=CC=C)C(=C)N1, 16 atoms, 16 bonds.      WITHOUT AROMATIZATION !!!!
///            "COCC1=C(N=CC2=C1C1=C(OC3=CC=C(Cl)C=C3)C=CC=C1N2)C(=O)OC(C)C",
///            "COCC1=CN=C(C(=O)OC(C)C)C2=C1C1=CC=C(OC3=CC=C(Cl)C=C3)C=C1N2",
// The SAME, but pre AROMATIZATED (else PRECONDITION Exception with Implicit Hs / 16 bonds only)
//            "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
//            "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12",
//
/* /TEST 4
"CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934",
"CN(C)c1ccc(CC(=O)NCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL152361",
"CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336",
"CN(C)c1ccc(CC(=O)NCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157429",
"CN(C)c1ccc(CC(=O)NCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL357551",
"CN(C)c1ccc(CC(=O)NCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL421974",
"CN(C)c1ccc(CC(NCCCCCC(NO)=O)=O)cc1 CHEMBL484488",
"CC(C)Cc1ccc(C(C)C(=O)NC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL564780",
"c1cc([N+]([O-])=O)ccc1CC(=O)NC1CCCCCC1 CHEMBL1553142",
"CC1(C)NC(C)(C)CC(NC(=O)Cc2ccccc2)C1 CHEMBL1703640",
//# 3 . 1 14 14 0.08 sec MCS: CCCCNC(=O)Cc1ccccc1
*/
/*
//TEST 5 FAILED
//.5 . 1 20 22 0.23 N(-c1:c:c:c:c:c:1)-C-c1:n(-C-c2:c:c:c:c:c:2):c:n:c:1
//MCS: ccccc-ccc(c)NCc1cncn1Cc1ccccc1 23 atoms, 24 bonds Time elapsed 6.39 seconds
"c1ccc(Cn2c(CNc3cc(-c4ccccc4)ccc3)cnc2)cc1 CHEMBL485450",   // QUERY
"Cc1ccc(Cn2c(CNc3cc(-c4ccccc4)ccc3)cnc2)cc1 CHEMBL498061",
"c1ccc(-c2ccc(Cn3c(CNc4cc(-c5ccccc5)ccc4)cnc3)cc2)cc1 CHEMBL485449",
"Clc1ccc(Cn2c(CNc3cc(-c4ccccc4)ccc3)cnc2)cc1 CHEMBL525178",
"Cc1cccc(-c2cc(NCc3cncn3Cc3ccccc3)ccc2)c1 CHEMBL482856",
"c1ccc(Cn2c(CNc3cc(-c4ccncc4)ccc3)cnc2)cc1 CHEMBL497405",
"Cc1c(-c2cc(NCc3cncn3Cc3ccccc3)ccc2)cccc1 CHEMBL520679",
"c1ccc(Cn2c(CNc3cc(-c4cnccc4)ccc3)cnc2)cc1 CHEMBL496834",
"O=[N+]([O-])c1ccc(Cn2c(CNc3cc(-c4ccccc4)ccc3)cnc2)cc1 CHEMBL497885",
"c1ccc(-c2ccc(Cn3c(CNc4ccc(-c5ccccc5)cc4)cnc3)cc2)cc1 CHEMBL497407",
*/
//"CCCCCCC1C23C4=c5c6c7c8c5c5c9c%10c%11c%12c(c%108)c8c7c7c%10c%13c%14c%15c%16c%17c%18c%19c%20c(c%21c%22c%23c(c9C(C25C[N+]1(C)C)C1c2c3c3c5c9c2-c(c%231)c(c%22%19)C%18C9C1(C5=C%13C(C43)c6%10)C%14%17C[N+](C)(C)C1CCCCCC)c%21%11)c%12c(c%16%20)c8c%157 CHEMBL439119"
/*
//WAS 8 sec test 
//NOW:
//Inspected Seeds      = 54250    Rejected by BestSize = 21360    Rejected by WrongComposition = 191( 134 generated)
//MatchCalls = 44607  MatchFound = 30334
//AtomCompareCalls = 103017020    BondCompareCalls = 25286076
//MCS : CCCC(NC(=O)CNC(=O)C(Cc(c)c)NC(=O)CNC(=O)CNC(=O)CNC(=O)C(C)NC=O)C(=O)NC(CCC)C(=O)NC(C)C 49 atoms, 48 bonds
//Time elapsed 20.19 seconds
//Time elapsed 35.65 seconds + FIX
//Time elapsed 123 seconds + FIX + can SMILES
"CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified QUERY",// CHEMBL439258
"CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccccc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258",// CHEMBL439258
"CC(C)CC(NC(=O)CNC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified",// CHEMBL439258

"CCCCC(NC(C(CCC(O)=O)NC(C(CC(C)C)NC(C(C(C)C)NC(=O)C(CCC(O)=O)NC(C(CCCN=C(N)N)NC(C(NC(=O)C(NC(C(NC(C1CCCNC(=O)CCC(N)C(=O)NC(CC(C)C)C(=O)NC(C(C)O)C(=O)N1)=O)Cc1c[nH]cn1)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)C(NC(C)C(NC(CCCN=C(N)N)C(NC(C)C(NC(CCC(O)=O)C(NC(CCC(N)=O)C(NC(CC(C)C)C(NC(C)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(C)C(NC(Cc1c[nH]cn1)C(NC(CO)C(NC(CC(N)=O)C(NC(CCCN=C(N)N)C(NC(CCCCN)C(NC(CC(C)C)C(NC(CCCC)C(NC(C(NC(C(C)CC)C(NC(C(N)=O)C(C)CC)=O)=O)CCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL438567",
"CCC(C)C(NC(CNC(=O)C(C)NC(=O)C(C)NC(C(Cc1nc[nH]c1)NC(C(CC(N)=O)NC(CNC(C(CO)NC(=O)C(C)NC(=O)C(CCC(N)=O)NC(C(NC(=O)C(NC(C(CCCN=C(N)N)NC(C(CCC(N)=O)NC(C(NC(C(CCCN=C(N)N)NC(CNC(C(CCC(N)=O)NC(C(CC(C)C)NC(C(C)N)=O)=O)=O)=O)=O)CC(C)C)=O)=O)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)=O)C(NC(CC(C)C)C(NC(C(O)C)C(NC(CCSC)C(O)=O)=O)=O)=O CHEMBL429374",
"CC(C)CC1NC(=O)C(CCCCN)NC(=O)C(Cc2ccc(O)cc2)NC(=O)CNC(=O)C2NC(=O)C(NC(C(C(C)C)NC(CNC(C3NC(=O)CC3)=O)=O)=O)CSSCC(C(O)=O)NC(=O)C3N(CCC3O)C(=O)C(Cc3ccccc3)NC(=O)C(CSSC2)NC1=O CHEMBL1076370",
*/

// SLOW
//test 45 #  Using CHEMBL551656 CHEMBL563796 CHEMBL561978 CHEMBL559467 CHEMBL550503 CHEMBL562866 CHEMBL552190 CHEMBL181547 CHEMBL359567 CHEMBL373316 
// 45 . 1 30 32 27.01 n12-C-c:c(-c:2:c:c2-C(-O)(-C-C)-C(-O-C-c:2:c:1=O)=O):n:c(:c:c:c-O):c(:c):c-C-C-C
//MCS: CCCCc(cc)c1Cn2-cc3COC(=O)[C](O)(CC):c3cc2-c1nccccO 30 atoms, 32 bonds
//Time elapsed 40.41 seconds
//MCS : CCCCc(cc)cCn1c(-cncccc)cc2c(COC(=O)[C]:2(O)CC)c-1=O 30 atoms, 31 bonds
//Time elapsed 9.70 seconds
"CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(OC)cc2c1C1CCCCC1 CHEMBL551656",
"CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(OC)cc2c1C1CCCC1 CHEMBL563796", //Q
"CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(OC)cc1c2C1CCCCCC1 CHEMBL561978",
"CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(OC)cc1c2C1CCC1 CHEMBL559467",
"CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(O)cc1c2C1CCCC1 CHEMBL550503",
"CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(O)cc2c1C1CCCCCC1 CHEMBL562866",
"CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(O)cc1c2C1CCCCC1 CHEMBL552190",
"CCC1(O)c2c(c(=O)n3c(c2)-c2nc4cc5c(cc4c(C4CCCC4)c2C3)OCO5)COC1=O CHEMBL181547",
"CCC1(O)c2c(c(=O)n3c(c2)-c2nc4c(c(C5CCCCC5)c2C3)cc2c(c4)OCO2)COC1=O CHEMBL359567",
"CCCc1c(OC)ccc2nc3c(c(CC)c21)Cn1c-3cc2c(c1=O)COC(=O)C2(O)CC CHEMBL373316",


/*
// # 190 TEST must < 0.27 sec
"COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)oc2cc1  CHEMBL1479679",
"COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1333382",
"Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)nc2cc1  CHEMBL1437584",
"COc1c(NC(=O)CSc2ccc(Cl)cc2)cc(-c2nc3ccccc3o2)cc1  CHEMBL1601350",
"Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008",
"Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)nc2cc1  CHEMBL1612903",
"COc1cc2nc(-c3cc(NC(=O)Cc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1316483",
"Cc1c(NC(=O)CSc2ccc(Cl)cc2)cccc1-c1nc2cc(Cl)ccc2o1  CHEMBL1568754",
"COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972",
"Cc1ccc(SCC(=O)Nc2cc(-c3nc4cc(C)ccc4o3)c(O)cc2)cc1  CHEMBL1611932",
//# 19 21 1.37 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1
//  19 21 2.36 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1 19 atoms, 21 bonds
*/
/*
//FIXED with //C1 - decrease excludeBonds   /// FAILED !!!!!!!!
//#  Using CHEMBL1515359 CHEMBL1590658 CHEMBL1447567 CHEMBL1384017 CHEMBL1456416 CHEMBL1308819 CHEMBL1703007 CHEMBL1707819 CHEMBL1500793 CHEMBL1334715
//32 . 1 31 33 0.82 S(-N1-C-C-O-C-C-1)(-c1:c:c:c(-N(-C-C)-C-C):c(-N-C(-C=C-c2:c:c:c:c:c:2)=O):c:1)(=O)=O
"O=C(Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCOCC1)C=Cc1ccc(Cl)cc1  CHEMBL1515359",
"c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)=O)cc1  CHEMBL1590658",
"Cc1ccc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)cc1  CHEMBL1447567",
"c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCCC2)=O)cc1  CHEMBL1384017",
"O=C(C=Cc1ccc(F)cc1)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1456416",
"c1cc(F)cc(C=CC(=O)Nc2c(N3CCCC3)ccc(S(N3CCOCC3)(=O)=O)c2)c1  CHEMBL1308819",
"CCN1CCN(c2ccc(S(N3CCOCC3)(=O)=O)cc2NC(=O)C=Cc2ccc(C)cc2)CC1  CHEMBL1703007",
"c1cc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)c([N+]([O-])=O)cc1  CHEMBL1707819",
"N#CC(=Cc1ccccc1)C(=O)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1500793",
"C(=Cc1ccc2c(c1)OCO2)C(Nc1cc(S(=O)(=O)N2CCOCC2)ccc1N1CCOCC1)=O  CHEMBL1334715",
// 31 31 0.05 sec MCS: CCOCCNS(=O)(=O)c1ccc(c(c1)NC(=O)C=Cc(c)cccc)N(CC)CC
// 31 32 0.15 sec MCS: CCOCCNS(=O)(=O)c1ccc(c(c1)NC(=O)C=Cc1ccccc1)N(CC)CC 31 atoms, 32 bonds
// 31 33 0.35 sec MCS: CCN(CC)c1ccc(cc1NC(=O)C=Cc1ccccc1)S(=O)(=O)N1CCOCC1
*/
    };

    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    {
//        mols.push_back(ROMOL_SPTR(SmartsToMol(smi[i]))); //it skips aromaticity
#if 0
// optional temporary AROMATIZATION - for this test ONLY
//(to avoid PRECONDITION Exception with Implicit Hs)
        std::auto_ptr<RWMol> mol(SmartsToMol(smi[i]));
        unsigned dummy;
        RDKit::MolOps::sanitizeMol(*mol, dummy, RDKit::MolOps::SANITIZE_ADJUSTHS|RDKit::MolOps::SANITIZE_SETAROMATICITY);
        std::string s = MolToSmiles(*mol);
        RWMol* m = SmartsToMol(s);
        RDKit::MolOps::sanitizeMol(*m, dummy, RDKit::MolOps::SANITIZE_ADJUSTHS|RDKit::MolOps::SANITIZE_SETAROMATICITY);
        mols.push_back(ROMOL_SPTR( m ));
        std::cout << "MOL : " << s <<"\n";
//------
#else
        std::string id;
        mols.push_back(ROMOL_SPTR(SmartsToMol( getSmilesOnly(smi[i], &id) )));
        std::cout << id << "\n";
#endif
    }

    t0 = nanoClock();

    MCSParameters p;
    //p.Threshold = 0.7;
    //p.Timeout = 9;
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}


void testRing1()
{
    std::cout << "\ntestRing1()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
        "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
//      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12", // original molecule
        "COCc1cnc(C(=O)OC(C)C)c2[nH]ccc(Oc4ccc(Cl)cc4)cccc12",   // ring 3 removed
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));   // with RING INFO

    MCSParameters p;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly   = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}


void test504()
{
    std::cout << "\ntest504()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
        //"O=C(-N-c1:c:c:c:c:c:1)-N-C-C-C-N(-C)-C-C-C-C-C-N-C(-C1-C-C-1-c1:c:c:c(-Cl):c(-Cl):c:1)=O", // python MCS == SMARTS !!!
        //"Clc1ccc(NC(=O)NC2CCN(CCCCCNC(=O)C3CC3[c:3]3[cH:2]cc(Cl)c(Cl)[cH:1]3)C2)cc1",
        //TEST 504
"C(CCNC(C1CC1[c:1]1[c:2]c(Cl)c(Cl)c[c:3]1)=O)CCN1CCC(NC(Nc2ccc(Cl)cc2)=O)C1 CHEMBL545864",  // - QUERY
//"C(CCNC(C1CC1[c:1]1[c:2]c(Cl)c(Cl)c[c:3]1)=O)CCN1CCC(NC(Nc2ccccc2)=O)C1 CHEMBL545864",  // - QUERY - Cl:30

        "FC(F)(F)c1cc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl CHEMBL528228",
        "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)ccc1Cl CHEMBL525875",
        "Fc1ccc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F CHEMBL527277",
        "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl CHEMBL537333",
        "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)cc1C(F)(F)F CHEMBL588077",
        "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)C2)=O)cc1 CHEMBL525307",
        "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F CHEMBL581847",
        "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)CC2)=O)cc1 CHEMBL579547",

//        "C(CCNC(C1CC1c1cc(Cl)c(Cl)cc1)=O)CCN1CCC(NC(Nc2ccc(Cl)cc2)=O)C1 CHEMBL545864",  // - QUERY
//        "Clc1ccc(NC(=O)NC2CCN(CCCCCNC(=O)C3CC3c3ccc(Cl)c(Cl)c3)C2)cc1",   // - the same QUERY with RIGHT MCS !!!
//        "Clc1ccc(NC(=O)NC2CCN(CCCCCNC(=O)C3CC3[c:1]3[cH:2]cc(Cl)c(Cl)[cH:3]3)C2)cc1",   // - the same QUERY with Atom MAP and RIGHT MCS !!!

        "N#Cc1cccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)c1 CHEMBL529994",
    };
    RWMol* qm = SmartsToMol( getSmilesOnly(smi[0]) );
    unsigned nq = qm->getNumAtoms();
    for(size_t ai = 0; ai < nq; ai++)
    {
        Atom* atom = qm->getAtomWithIdx(ai);
        atom->setProp("molAtomMapNumber", (int)ai);
    }
    std::cout<<"Query +MAP "<< MolToSmiles(*qm) <<"\n";
    mols.push_back(ROMOL_SPTR(qm));   // with RING INFO
    for(int i=1; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmartsToMol( getSmilesOnly(smi[i]) )));   // with RING INFO
    MCSParameters p;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void test18()
{
    std::cout << "\ntest18()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
    //TEST 18
    "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991", //-- QUERY
    "Cc1ccc(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682990",
    "Cc1cccnc1CN(C(C)c1ccccn1)CCCCN CHEMBL1682998",
    "CC(N(CCCCN)Cc1c(N)cccn1)c1ccccn1 CHEMBL1682987",
    "Cc1cc(C)c(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682992",
    "Cc1cc(C(C)N(CCCCN)Cc2c(C)cccn2)ncc1 CHEMBL1682993",
    "Cc1nc(C(C)N(CCCCN)Cc2nc3c([nH]2)cccc3)ccc1 CHEMBL1682878",
    "CC(c1ncccc1)N(CCCCN)Cc1nc2c([nH]1)cccc2 CHEMBL1682867",
    "CC(N(CCCCN)Cc1c(C(C)(C)C)cccn1)c1ccccn1 CHEMBL1682989",
    "CC(N(CCCCN)Cc1c(C(F)(F)F)cccn1)c1ccccn1 CHEMBL1682988",
    //# 18 .  20 20 0.04 sec. Python MCS: CC(c1ccccn1)N(CCCCN)Ccnccc
    };
    RWMol* qm = SmartsToMol( getSmilesOnly(smi[0]) );
    unsigned nq = qm->getNumAtoms();
    for(size_t ai = 0; ai < nq; ai++)
    {
        Atom* atom = qm->getAtomWithIdx(ai);
        atom->setProp("molAtomMapNumber", (int)ai);
    }
    std::cout<<"Query +MAP "<< MolToSmiles(*qm) <<"\n";
    mols.push_back(ROMOL_SPTR(qm));   // with RING INFO
    for(int i=1; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmartsToMol( getSmilesOnly(smi[i]) )));   // with RING INFO
    MCSParameters p;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testThreshold()
{
    std::cout << "\ntestThreshold()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
        "CCC", "CCCO", "CCCN", "CC",
        "CCC", "CCCO", "CCCN", "CC",
        "CCC", "CC",
/*
      "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12",

//        "CC(c1ccccn1)N(C)Ccnccc" // short "MCS" for test 18
    //TEST 18
    "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991", //-- QUERY
    "Cc1ccc(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682990",
    "Cc1cccnc1CN(C(C)c1ccccn1)CCCCN CHEMBL1682998",
    "CC(N(CCCCN)Cc1c(N)cccn1)c1ccccn1 CHEMBL1682987",
    "Cc1cc(C)c(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682992",
    "Cc1cc(C(C)N(CCCCN)Cc2c(C)cccn2)ncc1 CHEMBL1682993",
    "Cc1nc(C(C)N(CCCCN)Cc2nc3c([nH]2)cccc3)ccc1 CHEMBL1682878",
    "CC(c1ncccc1)N(CCCCN)Cc1nc2c([nH]1)cccc2 CHEMBL1682867",
    "CC(N(CCCCN)Cc1c(C(C)(C)C)cccn1)c1ccccn1 CHEMBL1682989",
    "CC(N(CCCCN)Cc1c(C(F)(F)F)cccn1)c1ccccn1 CHEMBL1682988",
    //# 18 .  20 20 0.04 sec MCS: CC(c1ccccn1)N(CCCCN)Ccnccc
*/
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmartsToMol( getSmilesOnly(smi[i]) )));
    findMCS(mols);
    MCSParameters p;
    p.Threshold = 0.7;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testSimpleFast()
{
    std::cout << "\ntestSimpleFast()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
   // SHORT TEST for 26 bonds.
// Python MCS = 26 bonds : COCc1cncc(c1):n:c1cccc(Oc2ccc(Cl)cc2)c1
// MCS 26: COCc1c-ncc(c1)nc1cccc(c1)Oc1ccc(Cl)cc1 24 atoms, 26 bonds
///            "COCC1=C(N=CC2=C1C1=C(OC3=CC=C(Cl)C=C3)C=CC=C1N2)C(=O)OC(C)C",
///            "COCC1=CN=C(C(=O)OC(C)C)C2=C1C1=CC=C(OC3=CC=C(Cl)C=C3)C=C1N2",
// The SAME, but pre AROMATIZATED (else PRECONDITION Exception with Implicit Hs / 16 bonds only)
        "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
        "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12",
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
    MCSParameters p;
    p.Threshold = 0.7;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly   = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testSimple()
{
    std::cout << "\ntestSimple()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] =
    {
        // LONG TIME TEST for performance analisis (about 30 sec)
        "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified QUERY",// CHEMBL439258
        "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccccc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258",// CHEMBL439258
        "CC(C)CC(NC(=O)CNC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified",// CHEMBL439258

        "CCCCC(NC(C(CCC(O)=O)NC(C(CC(C)C)NC(C(C(C)C)NC(=O)C(CCC(O)=O)NC(C(CCCN=C(N)N)NC(C(NC(=O)C(NC(C(NC(C1CCCNC(=O)CCC(N)C(=O)NC(CC(C)C)C(=O)NC(C(C)O)C(=O)N1)=O)Cc1c[nH]cn1)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)C(NC(C)C(NC(CCCN=C(N)N)C(NC(C)C(NC(CCC(O)=O)C(NC(CCC(N)=O)C(NC(CC(C)C)C(NC(C)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(C)C(NC(Cc1c[nH]cn1)C(NC(CO)C(NC(CC(N)=O)C(NC(CCCN=C(N)N)C(NC(CCCCN)C(NC(CC(C)C)C(NC(CCCC)C(NC(C(NC(C(C)CC)C(NC(C(N)=O)C(C)CC)=O)=O)CCC(O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL438567",
        "CCC(C)C(NC(CNC(=O)C(C)NC(=O)C(C)NC(C(Cc1nc[nH]c1)NC(C(CC(N)=O)NC(CNC(C(CO)NC(=O)C(C)NC(=O)C(CCC(N)=O)NC(C(NC(=O)C(NC(C(CCCN=C(N)N)NC(C(CCC(N)=O)NC(C(NC(C(CCCN=C(N)N)NC(CNC(C(CCC(N)=O)NC(C(CC(C)C)NC(C(C)N)=O)=O)=O)=O)=O)CC(C)C)=O)=O)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)=O)C(NC(CC(C)C)C(NC(C(O)C)C(NC(CCSC)C(O)=O)=O)=O)=O CHEMBL429374",
        "CC(C)CC1NC(=O)C(CCCCN)NC(=O)C(Cc2ccc(O)cc2)NC(=O)CNC(=O)C2NC(=O)C(NC(C(C(C)C)NC(CNC(C3NC(=O)CC3)=O)=O)=O)CSSCC(C(O)=O)NC(=O)C3N(CCC3O)C(=O)C(Cc3ccccc3)NC(=O)C(CSSC2)NC1=O CHEMBL1076370",
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
    MCSParameters p;
    //p.Threshold = 0.5;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly   = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testCmndLineSMILES(int argc, const char* argv[])
{
    std::vector<ROMOL_SPTR> mols;
    for(int i=1; i<argc; i++)
        mols.push_back(ROMOL_SPTR(SmartsToMol(argv[i])));
    MCSResult res = findMCS(mols);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}


double testFileSDF(const char* test)
{
    std::cout << "\ntestFileSDF(): " << test << "\n";
    std::vector<ROMOL_SPTR> mols;
    std::string fn(test);
    RDKit::SDMolSupplier suppl(fn);
    while(!suppl.atEnd())
    {
        ROMol *m=suppl.next();
        if(m)
            mols.push_back(ROMOL_SPTR(m));
    }
    MCSParameters p;
//    p.MaximizeBonds = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    double t = (nanoClock() - t0) / 1000000.;
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
    return t;
}

void testFileSMILES(const char* test)
{
    std::vector<ROMOL_SPTR> mols;
    char smiles[4096];
    unsigned n=0;
    FILE* f = fopen(test, "rt");
    std::cout<<"Loading SMILES ... \n";
    while(fgets(smiles, sizeof(smiles), f))
    {
        std::cout<<"\rLine: "<< ++n <<" ";
        if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) // commented to skip
            if(strlen(smiles) > 92) // minimal query size !!!
                mols.push_back(ROMOL_SPTR(SmartsToMol(getSmilesOnly(smiles))));
    }
    fclose(f);
    printTime();
    std::cout<<"FIND MCS in "<<mols.size()<<" molecules.\n\n";
    t0 = nanoClock();
    MCSParameters p;
    //p.Threshold = 0.7;
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

//====================================================================================================
//====================================================================================================

void testGregSDFFileSetFiltered()
{
    const std::string sdf_dir("Greg'sComparision/data/filtered/");
    const char* sdf[] =
    {
    //fmcs: 	0.11		27	29	beta2_adrenergic_aid465635.filtered.sdf	O=C(O)C1:C:C:C(C2:C:C:C(CCNCC(O)C3:C:C:C:C:C:3):C:C:2):C:C:1
    "beta2_adrenergic_aid465635.filtered.sdf",
    //fmcs: 	1.90		22	23	beta2_adrenergic_aid578239.filtered.sdf	CNC1:C:C:C(CCNCC(O)COC2:C:C:C:C:C:2):C:C:1
    "beta2_adrenergic_aid578239.filtered.sdf",
    //fmcs: 	1.88		22	23	beta2_adrenergic_aid578240.filtered.sdf	CNC1:C:C:C(CCNCC(O)COC2:C:C:C:C:C:2):C:C:1
    "beta2_adrenergic_aid578240.filtered.sdf",
    //fmcs: 	0.07		17	18	beta2_adrenergic_aid759384.filtered.sdf	CCNCC(O)C1:C:C:C(O):C2:N:C(=O):[SH]:C:2:1
    "beta2_adrenergic_aid759384.filtered.sdf",
    //fmcs: 	1.01		25	25	d3_aid329485.filtered.sdf	C:C:C:C:C:CC(=O)NCCCCN1CCN(C:C:C:C:C:C)CC1
    "d3_aid329485.filtered.sdf",
    //fmcs: 	0.13		20	21	d3_aid367038.filtered.sdf	C:C:C:C:N:C:CCN1CCN(C2:C:C:C:C:C:2)CC1
    "d3_aid367038.filtered.sdf",
    //fmcs: 	0.78		27	29	d3_aid563533.filtered.sdf	O=C(C:C:C1:C:C:C:C:C:1)NCCCCN1CCN(C2:C:C:C:C:C:2)CC1
    "d3_aid563533.filtered.sdf",
    //fmcs: 	0.09		14	14	d3_aid563770.filtered.sdf	CCCN(CC)C(=O)C1:C:C:C:C:C:1
    "d3_aid563770.filtered.sdf",
    //fmcs: 	0.14		22	24	d3_aid578980.filtered.sdf	C(:C:N:C1:C:C:N:C:N:1)CN1CCN(C2:C:C:C:C:C:2)CC1
    "d3_aid578980.filtered.sdf",
    //fmcs: 	0.20		15	17	d3_aid58783.filtered.sdf	CCN1CC2COC3:C:C:C:C:C:3C2C1
    "d3_aid58783.filtered.sdf",
    //fmcs: 	0.08		14	14	d3_aid62278.filtered.sdf	C:C(CNCC):N:CC1:C:C:C:C:C:1
    "d3_aid62278.filtered.sdf",
    //fmcs: 	0.05		7	7	d3_aid62281.filtered.sdf	CC1:C:C:C:C:C:1
    "d3_aid62281.filtered.sdf",
    //fmcs: 	33.13		26	27	d3_aid62457.filtered.sdf	CC(=O)NCCCCN1CCC2:C:C:C(OS(=O)(=O)C(F)(F)F):C:C:2C1
    "d3_aid62457.filtered.sdf",
    //fmcs: 	0.08		14	15	d3_aid62774.filtered.sdf	CCN1CCCC(C2:C:C:C:C:C:2)C1
    "d3_aid62774.filtered.sdf",
    //fmcs: 	0.33		14	14	d3_aid642590.filtered.sdf	C:C:C:C:C:C:CCN1CCCCC1
    "d3_aid642590.filtered.sdf",
    };

    double totalT=0.;
    for(int i=0; i<sizeof(sdf)/sizeof(sdf[0]); i++)
        totalT += testFileSDF((sdf_dir+sdf[i]).c_str());
    printf("\nTOTAL Time elapsed %.2f seconds\n================================================\n", totalT);
}
//====================================================================================================
//====================================================================================================


int main(int argc, const char* argv[])
{ 
    RDKit::FMCS::ConsoleOutputEnabled = true;
// use maximum CPU resoures to increase time measuring accuracy and stability in multi process environment
#ifdef WIN32
    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
    SetThreadPriority(GetCurrentThread (), THREAD_PRIORITY_HIGHEST );
#else
    setpriority(PRIO_PROCESS, getpid(), -20); 
#endif

    T0 = nanoClock();
    t0 = nanoClock();

#ifdef xxWIN32  // brief test set for testing and issue investigation
{
    //test18();
    //testThreshold();
    testSimple();

    //testGregSDFFileSetFiltered();
    //test18();
    //return 0;

//    testRing1();
//    return 0;

//    testSimpleFast();
//    testSimple();
//    return 0;

//    testGregSDFFileSetFiltered();
    //return 0;

    std::vector<unsigned> tc;
//    tc.push_back(10);
//    tc.push_back(992);
//992 PYTHON 20 21 N1(-C-C=C(-c2:c(:c:c:c):c:n:c:2)-C-C-1)-C-C-C-C-C-C
//992 . 1    27 28 1.10 CNcc(CCCCCN1CCC(=CC1)c1cncc1ccc)cccc
//now        25 26

/*
    tc.push_back(33);
    tc.push_back(59);
      tc.push_back(124);

    tc.push_back(345);
    tc.push_back(605);
    tc.push_back(619);
*/
//    testFileMCSB(argv[1], 300, tc);
    return 0;
}
#endif

    if(3 == argc && '-' == argv[1][0])
     switch(argv[1][1]) // ./test -s|m|b <filename with test files list>
    {
     case 's':  // smiles files list
         {
            char test[256];
            FILE* f = fopen(argv[2], "rt");
            while(fgets(test, sizeof(test), f))
                testFileSMILES(test);
            fclose(f);
         }
         break;
     case 'm':  // SDF mol files list
         {
            char test[256];
            FILE* f = fopen(argv[2], "rt");
            while(fgets(test, sizeof(test), f))
                testFileSDF(test);
            fclose(f);
         }
         break;
     case 'b':
         {
            std::vector<unsigned> tc; // empty -> all
            testFileMCSB(argv[2], 30, tc);    // .mcsb
         }
         break;
     default:
         break;
    }
    else if(2 == argc)   // .sdf /.smi file
    {
        if(0==strcmp(argv[1]+strlen(argv[1])-4, ".smi"))
            testFileSMILES(argv[1]);// .smi
        else if(0==strcmp(argv[1]+strlen(argv[1])-4, ".sdf"))
            testFileSDF(argv[1]);   // .sdf
        else if(0==strcmp(argv[1]+strlen(argv[1])-4, "mcsb"))
            testFileMCSB(argv[1], 30);  // .mcsb
        else
            printf("UNKNOWN File Extention.\n");
    }
    else if(argc > 1+2)
        testCmndLineSMILES(argc, argv);
    else
    {
        testSimpleFast();
        testSimple();
//        test1Basics();
//        testGregSDFFileSetFiltered();
    }
//  BOOST_LOG(rdInfoLog) << "*******************************************************\n";

    unsigned long long t1 = nanoClock();
    double sec = double(t1-T0) / 1000000.;
    printf("TOTAL Time elapsed %.2f seconds\n", sec);

    return 0;
}


