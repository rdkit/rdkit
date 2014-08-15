// $Id: testFMCS.cpp $
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
#include "../../RDKitBase.h"
#include "../../FileParsers/FileParsers.h" //MOL single molecule !
#include "../../FileParsers/MolSupplier.h" //SDF
#include "../../SmilesParse/SmilesParse.h"
#include "../../SmilesParse/SmilesWrite.h"
#include "../../SmilesParse/SmartsWrite.h"
#include "../FMCS.h"
#include "../DebugTrace.h" //#ifdef VERBOSE_STATISTICS_ON

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

void printTime() {
    unsigned long long t1 = nanoClock();
    double sec = double(t1-t0) / 1000000.;
    printf("Time elapsed %.3lf seconds\n", sec);
    t0 = nanoClock();
}

std::string getSmilesOnly(const char* smiles, std::string* id=0) { // remove label, because RDKit parse FAILED
    const char* sp = strchr(smiles,' ');
    unsigned n = (sp ? sp-smiles+1 : strlen(smiles));
    if(id)
        *id = std::string(smiles+n);
    return std::string(smiles, n);
}

std::string getSmilesOnlyTxt(const char* smiles, std::string* id=0) { // remove label from string like "CHEMBL90218 NS(=O)(=O)c1ccc(NC(=O)c2cccc(C(=O)O)n2)c(Cl)c1"
    const char* sp = strchr(smiles,' ');
    if(sp && '\0'!=*sp) {
        if(id)
            *id = std::string(smiles, sp-smiles);
        sp++;
        size_t i;
        for(i = strlen(sp)-1; sp[i] < ' ' && i>=0; i--)
            ;
        return std::string(sp, i+1);
    } else
        return smiles;
}


std::string getSmilesOnlyChEMBL(const char* smiles, std::string* id=0) { // remove label, because RDKit parse FAILED
    const char* sp = strchr(smiles, '\t');
    if(sp) {
        unsigned n = (sp ? sp-smiles+1 : strlen(smiles));
        if(id)
            *id = std::string(smiles, n);
        sp = strchr(++sp,'\t');
        if(sp)
            sp++;
    }
    return std::string(sp);
}

//====================================================================================================

MCSParameters p;


void testFileMCSB(const char* test, unsigned timeout=30, std::vector<unsigned> test_N=std::vector<unsigned>()) { // optional list of some tests for investigation
    p.Verbose = false;

    std::vector<ROMOL_SPTR> mols; // IT CAN OCCUPY A LOT OF MEMORY. store SMILES only to reduce memory usage.
    char str [4096];
    std::string molFile, id;
    std::map<std::string, size_t>         molIdMap;
    std::vector<std::string>            smilesList;
    std::list< std::vector<std::string> > testCase;
    std::string referenceOutFile(test);
    referenceOutFile += ".REF.out";
    std::string outFile(test);
    if(!test_N.empty()) {
        if(1==test_N.size())
            sprintf(str,".%u.out", test_N[0]);
        else
            sprintf(str,".%u-%u.out", test_N[0], test_N.back());
        outFile += str;
    } else {
        outFile += ".Cpp.out";
    }

    std::string outSmilesFile(test);
    outSmilesFile += ".smiles";
    unsigned n=0, passed=0, failed=0, failed_less=0, timedout=0;
    double secTotal = 0.;

    std::vector<MCSResult>  referenceResults;
    std::vector<float>      referenceResultsTime;

    FILE* f = fopen(referenceOutFile.c_str(), "rt");
    if(!f)
        perror("Could not open reference test result file");
    else {
        std::cout<<"Loading reference test results ... \n";
        while(fgets(str, sizeof(str), f))
            if('#' != str[0]) {
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
    if(!f) {
        perror("Could not open test case list MCSB file");
        exit(1);
    }
    {
        std::cout<<"Loading MCSB test list ... \n";
        if(fgets(str, sizeof(str), f))
            if(fgets(str, sizeof(str), f)) {
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
            if('#' != str[0]) {
                //str= "1 CHEMBL526291 CHEMBL498211 ..."
                char name [256];
                unsigned nn, len;
                n++;
                testCase.push_back(std::vector<std::string>());
                sscanf(str, "%u%n", &nn, &len);
                while('\0'!=*(str+len) && 1==sscanf(str+len, "%s%n", name, &nn)) {
                    len += nn;
                    testCase.back().push_back(std::string(name));
                }
            }
        std::cout<<n<<" Test cases loaded\n";
    }
    fclose(f);

    f = fopen(molFile.c_str(), "rt");
    if(!f) {
        perror("Could not open molecules file");
        exit(1);
    }
    std::cout<<"Loading SMILES ... \n";
    n = 0;
    while(fgets(str, sizeof(str), f)) {
        std::cout<<"\rLine: "<< ++n <<" ";
        if('#' != str[0] && ' ' != str[0] && '/' != str[0]) { // commented to skip
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
    if(!f) {
        perror("Could not create output file");
        exit(1);
    }
    setvbuf(f , 0, _IOFBF, 4*1024);
    FILE* fs = 0;  //fopen(outSmilesFile.c_str(), "wt");
    if(fs)
        setvbuf(fs, 0, _IOFBF, 4*1024);
#ifdef xxVERBOSE_STATISTICS_ON
    FILE* ft  = fopen((outFile+".stat.csv").c_str(), "wt");
    setvbuf(ft, 0, _IOFBF,  4*1024); // small file
    if(ft)
        fprintf(ft, "N; Status; dAtoms; dBonds; t(sec); ref.t; Seed; MatchCall; AtomCmp; BondCmp\n"); //CSV Header
#endif
    n = 0;
    p.Timeout = timeout;
    p.Threshold = 1.0;

    fprintf(f, "#software RDKit C++ FMCS \n#options  timeout=%u threshold=%g\n", p.Timeout, p.Threshold);
    std::cout<<"Perform test cases ... \n";
    for(std::list< std::vector<std::string> >::const_iterator tc = testCase.begin(); tc != testCase.end(); tc++, n++) {
        if(!test_N.empty() && test_N.end() == std::find(test_N.begin(), test_N.end(), n+1))
            continue;

        std::cout<<"\rTest: "<< n+1 <<" ";
        if(!test_N.empty()) // test case is listed
            std::cout<<"\n";

        std::vector<ROMOL_SPTR> tcmols;
        fprintf(f, "# %u Using ", n+1);
        if(fs)
            fprintf(fs, "\n//TEST %u\n", n+1);
        for(std::vector<std::string>::const_iterator mid = tc->begin(); mid != tc->end(); mid++) {
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

        if(!referenceResults.empty()) {
            fprintf(f, "# %u REFCMP: time %s %.2f sec.\n", n+1, fabs(referenceResultsTime[n] - sec)<sec/25.?"EQUAL"
                    :(referenceResultsTime[n] < sec ? (sec<.7?"slow":"SLOW")
                          :(sec<.7 && referenceResultsTime[n]<0.7 ? "fast":"FAST")), referenceResultsTime[n]);
            if(!referenceResults[n].Canceled) { // && !res.Canceled)
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
                else {
                    if(referenceResults[n].NumBonds > res.NumBonds)
                        failed_less++;
                    failed++;
                }
            } else
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
    if(f)
        fclose(f);
    if(fs)
        fclose(fs);
#ifdef xxVERBOSE_STATISTICS_ON
    if(ft)
        fclose(ft);
#endif
    printTime();
}

//=========================================================================

void test504() {
    std::cout << "\ntest504()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] = {
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

        "N#Cc1cccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)c1 CHEMBL529994",
    };
    RWMol* qm = SmilesToMol( getSmilesOnly(smi[0]) );
    unsigned nq = qm->getNumAtoms();
    for(size_t ai = 0; ai < nq; ai++) {
        Atom* atom = qm->getAtomWithIdx(ai);
        atom->setProp("molAtomMapNumber", (int)ai);
    }
    std::cout<<"Query +MAP "<< MolToSmiles(*qm) <<"\n";
    mols.push_back(ROMOL_SPTR(qm));   // with RING INFO
    for(int i=1; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));   // with RING INFO
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

std::string testChEMBL_Txt(const char* test, double th=1.0, const char* csv="chembl_II_sets.C++.res.csv") {
    std::cout << "\ntestChEMBL_Txt() "<<test<<"\n";
    std::vector<ROMOL_SPTR> mols;
    char smiles[4096];
    unsigned n=0;
    FILE* f = fopen(test, "rt");
    if(!f) {
        perror("fopen testChEMBL_Txt()");
        return "";
    }
    char testsmi[512];
    strcpy(testsmi, test);
    strcpy(testsmi + strlen(test)-3, "smi");
    FILE* fsmi = fopen(testsmi, "wt");
    FILE* fres = fopen(csv, "at");
    while(fgets(smiles, sizeof(smiles), f)) {
        if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) { // commented to skip
            std::string id;
            std::string smi = getSmilesOnlyTxt(smiles, &id);
            fprintf(fsmi, "%s\n", (smi+" "+ id).c_str());
            mols.push_back(ROMOL_SPTR(SmilesToMol(smi)));
        }
    }
    fclose(f);
    fclose(fsmi);

    p.Threshold = th;

    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    unsigned long long tc1 = nanoClock();
    double sec = double(tc1-t0) / 1000000.;

    std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
    fprintf(fres, "%s;%u;%.1f;%.4f;%u;%u;%s;", test, mols.size(), th, sec, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());

    p.BondTyper = MCSBondCompareOrderExact;
    t0 = nanoClock();
    res = findMCS(mols, &p);
    tc1 = nanoClock();
    sec = double(tc1-t0) / 1000000.;
    p.BondTyper = MCSBondCompareOrder;

    std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
    fprintf(fres, "%.1f;%.4f;%u;%u;%s\n", sec, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());

    fclose(fres);
    return testsmi;
}

void testChEMBL_TxtALL_chembl_II_sets(double th=1.0) {
    FILE* fres = fopen("chembl_II_sets.C.res.csv", "wt");
    fprintf(fres,"test;Nmols;threshold;t,sec;atoms;bonds;mcs;E t,sec;E atoms;E bonds;E mcs\n");
    fclose(fres);
    const char* test[] = {
        "Target_no_100_15113.txt",
        "Target_no_100_16613.txt",
        "Target_no_100_30745.txt",
        "Target_no_100_30861.txt",
        "Target_no_100_31495.txt",
        "Target_no_100_34869.txt",
        "Target_no_100_36158.txt",
        "Target_no_100_36559.txt",
        "Target_no_100_37434.txt",
        "Target_no_100_39147.txt",
        "Target_no_100_43754.txt",
        "Target_no_100_44011.txt",
        "Target_no_100_44158.txt",
        "Target_no_100_48076.txt",
        "Target_no_100_48961.txt",
        "Target_no_100_49648.txt",
        "Target_no_100_49703.txt",
        "Target_no_100_51128.txt",
        "Target_no_100_51343.txt",
        "Target_no_100_51842.txt",
        "Target_no_100_54554.txt",
        "Target_no_100_55749.txt",
        "Target_no_10188_17256.txt",
        "Target_no_10188_20244.txt",
        "Target_no_10188_20245.txt",
        "Target_no_10188_30149.txt",
        "Target_no_10188_30499.txt",
        "Target_no_10188_36186.txt",
        "Target_no_10188_39262.txt",
        "Target_no_10188_43966.txt",
        "Target_no_10188_45279.txt",
        "Target_no_10188_48029.txt",
        "Target_no_10188_49064.txt",
        "Target_no_10188_49801.txt",
        "Target_no_10188_50681.txt",
        "Target_no_10188_54039.txt",
        "Target_no_10188_58265.txt",
        "Target_no_10188_58358.txt",
        "Target_no_10188_58866.txt",
        "Target_no_10188_60810.txt",
        "Target_no_10193_21358.txt",
        "Target_no_10193_21597.txt",
        "Target_no_10193_30111.txt",
        "Target_no_10193_30349.txt",
        "Target_no_10193_31746.txt",
        "Target_no_10193_31822.txt",
        "Target_no_10193_35220.txt",
        "Target_no_10193_37603.txt",
        "Target_no_10193_37964.txt",
        "Target_no_10193_39508.txt",
        "Target_no_10193_42164.txt",
        "Target_no_10193_42577.txt",
        "Target_no_10193_44072.txt",
        "Target_no_10193_45522.txt",
        "Target_no_10193_46700.txt",
        "Target_no_10193_46719.txt",
        "Target_no_10193_46744.txt",
        "Target_no_10193_48517.txt",
        "Target_no_10193_49087.txt",
        "Target_no_10193_50034.txt",
        "Target_no_10193_51475.txt",
        "Target_no_10193_53483.txt",
        "Target_no_10193_53506.txt",
        "Target_no_10193_54598.txt",
        "Target_no_10193_55614.txt",
        "Target_no_10193_55737.txt",
        "Target_no_10193_57090.txt",
        "Target_no_10193_57705.txt",
        "Target_no_10193_57785.txt",
        "Target_no_10193_58463.txt",
        "Target_no_10193_58700.txt",
        "Target_no_10193_61003.txt",
        "Target_no_10193_61349.txt",
        "Target_no_10260_38526.txt",
        "Target_no_10260_50778.txt",
        "Target_no_10260_52179.txt",
        "Target_no_10260_54285.txt",
        "Target_no_10280_30224.txt",
        "Target_no_10280_31222.txt",
        "Target_no_10280_37349.txt",
        "Target_no_10280_37369.txt",
        "Target_no_10280_37465.txt",
        "Target_no_10280_37541.txt",
        "Target_no_10280_38371.txt",
        "Target_no_10280_44314.txt",
        "Target_no_10280_44578.txt",
        "Target_no_10280_45765.txt",
        "Target_no_10280_46281.txt",
        "Target_no_10280_47928.txt",
        "Target_no_10280_50624.txt",
        "Target_no_10280_51341.txt",
        "Target_no_10280_51662.txt",
        "Target_no_10280_51767.txt",
        "Target_no_10280_53081.txt",
        "Target_no_10280_53396.txt",
        "Target_no_10280_53407.txt",
        "Target_no_10280_58823.txt",
        "Target_no_10280_59190.txt",
        "Target_no_10280_60672.txt",
        "Target_no_10280_61176.txt",
        "Target_no_10434_20905.txt",
        "Target_no_10434_30074.txt",
        "Target_no_10434_30175.txt",
        "Target_no_10434_30423.txt",
        "Target_no_10434_31555.txt",
        "Target_no_10434_34878.txt",
        "Target_no_10434_35926.txt",
        "Target_no_10434_42168.txt",
        "Target_no_107_20859.txt",
        "Target_no_107_30060.txt",
        "Target_no_107_30061.txt",
        "Target_no_107_30606.txt",
        "Target_no_107_30866.txt",
        "Target_no_107_31249.txt",
        "Target_no_107_34921.txt",
        "Target_no_107_35188.txt",
        "Target_no_107_36636.txt",
        "Target_no_107_37466.txt",
        "Target_no_107_39202.txt",
        "Target_no_107_43952.txt",
        "Target_no_107_49345.txt",
        "Target_no_107_49591.txt",
        "Target_no_107_49854.txt",
        "Target_no_107_51416.txt",
        "Target_no_107_53537.txt",
        "Target_no_107_58879.txt",
        "Target_no_107_60679.txt",
        "Target_no_107_61269.txt",
        "Target_no_108_20600.txt",
        "Target_no_108_20859.txt",
        "Target_no_108_30060.txt",
        "Target_no_108_30606.txt",
        "Target_no_108_30866.txt",
        "Target_no_108_31558.txt",
        "Target_no_108_35188.txt",
        "Target_no_108_37466.txt",
        "Target_no_108_43875.txt",
        "Target_no_108_48911.txt",
        "Target_no_108_49591.txt",
        "Target_no_108_49854.txt",
        "Target_no_108_57887.txt",
        "Target_no_108_60679.txt",
        "Target_no_108_6857.txt",
        "Target_no_10980_15478.txt",
        "Target_no_10980_20853.txt",
        "Target_no_10980_30840.txt",
        "Target_no_10980_30862.txt",
        "Target_no_10980_30918.txt",
        "Target_no_10980_30994.txt",
        "Target_no_10980_31014.txt",
        "Target_no_10980_31508.txt",
        "Target_no_10980_31623.txt",
        "Target_no_10980_35003.txt",
        "Target_no_10980_35004.txt",
        "Target_no_10980_35281.txt",
        "Target_no_10980_35957.txt",
        "Target_no_10980_36520.txt",
        "Target_no_10980_36641.txt",
        "Target_no_10980_38586.txt",
        "Target_no_10980_44431.txt",
        "Target_no_10980_45093.txt",
        "Target_no_10980_51302.txt",
        "Target_no_10980_5247.txt",
        "Target_no_10980_52937.txt",
        "Target_no_10980_57198.txt",
        "Target_no_10980_5739.txt",
        "Target_no_10980_6535.txt",
        "Target_no_11140_20205.txt",
        "Target_no_11140_20986.txt",
        "Target_no_11140_30091.txt",
        "Target_no_11140_34698.txt",
        "Target_no_11140_34705.txt",
        "Target_no_11140_37038.txt",
        "Target_no_11140_37272.txt",
        "Target_no_11140_37447.txt",
        "Target_no_11140_38909.txt",
        "Target_no_11140_38953.txt",
        "Target_no_11140_39472.txt",
        "Target_no_11140_45825.txt",
        "Target_no_11140_47179.txt",
        "Target_no_11140_49028.txt",
        "Target_no_11140_49860.txt",
        "Target_no_11140_51984.txt",
        "Target_no_11140_53395.txt",
        "Target_no_11140_53826.txt",
        "Target_no_11140_54121.txt",
        "Target_no_11140_54530.txt",
        "Target_no_11140_55738.txt",
        "Target_no_11140_57725.txt",
        "Target_no_11140_58590.txt",
        "Target_no_11365_20432.txt",
        "Target_no_11365_30331.txt",
        "Target_no_11365_37956.txt",
        "Target_no_11365_48114.txt",
        "Target_no_11365_58566.txt",
        "Target_no_11489_34731.txt",
        "Target_no_11489_34818.txt",
        "Target_no_11489_37318.txt",
        "Target_no_11489_37339.txt",
        "Target_no_11489_37564.txt",
        "Target_no_11489_39504.txt",
        "Target_no_11489_39561.txt",
        "Target_no_11489_46082.txt",
        "Target_no_11489_46089.txt",
        "Target_no_11489_47349.txt",
        "Target_no_11489_47757.txt",
        "Target_no_11489_51593.txt",
        "Target_no_11489_54754.txt",
        "Target_no_11489_58262.txt",
        "Target_no_11489_58958.txt",
        "Target_no_114_16827.txt",
        "Target_no_114_20503.txt",
        "Target_no_114_20724.txt",
        "Target_no_114_30254.txt",
        "Target_no_114_30321.txt",
        "Target_no_114_30572.txt",
        "Target_no_114_31116.txt",
        "Target_no_114_31281.txt",
        "Target_no_114_31443.txt",
        "Target_no_114_34747.txt",
        "Target_no_114_34920.txt",
        "Target_no_114_35040.txt",
        "Target_no_114_36623.txt",
        "Target_no_114_37136.txt",
        "Target_no_114_37208.txt",
        "Target_no_114_38298.txt",
        "Target_no_114_38554.txt",
        "Target_no_114_39088.txt",
        "Target_no_114_39543.txt",
        "Target_no_114_44049.txt",
        "Target_no_114_44562.txt",
        "Target_no_114_45180.txt",
        "Target_no_114_45867.txt",
        "Target_no_114_46087.txt",
        "Target_no_114_46829.txt",
        "Target_no_114_47546.txt",
        "Target_no_114_48395.txt",
        "Target_no_114_50461.txt",
        "Target_no_114_50765.txt",
        "Target_no_114_51669.txt",
        "Target_no_114_52381.txt",
        "Target_no_114_57231.txt",
        "Target_no_114_58246.txt",
        "Target_no_11534_30999.txt",
        "Target_no_11534_31021.txt",
        "Target_no_11534_31650.txt",
        "Target_no_11534_35011.txt",
        "Target_no_11534_37014.txt",
        "Target_no_11534_44598.txt",
        "Target_no_11534_48446.txt",
        "Target_no_11534_50614.txt",
        "Target_no_11534_50615.txt",
        "Target_no_11534_50618.txt",
        "Target_no_11534_51950.txt",
        "Target_no_11575_31616.txt",
        "Target_no_11575_38381.txt",
        "Target_no_11575_45757.txt",
        "Target_no_11575_46816.txt",
        "Target_no_11575_50701.txt",
        "Target_no_11575_54049.txt",
        "Target_no_11575_57890.txt",
        "Target_no_11575_61045.txt",
        "Target_no_11631_30194.txt",
        "Target_no_11631_31421.txt",
        "Target_no_11631_35538.txt",
        "Target_no_11631_55727.txt",
        "Target_no_11631_60676.txt",
        "Target_no_11631_61130.txt",
        "Target_no_11631_6444.txt",
        "Target_no_121_15623.txt",
        "Target_no_121_20096.txt",
        "Target_no_121_21591.txt",
        "Target_no_121_30184.txt",
        "Target_no_121_30380.txt",
        "Target_no_121_30745.txt",
        "Target_no_121_30861.txt",
        "Target_no_121_31619.txt",
        "Target_no_121_34869.txt",
        "Target_no_121_35284.txt",
        "Target_no_121_36158.txt",
        "Target_no_121_36915.txt",
        "Target_no_121_37434.txt",
        "Target_no_121_37541.txt",
        "Target_no_121_37677.txt",
        "Target_no_121_38837.txt",
        "Target_no_121_39530.txt",
        "Target_no_121_43793.txt",
        "Target_no_121_44158.txt",
        "Target_no_121_44164.txt",
        "Target_no_121_44579.txt",
        "Target_no_121_44593.txt",
        "Target_no_121_44692.txt",
        "Target_no_121_44914.txt",
        "Target_no_121_45878.txt",
        "Target_no_121_48050.txt",
        "Target_no_121_48706.txt",
        "Target_no_121_48944.txt",
        "Target_no_121_48961.txt",
        "Target_no_121_49880.txt",
        "Target_no_121_50055.txt",
        "Target_no_121_50512.txt",
        "Target_no_121_51064.txt",
        "Target_no_121_53537.txt",
        "Target_no_121_55749.txt",
        "Target_no_121_58879.txt",
        "Target_no_121_6876.txt",
        "Target_no_12209_44664.txt",
        "Target_no_12209_46451.txt",
        "Target_no_12209_50640.txt",
        "Target_no_12209_56102.txt",
        "Target_no_12209_58015.txt",
        "Target_no_12209_59089.txt",
        "Target_no_12209_59092.txt",
        "Target_no_12209_61340.txt",
        "Target_no_12252_31458.txt",
        "Target_no_12252_38481.txt",
        "Target_no_12252_42104.txt",
        "Target_no_12252_49475.txt",
        "Target_no_12252_49933.txt",
        "Target_no_12252_50558.txt",
        "Target_no_126_37249.txt",
        "Target_no_126_39117.txt",
        "Target_no_126_48227.txt",
        "Target_no_126_48229.txt",
        "Target_no_126_53969.txt",
        "Target_no_126_54262.txt",
        "Target_no_126_55353.txt",
        "Target_no_126_61212.txt",
        "Target_no_12952_20250.txt",
        "Target_no_12952_20506.txt",
        "Target_no_12952_20626.txt",
        "Target_no_12952_20689.txt",
        "Target_no_12952_21358.txt",
        "Target_no_12952_30158.txt",
        "Target_no_12952_30349.txt",
        "Target_no_12952_34843.txt",
        "Target_no_12952_35220.txt",
        "Target_no_12952_36318.txt",
        "Target_no_12952_36800.txt",
        "Target_no_12952_39349.txt",
        "Target_no_12952_45385.txt",
        "Target_no_12952_45676.txt",
        "Target_no_12952_46700.txt",
        "Target_no_12952_46719.txt",
        "Target_no_12952_48517.txt",
        "Target_no_12952_52690.txt",
        "Target_no_12952_53483.txt",
        "Target_no_12952_55435.txt",
        "Target_no_12952_56102.txt",
        "Target_no_12952_58015.txt",
        "Target_no_12952_58375.txt",
        "Target_no_12952_58700.txt",
        "Target_no_13001_16724.txt",
        "Target_no_13001_2710.txt",
        "Target_no_13001_30128.txt",
        "Target_no_13001_31250.txt",
        "Target_no_13001_35519.txt",
        "Target_no_13001_36905.txt",
        "Target_no_13001_46788.txt",
        "Target_no_13001_4737.txt",
        "Target_no_13001_55736.txt",
        "Target_no_13001_57473.txt",
        "Target_no_130_16146.txt",
        "Target_no_130_30313.txt",
        "Target_no_130_30539.txt",
        "Target_no_130_30653.txt",
        "Target_no_130_30988.txt",
        "Target_no_130_31179.txt",
        "Target_no_130_31339.txt",
        "Target_no_130_3525.txt",
        "Target_no_130_36886.txt",
        "Target_no_130_37020.txt",
        "Target_no_130_38128.txt",
        "Target_no_130_43952.txt",
        "Target_no_130_45210.txt",
        "Target_no_130_45899.txt",
        "Target_no_130_47889.txt",
        "Target_no_130_51800.txt",
        "Target_no_130_53158.txt",
        "Target_no_130_5330.txt",
        "Target_no_130_55648.txt",
        "Target_no_130_56096.txt",
        "Target_no_130_57287.txt",
        "Target_no_130_57399.txt",
        "Target_no_130_57606.txt",
        "Target_no_130_58470.txt",
        "Target_no_130_60894.txt",
        "Target_no_15_15010.txt",
        "Target_no_15_20304.txt",
        "Target_no_15_20542.txt",
        "Target_no_15_20788.txt",
        "Target_no_15_21358.txt",
        "Target_no_15_30349.txt",
        "Target_no_15_31103.txt",
        "Target_no_15_31764.txt",
        "Target_no_15_31822.txt",
        "Target_no_15_35220.txt",
        "Target_no_15_37964.txt",
        "Target_no_15_42164.txt",
        "Target_no_15_44072.txt",
        "Target_no_15_45522.txt",
        "Target_no_15_46719.txt",
        "Target_no_15_46744.txt",
        "Target_no_15_48416.txt",
        "Target_no_15_49740.txt",
        "Target_no_15_50034.txt",
        "Target_no_15_50640.txt",
        "Target_no_15_53483.txt",
        "Target_no_15_53506.txt",
        "Target_no_15_53814.txt",
        "Target_no_15_55435.txt",
        "Target_no_15_56074.txt",
        "Target_no_15_56107.txt",
        "Target_no_15_57090.txt",
        "Target_no_15_57705.txt",
        "Target_no_15_57785.txt",
        "Target_no_15_57861.txt",
        "Target_no_15_58583.txt",
        "Target_no_15_59089.txt",
        "Target_no_15_59092.txt",
        "Target_no_165_20600.txt",
        "Target_no_165_30764.txt",
        "Target_no_165_34876.txt",
        "Target_no_165_37427.txt",
        "Target_no_165_38266.txt",
        "Target_no_165_45283.txt",
        "Target_no_165_46990.txt",
        "Target_no_165_47123.txt",
        "Target_no_165_47900.txt",
        "Target_no_165_48050.txt",
        "Target_no_165_48062.txt",
        "Target_no_165_48076.txt",
        "Target_no_165_50718.txt",
        "Target_no_165_50972.txt",
        "Target_no_165_51372.txt",
        "Target_no_165_51830.txt",
        "Target_no_165_53049.txt",
        "Target_no_165_53505.txt",
        "Target_no_165_54021.txt",
        "Target_no_165_54565.txt",
        "Target_no_165_54570.txt",
        "Target_no_165_57666.txt",
        "Target_no_165_57890.txt",
        "Target_no_165_58072.txt",
        "Target_no_165_58634.txt",
        "Target_no_165_58643.txt",
        "Target_no_165_59080.txt",
        "Target_no_165_59202.txt",
        "Target_no_165_59346.txt",
        "Target_no_165_60379.txt",
        "Target_no_165_60385.txt",
        "Target_no_165_61098.txt",
        "Target_no_17045_31154.txt",
        "Target_no_17045_35648.txt",
        "Target_no_17045_48914.txt",
        "Target_no_17045_51975.txt",
        "Target_no_17045_52435.txt",
        "Target_no_17045_55790.txt",
        "Target_no_17045_58958.txt",
        "Target_no_17045_61146.txt",
        "Target_no_19905_20191.txt",
        "Target_no_19905_20562.txt",
        "Target_no_19905_20606.txt",
        "Target_no_19905_20607.txt",
        "Target_no_19905_30003.txt",
        "Target_no_19905_30402.txt",
        "Target_no_19905_30726.txt",
        "Target_no_19905_30924.txt",
        "Target_no_19905_31457.txt",
        "Target_no_19905_31460.txt",
        "Target_no_19905_31600.txt",
        "Target_no_19905_31620.txt",
        "Target_no_19905_31640.txt",
        "Target_no_19905_31782.txt",
        "Target_no_19905_31797.txt",
        "Target_no_19905_31829.txt",
        "Target_no_19905_31856.txt",
        "Target_no_19905_34690.txt",
        "Target_no_19905_34691.txt",
        "Target_no_19905_34694.txt",
        "Target_no_19905_35712.txt",
        "Target_no_19905_37497.txt",
        "Target_no_19905_37539.txt",
        "Target_no_19905_48070.txt",
        "Target_no_19905_48219.txt",
        "Target_no_19905_48824.txt",
        "Target_no_19905_48902.txt",
        "Target_no_19905_52867.txt",
        "Target_no_19905_53534.txt",
        "Target_no_19905_54570.txt",
        "Target_no_19905_58727.txt",
        "Target_no_19905_60159.txt",
        "Target_no_19905_60661.txt",
        "Target_no_259_30037.txt",
        "Target_no_259_34734.txt",
        "Target_no_259_35365.txt",
        "Target_no_259_35787.txt",
        "Target_no_259_36793.txt",
        "Target_no_259_37845.txt",
        "Target_no_259_38328.txt",
        "Target_no_259_39374.txt",
        "Target_no_259_43936.txt",
        "Target_no_259_44059.txt",
        "Target_no_259_44405.txt",
        "Target_no_259_44424.txt",
        "Target_no_259_45166.txt",
        "Target_no_259_46158.txt",
        "Target_no_259_46425.txt",
        "Target_no_259_47187.txt",
        "Target_no_259_48435.txt",
        "Target_no_259_48817.txt",
        "Target_no_259_48849.txt",
        "Target_no_259_48957.txt",
        "Target_no_259_48965.txt",
        "Target_no_259_49519.txt",
        "Target_no_259_49921.txt",
        "Target_no_259_50392.txt",
        "Target_no_259_50972.txt",
        "Target_no_259_52405.txt",
        "Target_no_259_52442.txt",
        "Target_no_259_53169.txt",
        "Target_no_259_54640.txt",
        "Target_no_259_57207.txt",
        "Target_no_259_58537.txt",
        "Target_no_259_58920.txt",
        "Target_no_259_60358.txt",
        "Target_no_25_16077.txt",
        "Target_no_25_34943.txt",
        "Target_no_25_36391.txt",
        "Target_no_25_36574.txt",
        "Target_no_25_52906.txt",
        "Target_no_25_54642.txt",
        "Target_no_43_16385.txt",
        "Target_no_43_20333.txt",
        "Target_no_43_20774.txt",
        "Target_no_43_35825.txt",
        "Target_no_43_37637.txt",
        "Target_no_43_47741.txt",
        "Target_no_43_48567.txt",
        "Target_no_43_48841.txt",
        "Target_no_43_49589.txt",
        "Target_no_43_50564.txt",
        "Target_no_43_51201.txt",
        "Target_no_43_53150.txt",
        "Target_no_43_55747.txt",
        "Target_no_43_58276.txt",
        "Target_no_51_20813.txt",
        "Target_no_51_2331.txt",
        "Target_no_51_2391.txt",
        "Target_no_51_30374.txt",
        "Target_no_51_30512.txt",
        "Target_no_51_30796.txt",
        "Target_no_51_31058.txt",
        "Target_no_51_31450.txt",
        "Target_no_51_36609.txt",
        "Target_no_51_36636.txt",
        "Target_no_51_38266.txt",
        "Target_no_51_43800.txt",
        "Target_no_51_47621.txt",
        "Target_no_51_49433.txt",
        "Target_no_51_49550.txt",
        "Target_no_51_49865.txt",
        "Target_no_51_51005.txt",
        "Target_no_51_52373.txt",
        "Target_no_51_53878.txt",
        "Target_no_51_54554.txt",
        "Target_no_51_57508.txt",
        "Target_no_51_57618.txt",
        "Target_no_51_59484.txt",
        "Target_no_51_6563.txt",
        "Target_no_61_16845.txt",
        "Target_no_61_37106.txt",
        "Target_no_61_37272.txt",
        "Target_no_61_47947.txt",
        "Target_no_61_50754.txt",
        "Target_no_61_55800.txt",
        "Target_no_61_57380.txt",
        "Target_no_65_35328.txt",
        "Target_no_65_36137.txt",
        "Target_no_65_38826.txt",
        "Target_no_65_44596.txt",
        "Target_no_65_51275.txt",
        "Target_no_65_55671.txt",
        "Target_no_72_20155.txt",
        "Target_no_72_30616.txt",
        "Target_no_72_30653.txt",
        "Target_no_72_30988.txt",
        "Target_no_72_31339.txt",
        "Target_no_72_31527.txt",
        "Target_no_72_31631.txt",
        "Target_no_72_36047.txt",
        "Target_no_72_36568.txt",
        "Target_no_72_36636.txt",
        "Target_no_72_37020.txt",
        "Target_no_72_38128.txt",
        "Target_no_72_43952.txt",
        "Target_no_72_45210.txt",
        "Target_no_72_47889.txt",
        "Target_no_72_51064.txt",
        "Target_no_72_51800.txt",
        "Target_no_72_5330.txt",
        "Target_no_72_57508.txt",
        "Target_no_72_58470.txt",
        "Target_no_72_60894.txt",
        "Target_no_87_20730.txt",
        "Target_no_87_20753.txt",
        "Target_no_87_30097.txt",
        "Target_no_87_30663.txt",
        "Target_no_87_35365.txt",
        "Target_no_87_35375.txt",
        "Target_no_87_35918.txt",
        "Target_no_87_36010.txt",
        "Target_no_87_36780.txt",
        "Target_no_87_36793.txt",
        "Target_no_87_37845.txt",
        "Target_no_87_39110.txt",
        "Target_no_87_39546.txt",
        "Target_no_87_43936.txt",
        "Target_no_87_44405.txt",
        "Target_no_87_44424.txt",
        "Target_no_87_44683.txt",
        "Target_no_87_45166.txt",
        "Target_no_87_46648.txt",
        "Target_no_87_46986.txt",
        "Target_no_87_47070.txt",
        "Target_no_87_47187.txt",
        "Target_no_87_48435.txt",
        "Target_no_87_48817.txt",
        "Target_no_87_49138.txt",
        "Target_no_87_49460.txt",
        "Target_no_87_49519.txt",
        "Target_no_87_49532.txt",
        "Target_no_87_49829.txt",
        "Target_no_87_50608.txt",
        "Target_no_87_54640.txt",
        "Target_no_87_61082.txt",
        "Target_no_90_15727.txt",
        "Target_no_90_20155.txt",
        "Target_no_90_30407.txt",
        "Target_no_90_3525.txt",
        "Target_no_93_12377.txt",
        "Target_no_93_1308.txt",
        "Target_no_93_34759.txt",
        "Target_no_93_35447.txt",
        "Target_no_93_42214.txt",
        "Target_no_93_45877.txt",
        "Target_no_93_47082.txt",
        "Target_no_93_57082.txt",
    };
    FILE *fcmd = fopen("chembl_II_sets.bat", "wt");
    // commands for prepare Python test:
    fprintf(fcmd, "DEL %s\n", "chembl_II_sets.P.res.csv");   //clear before append results
    fprintf(fcmd, "SET PATH=%PATH%;C:/LIB\n");
    fprintf(fcmd, "SET PYTHONPATH=C:/Projects/RDKit/RDKit_2013_09_1\n");
    fprintf(fcmd, "ECHO P test;P Nmols;P status;P time,sec;P nAtoms;P nBonds;P MCS >%s\n", "chembl_II_sets.P.res.csv");

    p.Timeout = 60;
    p.Verbose = false;

    for(int i=0; i<sizeof(test)/sizeof(test[0]); i++) {
        std::string smiName = std::string("chembl_II_sets/") + test[i];
        std::string testsmi = testChEMBL_Txt(smiName.c_str(), th, "chembl_II_sets.C.res.csv");
        fprintf(fcmd, "fmcs_bench.py --timeout %u %s >>%s\n", p.Timeout, testsmi.c_str(), "chembl_II_sets.P.res.csv");   // command for the same Python test
    }
}

void testChEMBL_TxtSLOW_chembl_II_sets(double th=1.0) {
    FILE* fres = fopen("chembl_II_sets.SLOW.C++.res.csv", "wt");
    fprintf(fres,"test;threshold;t,sec;atoms;bonds;mcs\n");
    fclose(fres);
    const char* test[] = {
        "Target_no_10980_51302.txt",
        "Target_no_114_31443.txt",
        "Target_no_10260_54285.txt",
        "Target_no_11489_37339.txt",
        "Target_no_10980_52937.txt",
        "Target_no_114_48395.txt",
        "Target_no_10280_45765.txt",
        "Target_no_10188_50681.txt",
        "Target_no_114_46087.txt",
        "Target_no_10188_58358.txt",
        "Target_no_11140_54121.txt",
        "Target_no_114_37208.txt",
        "Target_no_10193_53483.txt",
        "Target_no_10980_31623.txt",
        "Target_no_10280_60672.txt",
        "Target_no_10188_45279.txt",
        "Target_no_11140_49860.txt",
        "Target_no_10280_53081.txt",
        "Target_no_11489_46089.txt",
        "Target_no_10188_54039.txt",
        "Target_no_114_34747.txt",
        "Target_no_11365_30331.txt",
        "Target_no_10188_48029.txt",
        "Target_no_114_16827.txt",
        "Target_no_100_30745.txt",
        "Target_no_108_49591.txt",
        "Target_no_107_49591.txt",
        "Target_no_10980_35957.txt",
        "Target_no_11489_37318.txt",
        "Target_no_10188_58265.txt",
        "Target_no_107_58879.txt",
        "Target_no_11140_37272.txt",
        "Target_no_10188_39262.txt",
        "Target_no_10280_47928.txt",
        "Target_no_107_49345.txt",
        "Target_no_10188_20244.txt",
        "Target_no_10980_20853.txt",
        "Target_no_11365_58566.txt",
        "Target_no_114_44562.txt",
        "Target_no_100_16613.txt",
        "Target_no_11489_39561.txt",
        "Target_no_10434_35926.txt",
        "Target_no_10980_31508.txt",
        "Target_no_11489_34818.txt",
        "Target_no_10188_58866.txt",
        "Target_no_11489_58958.txt",
        "Target_no_11140_39472.txt",
        "Target_no_10980_5739.txt",
        "Target_no_11489_54754.txt",
        "Target_no_10260_38526.txt",
        "Target_no_11489_58262.txt",
        "Target_no_107_34921.txt",
        "Target_no_107_30866.txt",
        "Target_no_108_30866.txt",
        "Target_no_10980_36641.txt",
        "Target_no_10980_30840.txt",
        "Target_no_107_31249.txt",
        "Target_no_11140_53395.txt",
        "Target_no_10193_46700.txt",
        "Target_no_10980_30994.txt",
        "Target_no_11140_37038.txt",
    };
    for(int i=0; i<sizeof(test)/sizeof(test[0]); i++)
        testChEMBL_Txt((std::string("chembl_II_sets/") + test[i]).c_str(), th,"chembl_II_sets.SLOW.C++.res.csv");
}

void testChEMBLdat(const char* test, double th=1.0) {
    std::cout << "\ntestChEMBLdat() "<<test<<"\n";
    std::vector<ROMOL_SPTR> mols;
    char smiles[4096];
    unsigned n=0;
    FILE* f = fopen(test, "rt");
    if(!f) {
        perror("fopen testChEMBLdat()");
        return;
    }
    FILE* fs = fopen((std::string(test)+".smi").c_str(), "wt");
    std::cout<<"Loading SMILES ... \n";
    t0 = nanoClock();
    while(fgets(smiles, sizeof(smiles), f)) {
        std::cout<<"\rLine: "<< ++n <<" ";
        if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) { // commented to skip
            mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnlyChEMBL(smiles))));
            fputs(getSmilesOnlyChEMBL(smiles).c_str(), fs);
        }
    }
    fclose(f);
    fclose(fs);
    printTime();
    std::cout<<"FIND MCS in "<<mols.size()<<" molecules.\n\n";
    t0 = nanoClock();
    p.Threshold = th;
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testChEMBLdatALL(double th=1.0) {
    const char* test[] = {
        "cmp_list_ChEMBL_100126_actives.dat",
        "cmp_list_ChEMBL_100579_actives.dat",
        "cmp_list_ChEMBL_100_actives.dat",
        "cmp_list_ChEMBL_10188_actives.dat",
        "cmp_list_ChEMBL_10193_actives.dat",
        "cmp_list_ChEMBL_10198_actives.dat",
        "cmp_list_ChEMBL_10260_actives.dat",
        "cmp_list_ChEMBL_10280_actives.dat",
        "cmp_list_ChEMBL_10378_actives.dat",
        "cmp_list_ChEMBL_10417_actives.dat",
        "cmp_list_ChEMBL_10434_actives.dat",
        "cmp_list_ChEMBL_10475_actives.dat",
        "cmp_list_ChEMBL_10498_actives.dat",
        "cmp_list_ChEMBL_104_actives.dat",
        "cmp_list_ChEMBL_10579_actives.dat",
        "cmp_list_ChEMBL_105_actives.dat",
        "cmp_list_ChEMBL_10752_actives.dat",
        "cmp_list_ChEMBL_10773_actives.dat",
        "cmp_list_ChEMBL_107_actives.dat",
        "cmp_list_ChEMBL_108_actives.dat",
        "cmp_list_ChEMBL_10927_actives.dat",
        "cmp_list_ChEMBL_10980_actives.dat",
        "cmp_list_ChEMBL_11085_actives.dat",
        "cmp_list_ChEMBL_11140_actives.dat",
        "cmp_list_ChEMBL_11225_actives.dat",
        "cmp_list_ChEMBL_11279_actives.dat",
        "cmp_list_ChEMBL_11336_actives.dat",
        "cmp_list_ChEMBL_11359_actives.dat",
        "cmp_list_ChEMBL_11365_actives.dat",
        "cmp_list_ChEMBL_11442_actives.dat",
        "cmp_list_ChEMBL_11488_actives.dat",
        "cmp_list_ChEMBL_11489_actives.dat",
        "cmp_list_ChEMBL_114_actives.dat",
        "cmp_list_ChEMBL_11534_actives.dat",
        "cmp_list_ChEMBL_11536_actives.dat",
        "cmp_list_ChEMBL_11575_actives.dat",
        "cmp_list_ChEMBL_11631_actives.dat",
        "cmp_list_ChEMBL_11682_actives.dat",
        "cmp_list_ChEMBL_116_actives.dat",
        "cmp_list_ChEMBL_121_actives.dat",
        "cmp_list_ChEMBL_12209_actives.dat",
        "cmp_list_ChEMBL_12252_actives.dat",
        "cmp_list_ChEMBL_12261_actives.dat",
        "cmp_list_ChEMBL_12670_actives.dat",
        "cmp_list_ChEMBL_12679_actives.dat",
        "cmp_list_ChEMBL_126_actives.dat",
        "cmp_list_ChEMBL_12840_actives.dat",
        "cmp_list_ChEMBL_12911_actives.dat",
        "cmp_list_ChEMBL_12952_actives.dat",
        "cmp_list_ChEMBL_12968_actives.dat",
        "cmp_list_ChEMBL_13001_actives.dat",
        "cmp_list_ChEMBL_130_actives.dat",
        "cmp_list_ChEMBL_134_actives.dat",
        "cmp_list_ChEMBL_15_actives.dat",
        "cmp_list_ChEMBL_165_actives.dat",
        "cmp_list_ChEMBL_17045_actives.dat",
        "cmp_list_ChEMBL_18061_actives.dat",
        "cmp_list_ChEMBL_19905_actives.dat",
        "cmp_list_ChEMBL_20014_actives.dat",
        "cmp_list_ChEMBL_20174_actives.dat",
        "cmp_list_ChEMBL_219_actives.dat",
        "cmp_list_ChEMBL_234_actives.dat",
        "cmp_list_ChEMBL_237_actives.dat",
        "cmp_list_ChEMBL_259_actives.dat",
        "cmp_list_ChEMBL_25_actives.dat",
        "cmp_list_ChEMBL_276_actives.dat",
        "cmp_list_ChEMBL_28_actives.dat",
        "cmp_list_ChEMBL_36_actives.dat",
        "cmp_list_ChEMBL_43_actives.dat",
        "cmp_list_ChEMBL_51_actives.dat",
        "cmp_list_ChEMBL_52_actives.dat",
        "cmp_list_ChEMBL_61_actives.dat",
        "cmp_list_ChEMBL_65_actives.dat",
        "cmp_list_ChEMBL_72_actives.dat",
        "cmp_list_ChEMBL_87_actives.dat",
        "cmp_list_ChEMBL_8_actives.dat",
        "cmp_list_ChEMBL_90_actives.dat",
        "cmp_list_ChEMBL_93_actives.dat",
        "cmp_list_ChEMBL_zinc_decoys.dat",
    };
    for(int i=0; i<sizeof(test)/sizeof(test[0]); i++)
        testChEMBLdat((std::string("benchmarking_platform-master/compounds/ChEMBL/") + test[i]).c_str(), th);
}


void testTarget_no_10188_30149() {
    std::cout << "\ntestTarget_no_10188_30149()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] = {
        //Target_no_10188_30149.txt // VERY SLOWER than Python
        "CN(C)CCNC(=O)c1ccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL399167",
        "O=C(O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL197613",
        "c1ccc(Nc2ccc3c(c2)[nH]nc3-c2ccccc2)cc1 CHEMBL383177",          /// == QUERY
        "NC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL199136",
        "Clc1ccccc1Nc1ccc2c(c1)n[nH]c2-c1ccccc1 CHEMBL440566",
        "O=C(NCCCN1CCOCC1)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198687",
        "O=C(O)c1ccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL197698",
        "O=C(NC1CCNCC1)c1cccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL194806",
        "COc1ccccc1Nc1ccc2c(c1)[nH]nc2-c1ccccc1 CHEMBL254443",
        "CN(C)CCNC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198821",
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
    t0 = nanoClock();
#ifdef _DEBUG   // check memory leaks
    _CrtMemState _ms;
    _CrtMemCheckpoint(&_ms);
#endif
    {
        MCSResult res = findMCS(mols, &p);
        std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    }
#ifdef _DEBUG   // check memory leaks
    _CrtMemDumpAllObjectsSince(&_ms);
#endif
    printTime();
}


void testTarget_no_10188_49064() {
    std::cout << "\ntestTarget_no_10188_49064()\n";
    std::vector<ROMOL_SPTR> mols;
    const char* smi[] = {
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(O)c3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(F)cc3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(OC4OC(CO)C(O)C(O)C4O)cc3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCl)cc3)nc21",
        "Cn1c2nc(Nc3ccc(NC(=O)CCN)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(COCC(O)CO)c3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(O)cc3)nc21",
        "CC(=O)Nc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
        "Cn1c2nc(Nc3ccc(N)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCNC(=O)OC(C)(C)C)cc3)nc21",
        "Cc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(CO)c3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NCC(O)CO)cc3)nc21",
        "CCc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
        "Cn1c2nc(Nc3cccc(N)c3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
        "CC(=O)Nc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(CCO)cc3)nc21",
        "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(I)cc3)nc21",
        "CN1CCN(C(=O)c2ccc(Nc3ncc4cc(-c5c(Cl)cccc5Cl)c(=O)n(C)c4n3)cc2)CC1",
    };
    for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

void testCmndLineSMILES(int argc, const char* argv[]) {
    std::vector<ROMOL_SPTR> mols;
    for(int i=1; i<argc; i++)
        mols.push_back(ROMOL_SPTR(SmilesToMol(argv[i])));
    MCSResult res = findMCS(mols);
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}


double testFileSDF(const char* test) {
    std::cout << "\ntestFileSDF(): " << test << "\n";
    std::vector<ROMOL_SPTR> mols;
    std::string fn(test);
    RDKit::SDMolSupplier suppl(fn);
    while(!suppl.atEnd()) {
        ROMol *m=suppl.next();
        if(m)
            mols.push_back(ROMOL_SPTR(m));
    }
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    double t = (nanoClock() - t0) / 1000000.;
    std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
    return t;
}


void testFileSDF_RandomSet_SMI(const char* path="benchmark", const char* test="chembl13-10000-random-pairs.sdf") {
    p.Timeout = 60;
    p.Verbose = false;
    unsigned long long To = nanoClock();

    std::cout << "\ntestFileSDF_RandomSet(): " << path<<"/"<<test << "\n";
    std::vector<ROMOL_SPTR> mols;
    FILE* fcsv = fopen((std::string(path)+"_"+test+".SMI.C.csv").c_str(), "wt");
    fprintf(fcsv, "test;Nmols;status;t,sec;nAtoms;nBonds;C++ MCS; E status;E t,sec;E nAtoms;E nBonds;E C++ MCS\n");

    bool fileExist = true;
    unsigned n = 1;
    for( ; fileExist ; n++) {
        char smiName[256];
        sprintf(smiName,"%s/smiles/%s.%u.smi", path, test, n);
        FILE* fsmi = fopen(smiName, "rt");
        if(!fsmi) {
            fileExist = false;
            break;
        }
        char smiles[4096];
        while(fgets(smiles, sizeof(smiles), fsmi)) {
            mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(smiles))));
        }
        fclose(fsmi);
        if(mols.size()>1) {
            t0 = nanoClock();
            MCSResult res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            fprintf(fcsv, "%u;%u;%s;%.3f;%u;%u;%s;", n, mols.size(), res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
        }
        if(mols.size()>1) {
            p.BondTyper = MCSBondCompareOrderExact;
            t0 = nanoClock();
            MCSResult res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            fprintf(fcsv, "%u;%u;%s;%.3f;%u;%u;%s\n", n, mols.size(), res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
            p.BondTyper = MCSBondCompareOrder;
        }
        mols.clear();
    }
    fclose(fcsv);
}

void testFileSDF_RandomSet(const char* test="chembl13-10000-random-pairs.sdf", const char* path="benchmark") {
    p.Timeout = 60;
    p.Verbose = false;

    unsigned long long To = nanoClock();

    std::cout << "\ntestFileSDF_RandomSet(): " << path<<"/"<<test << "\n";
    std::vector<ROMOL_SPTR> mols;
    std::vector<ROMOL_SPTR> all_mols;
    FILE* fall = fopen((std::string(path)+"_ALL.P.bat").c_str(), "at");
    fprintf(fall, "CALL %s\n", (std::string(path)+"_"+test+".P.bat").c_str());
    fclose(fall);
    FILE* fcmd = fopen((std::string(path)+"_"+test+".P.bat").c_str(), "wt");
    FILE* fcsv = fopen((std::string(path)+"_"+test+".C.csv").c_str(), "wt");
    fprintf(fcsv, "test;Nmols;status;t,sec;nAtoms;nBonds;C++ MCS; E status;E t,sec;E nAtoms;E nBonds;E C++ MCS\n");
    std::string fn(std::string(path)+"/"+test);
    RDKit::MolSupplier* suppl = 0;
    try {
        if('f' == test[strlen(test)-1])   // sdf file
            suppl = new RDKit::SDMolSupplier(fn);
        else if('i' == test[strlen(test)-1])   // smi file
            suppl = new RDKit::SmilesMolSupplier(fn);
    } catch(...) {
        std::cout << "ERROR: RDKit could not load input file" << "\n";
        return;
    }
    if(!suppl) {
        std::cout << "ERROR: unsupported input file format" << "\n";
        return;
    }
    // commands for prepare Python test:
    fprintf(fcmd, "DEL %s\n", (std::string(path)+"_"+test+".P.csv").c_str());   //clear before append results
    fprintf(fcmd, "SET PATH=%PATH%;C:/LIB\n");
    fprintf(fcmd, "SET PYTHONPATH=C:/Projects/RDKit/RDKit_2013_09_1\n");
    fprintf(fcmd, "ECHO P test;P Nmols;P status;P time,sec;P nAtoms;P nBonds;P MCS >%s\n", (std::string(path)+"_"+test+".P.csv").c_str());

    unsigned n = 1;
    std::cout << "\n****** Load all molecules from SFD and PAIR SET test *********\n\n";
    for( ; !suppl->atEnd() ; n++) {
        char smiName[256];
        sprintf(smiName,"%s/smiles/%s.%u.smi", path, test, n);
        FILE* fsmi = fopen(smiName, "wt");
        if(!fsmi) {
            std::cout << "ERROR: could not create SMI file " << smiName <<"\n";
            return;
        }
        ROMol *m=0;
        for(int i=0; i<2 && !suppl->atEnd(); i++) {  // load sequential pair
            m = suppl->next();
            if(m) {
                mols.push_back(ROMOL_SPTR(m));
                all_mols.push_back(mols.back());
                fprintf(fsmi,"%s Mol%u\n", MolToSmiles(*m).c_str(), n+i);
            }
        }
        fclose(fsmi);
        if(mols.size()<2) {
            unlink(smiName);
            n--;
        } else {
            fprintf(fcmd, "fmcs_bench.py --id %u --timeout %u --threshold %.2f %s >>%s\n", n, p.Timeout, p.Threshold, smiName, (std::string(path)+"_"+test+".P.csv").c_str());   // command for the same Python test

            if(mols.size()>1) {
                t0 = nanoClock();
                MCSResult res = findMCS(mols, &p);
                double t = (nanoClock() - t0) / 1000000.;
                if(t < 0.00001)
                    t = 0.00001; // avoid division by zero
                printTime();
                std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
                fprintf(fcsv, "%u;%u;%s;%.5f;%u;%u;%s;", n, mols.size(), res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
            }
            if(mols.size()>1) {
                p.BondTyper = MCSBondCompareOrderExact;
                t0 = nanoClock();
                MCSResult res = findMCS(mols, &p);
                double t = (nanoClock() - t0) / 1000000.;
                if(t < 0.00001)
                    t = 0.00001; // avoid division by zero
                printTime();
                std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
                fprintf(fcsv, "%s;%.5f;%u;%u;%s\n", res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
                p.BondTyper = MCSBondCompareOrder;
            }
        }
        mols.clear();
    }
    delete(suppl);

    std::cout << "\n****** RANDOM SET test *********\n\n";
    const unsigned N_RandomTests = all_mols.size() * 7;
    srand(1);   // make stable pseudorandom sequence
    for(unsigned jn; n <= N_RandomTests; jn++, n++) {
        char smiName[256];
        sprintf(smiName,"%s/smilesRAND/%s.%u.smi", path, test, n);
        FILE* fsmi = fopen(smiName, "wt");
        fprintf(fcmd, "fmcs_bench.py --id %u --timeout %u --threshold %.2f %s >>%s\n", n, p.Timeout, p.Threshold, smiName, (std::string(path)+"_"+test+".P.csv").c_str());   // command for the same Python test
        ROMol *m=0;
        unsigned iN = 3+rand()%32;
        for(int i=0; i < iN; i++) {  // load random set
            mols.push_back(all_mols[rand()%(all_mols.size()-1)]);
            fprintf(fsmi,"%s Mol%u\n", MolToSmiles(*mols.back()).c_str(), n+i);
        }
        fclose(fsmi);
        if(mols.size()>1) {
            t0 = nanoClock();
            MCSResult res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            if(t < 0.00001)
                t = 0.00001; // avoid division by zero
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            fprintf(fcsv, "%u;%u;%s;%.5f;%u;%u;%s;", n, mols.size(), res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
        }
        if(mols.size()>1) {
            p.BondTyper = MCSBondCompareOrderExact;
            t0 = nanoClock();
            MCSResult res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            if(t < 0.00001)
                t = 0.00001; // avoid division by zero
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            fprintf(fcsv, "%s;%.5f;%u;%u;%s\n", res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
            p.BondTyper = MCSBondCompareOrder;
        }
        mols.clear();
    }
    fclose(fcsv);
    fclose(fcmd);

    std::cout << "\n****** BIG MCS RANDOM SET test *********\n\n";

    unsigned maxMol = 0;
    for(unsigned i=0; i < all_mols.size(); i++) {
        if(maxMol < all_mols[i]->getNumBonds())
            maxMol = all_mols[i]->getNumBonds();
    }
    const unsigned N_BigRandomTests = all_mols.size() > 2000 ? all_mols.size()/2 : all_mols.size()*2;
    const unsigned N_BigRandomTestsAttempts = all_mols.size()*130;
    const unsigned SizeOfBigMCS_ForBigRandomTests = maxMol < 24 ? maxMol*2/3 : 21;

    fall = fopen((std::string(path)+"_ALL_BIG.P.bat").c_str(), "at");
    fprintf(fall, "CALL %s\n", (std::string(path)+"_"+test+".BIG_MCS.P.bat").c_str());
    fclose(fall);
    fcmd = fopen((std::string(path)+"_"+test+".BIG_MCS.P.bat").c_str(), "wt");
    fcsv = fopen((std::string(path)+"_"+test+".BIG_MCS.C.csv").c_str(), "wt");
    fprintf(fcsv, "test;Nmols;status;t,sec;nAtoms;nBonds;C++ MCS; E status;E t,sec;E nAtoms;E nBonds;E C++ MCS\n");

    const unsigned n1 = n;
    fprintf(fcmd, "DEL %s\n", (std::string(path)+"_"+test+".BIG_MCS.P.csv").c_str());   //clear before append results
    fprintf(fcmd, "SET PATH=%PATH%;C:/LIB\n");
    fprintf(fcmd, "SET PYTHONPATH=C:/Projects/RDKit/RDKit_2013_09_1\n");
    fprintf(fcmd, "ECHO P test;P Nmols;P status;P time,sec;P nAtoms;P nBonds;P MCS >%s\n", (std::string(path)+"_"+test+".BIG_MCS.P.csv").c_str());
    for(size_t jn=0; jn < N_BigRandomTestsAttempts && n-n1 <= N_BigRandomTests; jn++) {
        char smiName[512];
        sprintf(smiName,"%s/smilesBIG/%s.%u.smi", path, test, n);
        ROMol *m=0;
        unsigned iN = 2+rand()%24;
        mols.clear();
        for(size_t i=0; i < iN; i++)    // load random set
            for(size_t ij=0; ij < all_mols.size()/2; ij++) {
                unsigned mi = rand()%(all_mols.size()-1);
                if(all_mols[mi]->getNumBonds() < SizeOfBigMCS_ForBigRandomTests)
                    continue;
                mols.push_back(all_mols[mi]);
                break;
            }
        MCSResult res;
        if(mols.size()>1) {
            t0 = nanoClock();
            res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            if(t < 0.00001)
                t = 0.00001; // avoid division by zero
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            if(res.NumBonds >= SizeOfBigMCS_ForBigRandomTests && res.isCompleted()) {
                fprintf(fcsv, "%u;%u;%s;%.5f;%u;%u;%s;", n, mols.size(), res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
                fprintf(fcmd, "fmcs_bench.py --id %u --timeout %u --threshold %.2f %s >>%s\n", n, p.Timeout, p.Threshold, smiName, (std::string(path)+"_"+test+".BIG_MCS.P.csv").c_str());   // command for the same Python test
                n++;
            }
        }
        if(res.NumBonds >= SizeOfBigMCS_ForBigRandomTests && res.isCompleted()) {
            FILE* fsmi = fopen(smiName, "wt");
            for(int i=0; i < mols.size(); i++) {  // load random set
                fprintf(fsmi,"%s Mol%u\n", MolToSmiles(*mols[i]).c_str(), n+i);
            }
            fclose(fsmi);
        }
        if(res.NumBonds >= SizeOfBigMCS_ForBigRandomTests && res.isCompleted()) {
            p.BondTyper = MCSBondCompareOrderExact;
            t0 = nanoClock();
            res = findMCS(mols, &p);
            double t = (nanoClock() - t0) / 1000000.;
            if(t < 0.00001)
                t = 0.00001; // avoid division by zero
            printTime();
            std::cout << n <<" MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
            fprintf(fcsv, "%s;%.5f;%u;%u;%s\n", res.isCompleted()?" ":"TIMEOUT", t, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
            p.BondTyper = MCSBondCompareOrder;
        }
        mols.clear();
    }
    fclose(fcsv);
    fclose(fcmd);

    t0 = T0;
    std::cout << n << " tests processed\n";
    printTime();
}

void testFileSMILES(const char* test) {
    std::vector<ROMOL_SPTR> mols;
    char smiles[4096];
    unsigned n=0;
    FILE* f = fopen(test, "rt");
    std::cout<<"Loading SMILES ... \n";
    while(fgets(smiles, sizeof(smiles), f)) {
        std::cout<<"\rLine: "<< ++n <<" ";
        if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) // commented to skip
//            if(strlen(smiles) > 92) // minimal query size !!!
            mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(smiles))));
    }
    fclose(f);
    printTime();
    std::cout<<"FIND MCS in "<<mols.size()<<" molecules.\n\n";
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
    printTime();
}

//====================================================================================================
//====================================================================================================

void testGregSDFFileSetFiltered() {
    const std::string sdf_dir("Greg'sComparision/data/filtered/");
    const char* sdf[] = {
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
    printf("\nTOTAL Time elapsed %.2lf seconds\n================================================\n", totalT);
}
//====================================================================================================
//====================================================================================================


int main(int argc, const char* argv[]) {
    p.Verbose = true;

// use maximum CPU resoures to increase time measuring accuracy and stability in multi process environment
#ifdef WIN32
//    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
    SetThreadPriority(GetCurrentThread (), THREAD_PRIORITY_HIGHEST );
#else
    setpriority(PRIO_PROCESS, getpid(), -20);
#endif

    T0 = nanoClock();
    t0 = nanoClock();

//testFileSDF_RandomSet("cmp_list_ChEMBL_12209_actives.dat.smi");   //FAST Test
//a lot of benchmark tests
    testChEMBL_TxtALL_chembl_II_sets();
    testFileSDF_RandomSet("chembl13-10000-random-pairs.sdf");
    testFileSDF_RandomSet("chembl13_pairs.smi");
    testFileSDF_RandomSet("chembl13_knearest_2.smi");
    testFileSDF_RandomSet("chembl13_knearest_10.smi");
    testFileSDF_RandomSet("chembl13_knearest_100.smi");
    p.Threshold = 0.8;
    testFileSDF_RandomSet("chembl13_threshold_80.smi");
    p.Threshold = 0.9;
    testFileSDF_RandomSet("chembl13_threshold_90.smi");
////    testFileSDF_RandomSet("");
    return 0;
//------------------------
    /*
    //    p.BondTyper = MCSBondCompareOrderExact;
        testChEMBL_Txt("chembl_II_sets/Target_no_10980_51302.txt"); // 271 sec  SLOW !!!
    //    return 0;

     //   testChEMBL_TxtALL_chembl_II_sets();
     //   testTarget_no_10188_30149();
    //    return 0;
        // SLOW tests
        test330();
        testChEMBL_Txt("chembl_II_sets/Target_no_10980_52937.txt");
        testChEMBL_Txt("chembl_II_sets/Target_no_11489_37339.txt");
        testChEMBL_Txt("chembl_II_sets/Target_no_10260_54285.txt");
        testChEMBL_Txt("chembl_II_sets/Target_no_10980_51302.txt"); // 271 sec  SLOW !!!
    */

#ifdef xxWIN32  // brief test set for testing and issue investigation
#ifdef _DEBUG   // check memory leaks
    _CrtMemState _ms;
    _CrtMemCheckpoint(&_ms);
#endif
    //while(1)   // check memory leaks in TaskManager or 'top -p ...'
    {
    }
#ifdef _DEBUG   // check memory leaks
    _CrtMemDumpAllObjectsSince(&_ms);
#endif
//testChEMBL_TxtALL_chembl_II_sets();
//testChEMBL_Txt("chembl_II_sets/Target_no_10980_51302.txt"); // 271 sec  SLOW !!!
    return 0;

    {
        RWMol* m = SmilesToMol("[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#7]-[#6]2:[#6]:[#6]:[#6]3:[#6](:[#6]:2):[#7]:[#7]:[#6]:3-[#6]4:[#6]:[#6]:[#6]:[#6]:[#6]:4");
#ifdef _DEBUG   // check memory leaks
        _CrtMemState _ms;
        _CrtMemCheckpoint(&_ms);
#endif
        std::cout << MolToSmarts(*m, true) <<"\n";
//    while(1)   // check memory leaks
        testTarget_no_10188_30149();
        testTarget_no_10188_49064();

#ifdef _DEBUG   // check memory leaks
        _CrtMemDumpAllObjectsSince(&_ms);
#endif
        delete m;
        return 0;

//    while(true) // MT test
//        test330();
        testChEMBLdatALL();
        testGregSDFFileSetFiltered();
        /*
            testChEMBLdat("cmp_list_ChEMBL_100126_actives.dat");
            testChEMBLdat("cmp_list_ChEMBL_100126_actives.dat", 0.8);
            testChEMBLdat("cmp_list_ChEMBL_12209_actives.dat");
        */
        return 0;

//    testFileSDF("Greg'sComparision/data/filtered/d3_aid58783.filtered.sdf");

//    testGregSDFFileSetFiltered();
        //return 0;

        return 0;

//    testGregSDFFileSetFiltered();
        //return 0;

        std::vector<unsigned> tc;
        tc.push_back(90);
        tc.push_back(326);
        tc.push_back(330);
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
        testFileMCSB(argv[1], 300, tc);
        return 0;
    }
#endif

    if(3 == argc && '-' == argv[1][0])
        switch(argv[1][1]) { // ./test -s|m|b <filename with test files list>
        case 's': { // smiles files list
            char test[256];
            FILE* f = fopen(argv[2], "rt");
            while(fgets(test, sizeof(test), f))
                testFileSMILES(test);
            fclose(f);
        }
        break;
        case 'm': { // SDF mol files list
            char test[256];
            FILE* f = fopen(argv[2], "rt");
            while(fgets(test, sizeof(test), f))
                testFileSDF(test);
            fclose(f);
        }
        break;
        case 'b': {
            std::vector<unsigned> tc; // empty -> all
            testFileMCSB(argv[2], 30, tc);    // .mcsb
        }
        break;
        default:
            break;
        }
    else if(2 == argc) { // .sdf /.smi file
        if(0==strcmp(argv[1]+strlen(argv[1])-4, ".smi"))
            testFileSMILES(argv[1]);// .smi
        else if(0==strcmp(argv[1]+strlen(argv[1])-4, ".sdf"))
            testFileSDF(argv[1]);   // .sdf
        else if(0==strcmp(argv[1]+strlen(argv[1])-4, "mcsb"))
            testFileMCSB(argv[1], 30);  // .mcsb
        else if(0==strcmp(argv[1]+strlen(argv[1])-4, ".dat"))
            testChEMBLdat(argv[1]);   // .sdf
        else
            printf("UNKNOWN File Extention.\n");
    } else if(argc > 1+2)
        testCmndLineSMILES(argc, argv);
    else {
        testGregSDFFileSetFiltered();
    }
//  BOOST_LOG(rdInfoLog) << "*******************************************************\n";

    unsigned long long t1 = nanoClock();
    double sec = double(t1-T0) / 1000000.;
    printf("TOTAL Time elapsed %.2lf seconds\n", sec);

    return 0;
}


