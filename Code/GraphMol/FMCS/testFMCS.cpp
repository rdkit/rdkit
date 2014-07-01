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
//#include <RDGeneral/RDLog.h>
//#include <RDGeneral/utils.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h" //MOL single molecule !
#include "../FileParsers/MolSupplier.h" //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "FMCS.h"
#include "DebugTrace.h" //#ifdef VERBOSE_STATISTICS_ON

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

void printTime()
{
  unsigned long long t1 = nanoClock();
  double sec = double(t1-t0) / 1000000.;
  printf("\nTime elapsed %.3lf seconds\n", sec);
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

std::string getSmilesOnlyTxt(const char* smiles, std::string* id=0) // remove label from string like "CHEMBL90218 NS(=O)(=O)c1ccc(NC(=O)c2cccc(C(=O)O)n2)c(Cl)c1"
{
  const char* sp = strchr(smiles,' ');
  if(sp && '\0'!=*sp)
    return std::string(++sp);
  else
    return smiles;
}


std::string getSmilesOnlyChEMBL(const char* smiles, std::string* id=0) // remove label, because RDKit parse FAILED
{
  const char* sp = strchr(smiles, '\t');
  if(sp)
    {
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


void testFileMCSB(const char* test, unsigned timeout=30, std::vector<unsigned> test_N=std::vector<unsigned>())  // optional list of some tests for investigation
{
  p.Verbose = false;

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
      mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i], &id) )));
      std::cout << id << "\n";
#endif
    }

  t0 = nanoClock();

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
  RWMol* qm = SmilesToMol( getSmilesOnly(smi[0]) );
  unsigned nq = qm->getNumAtoms();
  for(size_t ai = 0; ai < nq; ai++)
    {
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
  RWMol* qm = SmilesToMol( getSmilesOnly(smi[0]) );
  unsigned nq = qm->getNumAtoms();
  for(size_t ai = 0; ai < nq; ai++)
    {
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
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  findMCS(mols);
  p.Threshold = 0.7;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}

void test330()
{
  std::cout << "\ntest330()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
      //TEST 330  40 sec
      "CCC(C)C(NC(=O)C(NC(C(CCC(O)=O)NC(=O)C(NC(=O)C(NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(CN)=O)=O)=O)=O)C(C)CC)C(C)CC)=O)CCCCN)C(NC(C)C(NC(CCCCN)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O)=O)=O CHEMBL1240742",
      "CCC(C)C(NC(=O)C(NC(C(CCCCN)NC(=O)C(NC(=O)C(NC(C(CC(O)=O)NC(C(Cc1ccccc1)NC(C(CC(C)C)NC(CN)=O)=O)=O)=O)C(C)CC)C(C)CC)=O)CCCCN)C(NC(C)C(NC(CCC(O)=O)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O)=O)=O CHEMBL1240736",
      "CCC(C)C(NC(CN)=O)C(NC(C(NC(CC(O)=O)C(NC(C(NC(C)C(NC(CCCCN)C(NC(C(NC(CC(C)C)C(NC(Cc1ccccc1)C(NC(CCC(O)=O)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O)=O)=O)=O)CCCCN)=O)=O)=O)C(C)CC)=O)=O)C(C)CC)=O CHEMBL1240738",
      "CCC(C)C(NC(CN)=O)C(NC(Cc1ccccc1)C(NC(CC(O)=O)C(NC(CCCCN)C(NC(CC(C)C)C(NC(C)C(NC(CCCCN)C(NC(CCC(O)=O)C(NC(C(NC(CO)C(NC(C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)C(C)CC)=O)=O)C(C)CC)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL1240740",
      "CCC(C)C(NC(CN)=O)C(NC(Cc1c[nH]c2c1cccc2)C(NC(CO)C(NC(CC(O)=O)C(NC(CC(C)C)C(NC(C)C(NC(CCC(O)=O)C(NC(C(NC(C(NC(CCCCN)C(NC(CCCCN)C(NC(Cc1ccccc1)C(O)=O)=O)=O)=O)C(C)CC)=O)C(C)CC)=O)=O)=O)=O)=O)=O)=O CHEMBL1240741",
      "CCC(C)C(NC(=O)C(NC(=O)C(CCCCN)NC(C(CC(C)C)NC(C(Cc1c[nH]c2c1cccc2)NC(CN)=O)=O)=O)CCCCN)C(NC(CCC(O)=O)C(NC(CO)C(=O)NC(C(NC(C(NC(CC(O)=O)C(NC(C)C(NC(Cc1ccccc1)C(O)=O)=O)=O)=O)C(C)CC)=O)C(C)CC)=O)=O CHEMBL1240743",
      "CCC(C)C(NC(C(CCC(O)=O)NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)=O)C(NC(Cc1c[nH]c2ccccc12)C(O)=O)=O CHEMBL431874",
      "CCC(C)C(NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)C(NC(CCC(O)=O)C(NC(Cc1c[nH]c2ccccc12)C(O)=O)=O)=O CHEMBL262166",
      "CCC(C)C(NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)C(NC(CCCCN)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O CHEMBL313122",
      "CCC(C)C(NC(C(CCCCN)NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)=O)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O CHEMBL314239",
      //# 330 F  42 41 30.93 sec MCS: [#6]-[#6](-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6](-[#6]-[#6]-[#6])-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6])=[#8])=[#8])=[#8])=[#8])=[#8])-[#6](-[#7]-[#6](-[#6]-[#6](:[#6]):[#6]:[#6]:[#6]:[#6])-[#6](-[#8])=[#8])=[#8]
    };
  for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}


void testChEMBL_Txt(const char* test, double th=1.0, const char* csv="chembl_II_sets.C++.res.csv")
{
  std::cout << "\ntestChEMBL_Txt() "<<test<<"\n";
  std::vector<ROMOL_SPTR> mols;
  char smiles[4096];
  unsigned n=0;
  FILE* f = fopen(test, "rt");
  if(!f)
    {
      perror("fopen testChEMBL_Txt()");
      return;
    }
  FILE* fres = fopen(csv, "at");
  while(fgets(smiles, sizeof(smiles), f))
    {
      if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) // commented to skip
        {
          mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnlyTxt(smiles))));
        }
    }
  fclose(f);
  t0 = nanoClock();
  p.Threshold = th;
  MCSResult res = findMCS(mols, &p);
  unsigned long long tc1 = nanoClock();
  double sec = double(tc1-t0) / 1000000.;
  std::cout << "MCS : "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
  fprintf(fres, "%s; %.1f; %.2f; %u; %u; %s\n", test, th, sec, res.NumAtoms, res.NumBonds, res.SmartsString.c_str());
  fclose(fres);
}

void testChEMBL_TxtALL_chembl_II_sets(double th=1.0)
{
  FILE* fres = fopen("chembl_II_sets.C++.res.csv", "at");
  fprintf(fres,"test;threshold;t,sec;atoms;bonds;mcs\n");
  fclose(fres);
  const char* test[] =
    {
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
  for(int i=0; i<sizeof(test)/sizeof(test[0]); i++)
    testChEMBL_Txt((std::string("chembl_II_sets/") + test[i]).c_str(), th, "chembl_II_sets.C++.res.csv");
}

void testChEMBL_TxtSLOW_chembl_II_sets(double th=1.0)
{
  FILE* fres = fopen("chembl_II_sets.SLOW.C++.res.csv", "wt");
  fprintf(fres,"test;threshold;t,sec;atoms;bonds;mcs\n");
  fclose(fres);
  const char* test[] =
    {
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

void testChEMBLdat(const char* test, double th=1.0)
{
  std::cout << "\ntestChEMBLdat() "<<test<<"\n";
  std::vector<ROMOL_SPTR> mols;
  char smiles[4096];
  unsigned n=0;
  FILE* f = fopen(test, "rt");
  if(!f)
    {
      perror("fopen testChEMBLdat()");
      return;
    }
  FILE* fs = fopen((std::string(test)+".smi").c_str(), "wt");
  std::cout<<"Loading SMILES ... \n";
  t0 = nanoClock();
  while(fgets(smiles, sizeof(smiles), f))
    {
      std::cout<<"\rLine: "<< ++n <<" ";
      if('#' != smiles[0] && ' ' != smiles[0] && '/' != smiles[0]) // commented to skip
        {
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

void testChEMBLdatALL(double th=1.0)
{
  const char* test[] =
    {
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


void testTarget_no_10188_30149()
{
  std::cout << "\ntestTarget_no_10188_30149()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
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


void testTarget_no_10188_49064()
{
  std::cout << "\ntestTarget_no_10188_49064()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
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
  p.Threshold = 0.7;
  p.BondCompareParameters.RingMatchesRingOnly = true;
  p.BondCompareParameters.CompleteRingsOnly   = true;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}

void testSegFault()
{
  std::cout << "\ntestSegFault()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
      "CN(CCCN(C)CCc1ccccc1)CCOC(c1ccccc1)c1ccccc1",
      "CN(CCCc1ccccc1)CCCN(C)CCOC(c1ccccc1)c1ccccc1",
      "Fc1ccc(C(OCCNCCCNCCc2ccccc2)c2ccc(F)cc2)cc1",
      "O=C(Cc1ccccc1)NCCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(Br)cc1)NC=CNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(F)cc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccc(F)cc1)c1ccc(F)cc1)CCN(C)CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "COC(=O)C1C2CCC(CC1C(=O)Oc1ccccc1)N2C",
      "O=C1CN(CCc2ccccc2)CCN1CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1",
    };
  for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}


void testAtomCompareIsotopes()
{
  std::cout << "\ntestAtomCompareIsotopes()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
      "CC[13NH2]",
      "CC[13CH3]",
    };
  for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  t0 = nanoClock();
  p.AtomTyper = MCSAtomCompareIsotopes;
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}


void testAtomCompareAnyAtom()
{
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
      "c1ccccc1C",
      "c1ccccc1O",
      "c1ccccc1Cl",
      "c1ccccc1F",
      "c1ccccc1N",
    };
  for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  t0 = nanoClock();
  p.AtomTyper = MCSAtomCompareAny;
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: "<<res.SmartsString<<" "<< res.NumAtoms<<" atoms, "<<res.NumBonds<<" bonds\n";
  printTime();
}

void testAtomCompareAnyAtomBond()
{
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] =
    {
      "C1CCCCC1=C",
      "c1ccccc1O",
      "c1ccccc1Cl",
      "c1ccccc1F",
      "c1ccccc1N",
    };
  for(int i=0; i<sizeof(smi)/sizeof(smi[0]); i++)
    mols.push_back(ROMOL_SPTR(SmilesToMol( getSmilesOnly(smi[i]) )));
  t0 = nanoClock();
  p.AtomTyper = MCSAtomCompareAny;
  p.BondTyper = MCSBondCompareAny;
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
    mols.push_back(ROMOL_SPTR(SmilesToMol(argv[i])));
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

void testGregSDFFileSetFiltered()
{
  const std::string sdf_dir("Greg'sComparision/data/filtered/");
  const char* sdf[] =
    {
      //fmcs:     0.11            27      29      beta2_adrenergic_aid465635.filtered.sdf O=C(O)C1:C:C:C(C2:C:C:C(CCNCC(O)C3:C:C:C:C:C:3):C:C:2):C:C:1
      "beta2_adrenergic_aid465635.filtered.sdf",
      //fmcs:     1.90            22      23      beta2_adrenergic_aid578239.filtered.sdf CNC1:C:C:C(CCNCC(O)COC2:C:C:C:C:C:2):C:C:1
      "beta2_adrenergic_aid578239.filtered.sdf",
      //fmcs:     1.88            22      23      beta2_adrenergic_aid578240.filtered.sdf CNC1:C:C:C(CCNCC(O)COC2:C:C:C:C:C:2):C:C:1
      "beta2_adrenergic_aid578240.filtered.sdf",
      //fmcs:     0.07            17      18      beta2_adrenergic_aid759384.filtered.sdf CCNCC(O)C1:C:C:C(O):C2:N:C(=O):[SH]:C:2:1
      "beta2_adrenergic_aid759384.filtered.sdf",
      //fmcs:     1.01            25      25      d3_aid329485.filtered.sdf       C:C:C:C:C:CC(=O)NCCCCN1CCN(C:C:C:C:C:C)CC1
      "d3_aid329485.filtered.sdf",
      //fmcs:     0.13            20      21      d3_aid367038.filtered.sdf       C:C:C:C:N:C:CCN1CCN(C2:C:C:C:C:C:2)CC1
      "d3_aid367038.filtered.sdf",
      //fmcs:     0.78            27      29      d3_aid563533.filtered.sdf       O=C(C:C:C1:C:C:C:C:C:1)NCCCCN1CCN(C2:C:C:C:C:C:2)CC1
      "d3_aid563533.filtered.sdf",
      //fmcs:     0.09            14      14      d3_aid563770.filtered.sdf       CCCN(CC)C(=O)C1:C:C:C:C:C:1
      "d3_aid563770.filtered.sdf",
      //fmcs:     0.14            22      24      d3_aid578980.filtered.sdf       C(:C:N:C1:C:C:N:C:N:1)CN1CCN(C2:C:C:C:C:C:2)CC1
      "d3_aid578980.filtered.sdf",
      //fmcs:     0.20            15      17      d3_aid58783.filtered.sdf        CCN1CC2COC3:C:C:C:C:C:3C2C1
      "d3_aid58783.filtered.sdf",
      //fmcs:     0.08            14      14      d3_aid62278.filtered.sdf        C:C(CNCC):N:CC1:C:C:C:C:C:1
      "d3_aid62278.filtered.sdf",
      //fmcs:     0.05            7       7       d3_aid62281.filtered.sdf        CC1:C:C:C:C:C:1
      "d3_aid62281.filtered.sdf",
      //fmcs:     33.13           26      27      d3_aid62457.filtered.sdf        CC(=O)NCCCCN1CCC2:C:C:C(OS(=O)(=O)C(F)(F)F):C:C:2C1
      "d3_aid62457.filtered.sdf",
      //fmcs:     0.08            14      15      d3_aid62774.filtered.sdf        CCN1CCCC(C2:C:C:C:C:C:2)C1
      "d3_aid62774.filtered.sdf",
      //fmcs:     0.33            14      14      d3_aid642590.filtered.sdf       C:C:C:C:C:C:CCN1CCCCC1
      "d3_aid642590.filtered.sdf",
    };

  double totalT=0.;
  for(int i=0; i<sizeof(sdf)/sizeof(sdf[0]); i++)
    totalT += testFileSDF((sdf_dir+sdf[i]).c_str());
  printf("\nTOTAL Time elapsed %.2lf seconds\n================================================\n", totalT);
}
//====================================================================================================
//====================================================================================================


int main(int argc, const char* argv[])
{
  //    p.Verbose = true;

  // use maximum CPU resoures to increase time measuring accuracy and stability in multi process environment
#ifdef WIN32
  //    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
  SetThreadPriority(GetCurrentThread (), THREAD_PRIORITY_HIGHEST );
#else
  setpriority(PRIO_PROCESS, getpid(), -20); 
#endif

  T0 = nanoClock();
  t0 = nanoClock();

  testSimpleFast();
  //testChEMBL_TxtSLOW_chembl_II_sets();    // python total time is about 55/2 sec
  return 0;
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
  return 0;
  */
#ifdef xxWIN32  // brief test set for testing and issue investigation
#ifdef _DEBUG   // check memory leaks
  _CrtMemState _ms;
  _CrtMemCheckpoint(&_ms);
#endif
  //while(1)   // check memory leaks in TaskManager or 'top -p ...'
  {
    testTarget_no_10188_30149();
    testAtomCompareAnyAtom();
    testAtomCompareAnyAtomBond();
  }
#ifdef _DEBUG   // check memory leaks
  _CrtMemDumpAllObjectsSince(&_ms);
#endif
  //testAtomCompareIsotopes();
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

    testSimple();
    //while(true) // MT test
    //    test330();
    //testChEMBLdatALL();
    //testGregSDFFileSetFiltered();
    /*
      testChEMBLdat("cmp_list_ChEMBL_100126_actives.dat");
      testChEMBLdat("cmp_list_ChEMBL_100126_actives.dat", 0.8);
      testChEMBLdat("cmp_list_ChEMBL_12209_actives.dat"); 
    */
    return 0;

    //    testSimpleFast();
    testSimple();
    //    testFileSDF("Greg'sComparision/data/filtered/d3_aid58783.filtered.sdf");
    //    testSimple();

    //    testGregSDFFileSetFiltered();
    //test18();
    //return 0;

    //testThreshold();
    //    testRing1();
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
      else if(0==strcmp(argv[1]+strlen(argv[1])-4, ".dat"))
        testChEMBLdat(argv[1]);   // .sdf
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
  printf("TOTAL Time elapsed %.2lf seconds\n", sec);

  return 0;
}


